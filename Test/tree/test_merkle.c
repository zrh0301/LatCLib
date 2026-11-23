#include "tree/merkle.h"
#include "sym/sm3.h"
#include "util/bytes.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "line %d: %s\n", __LINE__, #condition); \
        return 1; \
    } \
} while (0)

enum { BYTES = 32, MAX_LEAVES = 16, MAX_HEIGHT = 4 };

typedef struct {
    uint8_t values[MAX_HEIGHT + 1][MAX_LEAVES][BYTES];
    pqc_merkle_node proof[MAX_LEAVES * MAX_HEIGHT];
    size_t proof_count;
    uint32_t deep_leaves;
    uint32_t selected_mask;
    uint32_t offset;
} reference_tree;

static void make_leaf(uint8_t *out, uint32_t index, void *ctx)
{
    (void)ctx;
    uint8_t encoded[4];
    pqc_store_be(encoded, sizeof(encoded), index);
    sm3_hash(out, encoded, sizeof(encoded));
}

static void hash_pair(uint8_t *out, const uint8_t *children,
                       uint32_t height, uint32_t index, void *ctx)
{
    (void)ctx;
    uint8_t input[8 + 2 * BYTES];
    pqc_store_be(input, 4, height);
    pqc_store_be(input + 4, 4, index);
    memcpy(input + 8, children, 2 * BYTES);
    sm3_hash(out, input, sizeof(input));
}

/* Recursive construction visits actual left and right children. Each return
 * reports whether that subtree contains any selected leaf. */
static int reference_node(reference_tree *tree, uint32_t height, uint32_t index)
{
    if (height == 0 || (height == 1 && index >= tree->deep_leaves / 2)) {
        const uint32_t leaf = height == 0 ? index : tree->deep_leaves / 2 + index;
        make_leaf(tree->values[height][index], leaf + tree->offset, NULL);
        return (tree->selected_mask >> leaf) & 1u;
    }
    const int left = reference_node(tree, height - 1, 2 * index);
    const int right = reference_node(tree, height - 1, 2 * index + 1);
    uint8_t children[2 * BYTES];
    memcpy(children, tree->values[height - 1][2 * index], BYTES);
    memcpy(children + BYTES, tree->values[height - 1][2 * index + 1], BYTES);
    hash_pair(tree->values[height][index], children, height,
               index + (tree->offset >> height), NULL);
    if (left != right)
        tree->proof[tree->proof_count++] =
            (pqc_merkle_node){height - 1, 2 * index + (left ? 1u : 0u)};
    return left || right;
}

static int compare_positions(const void *a, const void *b)
{
    const pqc_merkle_node *left = a, *right = b;
    if (left->height != right->height)
        return left->height < right->height ? -1 : 1;
    return (left->index > right->index) - (left->index < right->index);
}

static uint32_t tree_height(uint32_t leaves)
{
    uint32_t height = 0, width = 1;
    while (width < leaves) {
        width *= 2;
        ++height;
    }
    return height;
}

static int test_multiproofs(uint32_t leaf_count)
{
    const uint32_t height = tree_height(leaf_count);
    const uint32_t deep = height == 0 ? 1 : 2 * (leaf_count - (1u << (height - 1)));
    for (uint32_t mask = 1; mask < (1u << leaf_count); ++mask) {
        reference_tree expected = {.deep_leaves = deep, .selected_mask = mask};
        reference_node(&expected, height, 0);
        qsort(expected.proof, expected.proof_count, sizeof(*expected.proof), compare_positions);
        uint32_t selected[MAX_LEAVES];
        size_t count = 0;
        uint8_t leaves[MAX_LEAVES * BYTES];
        for (uint32_t i = 0; i < leaf_count; ++i) {
            if ((mask >> i) & 1u) {
                selected[count] = i;
                make_leaf(leaves + count++ * BYTES, i, NULL);
            }
        }

        pqc_merkle_node actual[MAX_LEAVES * MAX_HEIGHT];
        uint8_t auth[MAX_LEAVES * MAX_HEIGHT * BYTES];
        uint8_t root[BYTES], rebuilt[BYTES];
        size_t auth_count = SIZE_MAX;
        CHECK(pqc_merkle_auth_nodes(actual, &auth_count, MAX_LEAVES * MAX_HEIGHT,
                                    selected, count, leaf_count) == 0);
        CHECK(auth_count == expected.proof_count);
        CHECK(memcmp(actual, expected.proof, auth_count * sizeof(*actual)) == 0);
        CHECK(pqc_merkle_build(root, auth, BYTES, leaf_count, actual, auth_count,
                               make_leaf, NULL, hash_pair, NULL) == 0);
        CHECK(memcmp(root, expected.values[height][0], BYTES) == 0);
        for (size_t i = 0; i < auth_count; ++i)
            CHECK(memcmp(auth + i * BYTES,
                         expected.values[actual[i].height][actual[i].index], BYTES) == 0);

        CHECK(pqc_merkle_root_from_proof(rebuilt, selected, leaves, count, auth,
                                         auth_count, BYTES, leaf_count, hash_pair, NULL) == 0);
        CHECK(memcmp(root, rebuilt, BYTES) == 0);
        if (auth_count != 0) {
            memset(rebuilt, 0x5a, sizeof(rebuilt));
            CHECK(pqc_merkle_root_from_proof(rebuilt, selected, leaves, count, auth,
                                             auth_count - 1, BYTES, leaf_count,
                                             hash_pair, NULL) == -1);
            for (size_t i = 0; i < BYTES; ++i)
                CHECK(rebuilt[i] == 0x5a);
            auth[0] ^= 1;
            CHECK(pqc_merkle_root_from_proof(rebuilt, selected, leaves, count, auth,
                                             auth_count, BYTES, leaf_count,
                                             hash_pair, NULL) == 0);
            CHECK(memcmp(root, rebuilt, BYTES) != 0);
        }
    }
    return 0;
}

static int test_complete(uint32_t height, uint32_t offset)
{
    const uint32_t leaves = 1u << height;
    reference_tree expected = {.deep_leaves = leaves, .offset = offset};
    reference_node(&expected, height, 0);
    uint8_t root[BYTES], path[MAX_HEIGHT * BYTES], leaf[BYTES];
    for (uint32_t selected = 0; selected < leaves; ++selected) {
        CHECK(pqc_merkle_treehash(root, path, BYTES, height, selected, offset,
                                  make_leaf, NULL, hash_pair, NULL) == 0);
        CHECK(memcmp(root, expected.values[height][0], BYTES) == 0);
        for (uint32_t h = 0; h < height; ++h)
            CHECK(memcmp(path + h * BYTES,
                         expected.values[h][(selected >> h) ^ 1u], BYTES) == 0);
        make_leaf(leaf, selected + offset, NULL);
        CHECK(pqc_merkle_root_from_path(leaf, leaf, path, BYTES, height, selected,
                                        offset, hash_pair, NULL) == 0);
        CHECK(memcmp(root, leaf, BYTES) == 0);
    }
    CHECK(pqc_merkle_treehash(root, NULL, BYTES, height, UINT32_MAX, offset,
                              make_leaf, NULL, hash_pair, NULL) == 0);
    CHECK(memcmp(root, expected.values[height][0], BYTES) == 0);
    return 0;
}

static int test_invalid(void)
{
    uint8_t root[BYTES], auth[BYTES];
    memset(root, 0x5a, sizeof(root));
    const uint32_t duplicate[] = {1, 1};
    const uint32_t out_of_range[] = {2};
    const uint32_t selected[] = {0};
    pqc_merkle_node position = {2, 0};
    size_t count = 7;
    CHECK(pqc_merkle_auth_nodes(&position, &count, 2, duplicate, 2, 2) == -1);
    CHECK(pqc_merkle_auth_nodes(&position, &count, 2, out_of_range, 1, 2) == -1);
    CHECK(pqc_merkle_auth_nodes(NULL, &count, 0, selected, 1, 2) == -1);
    CHECK(count == 7);
    pqc_merkle_node wide[32];
    CHECK(pqc_merkle_auth_nodes(wide, &count, 32, selected, 1, UINT32_MAX) == 0);
    CHECK(count == 32 && wide[0].height == 0 && wide[31].height == 31);
    const uint32_t last_leaf[] = {UINT32_MAX - 1};
    CHECK(pqc_merkle_auth_nodes(wide, &count, 32, last_leaf, 1, UINT32_MAX) == 0);
    CHECK(count == 31 && wide[0].height == 1 && wide[30].height == 31);
    CHECK(pqc_merkle_build(root, auth, BYTES, 2, &position, 1,
                           make_leaf, NULL, hash_pair, NULL) == -1);
    CHECK(pqc_merkle_treehash(root, auth, BYTES, 32, 0, 0,
                              make_leaf, NULL, hash_pair, NULL) == -1);
    CHECK(pqc_merkle_treehash(root, auth, BYTES, 1, 0, 1,
                              make_leaf, NULL, hash_pair, NULL) == -1);
    CHECK(pqc_merkle_treehash(root, auth, SIZE_MAX, 0, 0, 0,
                              make_leaf, NULL, hash_pair, NULL) == -1);
    CHECK(pqc_merkle_root_from_path(root, root, NULL, BYTES, 1, 0, 0,
                                    hash_pair, NULL) == -1);
    CHECK(pqc_merkle_root_from_proof(root, selected, root, 1, auth, 1, BYTES, 1,
                                     hash_pair, NULL) == -1);
    for (size_t i = 0; i < sizeof(root); ++i)
        CHECK(root[i] == 0x5a);
    return 0;
}

int main(void)
{
    for (uint32_t height = 0; height <= MAX_HEIGHT; ++height) {
        CHECK(test_complete(height, 0) == 0);
        CHECK(test_complete(height, 3u << height) == 0);
    }
    for (uint32_t leaves = 1; leaves <= 9; ++leaves)
        CHECK(test_multiproofs(leaves) == 0);
    CHECK(test_invalid() == 0);
    return 0;
}
