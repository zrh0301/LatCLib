#include "tree/merkle.h"

#include <stdlib.h>
#include <string.h>

static int complete_dimensions(size_t bytes, uint32_t height, uint32_t leaf,
                                uint32_t offset)
{
    if (bytes == 0 || bytes > SIZE_MAX / 35 || height >= 32)
        return 0;
    const uint32_t width = UINT32_C(1) << height;
    return leaf < width && offset % width == 0 && offset <= UINT32_MAX - (width - 1);
}

int pqc_merkle_treehash(uint8_t *root, uint8_t *auth_path, size_t bytes,
                        uint32_t height, uint32_t leaf_index, uint32_t offset,
                        pqc_merkle_leaf_fn leaf, void *leaf_ctx,
                        pqc_merkle_hash_fn hash, void *hash_ctx)
{
    const uint32_t selected = leaf_index == UINT32_MAX ? 0 : leaf_index;
    if (!root || !leaf || !hash ||
        !complete_dimensions(bytes, height, selected, offset))
        return -1;
    uint8_t *stack = malloc((height + 3u) * bytes);
    if (!stack)
        return -1;
    uint8_t *children = stack + height * bytes;
    uint8_t *node = children + bytes;
    const uint32_t last = (UINT32_C(1) << height) - 1;

    for (uint32_t index = 0;; ++index) {
        leaf(node, index + offset, leaf_ctx);
        uint32_t position = index;
        uint32_t selected_position = leaf_index;
        uint32_t parent_offset = offset;
        for (uint32_t level = 0;; ++level) {
            if (level == height) {
                memcpy(root, node, bytes);
                free(stack);
                return 0;
            }
            if (auth_path && leaf_index != UINT32_MAX &&
                (position ^ selected_position) == 1)
                memcpy(auth_path + level * bytes, node, bytes);
            if ((position & 1u) == 0) {
                memcpy(stack + level * bytes, node, bytes);
                break;
            }
            memcpy(children, stack + level * bytes, bytes);
            position >>= 1;
            selected_position >>= 1;
            parent_offset >>= 1;
            hash(node, children, level + 1, position + parent_offset, hash_ctx);
        }
        if (index == last)
            break;
    }
    free(stack);
    return -1;
}

int pqc_merkle_root_from_path(uint8_t *root, const uint8_t *leaf,
                             const uint8_t *auth_path, size_t bytes,
                             uint32_t height, uint32_t leaf_index, uint32_t offset,
                             pqc_merkle_hash_fn hash, void *hash_ctx)
{
    if (!root || !leaf || !hash || (height != 0 && !auth_path) ||
        !complete_dimensions(bytes, height, leaf_index, offset))
        return -1;
    uint8_t *children = malloc(3 * bytes);
    if (!children)
        return -1;
    uint8_t *node = children + 2 * bytes;
    memcpy(node, leaf, bytes);
    for (uint32_t level = 0; level < height; ++level) {
        const unsigned side = leaf_index & 1u;
        memcpy(children + side * bytes, node, bytes);
        memcpy(children + (1u - side) * bytes, auth_path + level * bytes, bytes);
        leaf_index >>= 1;
        offset >>= 1;
        hash(node, children, level + 1, leaf_index + offset, hash_ctx);
    }
    memcpy(root, node, bytes);
    free(children);
    return 0;
}

typedef struct {
    uint32_t height;
    uint32_t deep_leaves;
    uint32_t base;
} tree_shape;

static tree_shape shape_for(uint32_t leaves)
{
    tree_shape shape = {0, 1, 1};
    uint64_t width = 1;
    while (width < leaves) {
        width <<= 1;
        ++shape.height;
    }
    if (shape.height != 0) {
        shape.base = (uint32_t)(width / 2);
        shape.deep_leaves = 2u * (leaves - shape.base);
    }
    return shape;
}

static uint32_t level_width(const tree_shape *shape, uint32_t height)
{
    if (height == 0)
        return shape->deep_leaves;
    return shape->base >> (height - 1);
}

static int valid_selection(const uint32_t *selected, size_t count, uint32_t leaves)
{
    if (!selected || count == 0 || count > leaves)
        return 0;
    for (size_t i = 0; i < count; ++i)
        if (selected[i] >= leaves || (i != 0 && selected[i] <= selected[i - 1]))
            return 0;
    return 1;
}

static size_t merge_indices(uint32_t *out, const uint32_t *left, size_t left_count,
                             const uint32_t *right, size_t right_count)
{
    size_t n = 0, i = 0, j = 0;
    while (i < left_count || j < right_count) {
        if (j == right_count || (i < left_count && left[i] < right[j]))
            out[n++] = left[i++];
        else
            out[n++] = right[j++];
    }
    return n;
}

int pqc_merkle_auth_nodes(pqc_merkle_node *out, size_t *auth_count,
                          size_t capacity, const uint32_t *selected,
                          size_t count, uint32_t leaves)
{
    if (!auth_count || !valid_selection(selected, count, leaves) ||
        count > SIZE_MAX / (4 * sizeof(uint32_t)))
        return -1;
    const tree_shape shape = shape_for(leaves);
    if ((shape.height != 0 && (!out || count > capacity / shape.height)))
        return -1;
    uint32_t *work = malloc(4 * count * sizeof(*work));
    if (!work)
        return -1;
    uint32_t *current = work, *upper = work + count;
    uint32_t *parents = work + 2 * count, *merged = work + 3 * count;
    size_t current_count = 0, upper_count = 0, total = 0;

    for (size_t i = 0; i < count; ++i) {
        if (selected[i] < shape.deep_leaves)
            current[current_count++] = selected[i];
        else
            upper[upper_count++] = selected[i] - shape.deep_leaves / 2;
    }
    for (uint32_t height = 0; height < shape.height; ++height) {
        size_t parent_count = 0;
        for (size_t i = 0; i < current_count;) {
            const uint32_t node = current[i];
            parents[parent_count++] = node >> 1;
            if (i + 1 < current_count && current[i + 1] == (node ^ 1u)) {
                i += 2;
            } else {
                out[total++] = (pqc_merkle_node){height, node ^ 1u};
                ++i;
            }
        }
        current_count = merge_indices(merged, parents, parent_count, upper, upper_count);
        memcpy(current, merged, current_count * sizeof(*current));
        upper_count = 0;
    }
    free(work);
    *auth_count = total;
    return 0;
}

int pqc_merkle_build(uint8_t *root, uint8_t *auth_values, size_t bytes,
                     uint32_t leaves, const pqc_merkle_node *auth,
                     size_t auth_count, pqc_merkle_leaf_fn leaf, void *leaf_ctx,
                     pqc_merkle_hash_fn hash, void *hash_ctx)
{
    if (!root || leaves == 0 || bytes == 0 || !leaf || !hash ||
        auth_count > SIZE_MAX / bytes ||
        (auth_count != 0 && (!auth || !auth_values)))
        return -1;
    const tree_shape shape = shape_for(leaves);
    size_t offsets[33], total = 0;
    for (uint32_t h = 0; h <= shape.height; ++h) {
        const size_t width = level_width(&shape, h);
        if (width > (SIZE_MAX - total) / bytes)
            return -1;
        offsets[h] = total;
        total += width * bytes;
    }
    for (size_t i = 0; i < auth_count; ++i)
        if (auth[i].height > shape.height ||
            auth[i].index >= level_width(&shape, auth[i].height))
            return -1;
    uint8_t *tree = malloc(total);
    if (!tree)
        return -1;

    for (uint32_t i = 0; i < shape.deep_leaves; ++i)
        leaf(tree + (size_t)i * bytes, i, leaf_ctx);
    if (shape.height != 0) {
        for (uint32_t i = 0; i < shape.base; ++i) {
            uint8_t *node = tree + offsets[1] + (size_t)i * bytes;
            if (i < shape.deep_leaves / 2)
                hash(node, tree + (size_t)(2u * i) * bytes, 1, i, hash_ctx);
            else
                leaf(node, shape.deep_leaves / 2 + i, leaf_ctx);
        }
    }
    for (uint32_t h = 2; h <= shape.height; ++h) {
        for (uint32_t i = 0; i < level_width(&shape, h); ++i)
            hash(tree + offsets[h] + (size_t)i * bytes,
                 tree + offsets[h - 1] + (size_t)(2u * i) * bytes, h, i, hash_ctx);
    }
    for (size_t i = 0; i < auth_count; ++i)
        memcpy(auth_values + i * bytes,
               tree + offsets[auth[i].height] + (size_t)auth[i].index * bytes, bytes);
    memcpy(root, tree + offsets[shape.height], bytes);
    free(tree);
    return 0;
}

int pqc_merkle_root_from_proof(uint8_t *root, const uint32_t *selected,
                              const uint8_t *leaves, size_t count,
                              const uint8_t *auth_values, size_t auth_count,
                              size_t bytes, uint32_t leaf_count,
                              pqc_merkle_hash_fn hash, void *hash_ctx)
{
    if (!root || !leaves || !hash || bytes == 0 ||
        auth_count > SIZE_MAX / bytes ||
        (auth_count != 0 && !auth_values) ||
        !valid_selection(selected, count, leaf_count) ||
        count > SIZE_MAX / (4 * sizeof(uint32_t)) ||
        bytes > SIZE_MAX / 4 || count > (SIZE_MAX - 2 * bytes) / (4 * bytes))
        return -1;
    const tree_shape shape = shape_for(leaf_count);
    uint32_t *positions = malloc(4 * count * sizeof(*positions));
    uint8_t *values = malloc((4 * count + 2) * bytes);
    if (!positions || !values) {
        free(positions);
        free(values);
        return -1;
    }

    uint32_t *current = positions, *upper = positions + count;
    uint32_t *parents = positions + 2 * count, *merged = positions + 3 * count;
    uint8_t *current_values = values, *upper_values = values + count * bytes;
    uint8_t *parent_values = values + 2 * count * bytes;
    uint8_t *merged_values = values + 3 * count * bytes;
    uint8_t *children = values + 4 * count * bytes;
    size_t current_count = 0, upper_count = 0, cursor = 0;
    int result = -1;
    for (size_t i = 0; i < count; ++i) {
        if (selected[i] < shape.deep_leaves) {
            current[current_count] = selected[i];
            memcpy(current_values + current_count++ * bytes, leaves + i * bytes, bytes);
        } else {
            upper[upper_count] = selected[i] - shape.deep_leaves / 2;
            memcpy(upper_values + upper_count++ * bytes, leaves + i * bytes, bytes);
        }
    }
    for (uint32_t height = 0; height < shape.height; ++height) {
        size_t parent_count = 0;
        for (size_t i = 0; i < current_count;) {
            const uint32_t node = current[i];
            if (i + 1 < current_count && current[i + 1] == (node ^ 1u)) {
                memcpy(children, current_values + i * bytes, 2 * bytes);
                i += 2;
            } else {
                if (cursor >= auth_count)
                    goto done;
                const unsigned side = node & 1u;
                memcpy(children + side * bytes, current_values + i * bytes, bytes);
                memcpy(children + (1u - side) * bytes, auth_values + cursor++ * bytes, bytes);
                ++i;
            }
            parents[parent_count] = node >> 1;
            hash(parent_values + parent_count * bytes, children, height + 1,
                 parents[parent_count], hash_ctx);
            ++parent_count;
        }
        size_t n = 0, i = 0, j = 0;
        while (i < parent_count || j < upper_count) {
            if (j == upper_count || (i < parent_count && parents[i] < upper[j])) {
                merged[n] = parents[i];
                memcpy(merged_values + n++ * bytes, parent_values + i++ * bytes, bytes);
            } else {
                merged[n] = upper[j];
                memcpy(merged_values + n++ * bytes, upper_values + j++ * bytes, bytes);
            }
        }
        memcpy(current, merged, n * sizeof(*current));
        memcpy(current_values, merged_values, n * bytes);
        current_count = n;
        upper_count = 0;
    }
    if (current_count == 1 && current[0] == 0 && cursor == auth_count) {
        memcpy(root, current_values, bytes);
        result = 0;
    }
done:
    free(values);
    free(positions);
    return result;
}
