#include "poly.h"

/* Scheme presets live in the compatibility/algorithm target. Algebra's
 * ntt_plan_init_256 knows only the numerical configurations it supports. */
void ntt_init_dili(lp_ctx *ctx)
{
    const ring_ctx ring = {
        {8380417, UINT32_C(58728449), 4193792, 2365951},
        256, RING_XN_PLUS_ONE
    };
    (void)ntt_plan_init_256(ctx, &ring, 8);
}

void ntt_init_kyber(lp_ctx *ctx)
{
    const ring_ctx ring = {
        {3329, UINT32_C(1806234369), 1353, 2988},
        256, RING_XN_PLUS_ONE
    };
    (void)ntt_plan_init_256(ctx, &ring, 7);
}
