//Barrett.c
#include "modmul.h"




int32_t barrett_mul(int32_t a, int32_t b,
                                   int32_t p, int32_t mu)
{
#if BARRETT_USE_ARM64_ASM
    uint32_t lo, hi, q;
    __asm__ volatile(
        "umull   %w[lo], %w[hi], %w[a],  %w[b] \n\t"
        "umulh   %w[q],  %w[lo], %w[mu]       \n\t"
        "msub    %w[lo], %w[q],  %w[p],  %w[lo]\n\t"
        "cmp     %w[lo], %w[p]                \n\t"
        "sub     %w[lo], %w[lo], %w[p], hs    \n\t"
        : [lo]"=&r"(lo), [hi]"=&r"(hi), [q]"=&r"(q)
        : [a]"r"(a),  [b]"r"(b),
          [p]"r"(p),  [mu]"r"(mu)
        : "cc"
    );
    return lo;
#else
    /* 可攜式 baseline  */
    uint64_t t  = (uint64_t)a * b;
    uint32_t q  = (uint32_t)((t * mu) >> 32);
    uint32_t r  = (uint32_t)(t - (uint64_t)q * p);
    if (r >= p) r -= p;
    return r;
#endif
}



int32_t barrett_mul_const(int32_t a,   /* variable */
                                         int32_t b,   /* compile-time constant */
                                         int32_t p,
                                         int32_t mu_b /* 事先算好 */)
{
#if BARRETT_USE_ARM64_ASM
    uint32_t lo, q;
    __asm__ volatile(
        "umull   %w[lo], wzr, %w[a], %w[b]   \n\t"
        "umulh   %w[q],  %w[a], %w[mu_b]     \n\t"
        "msub    %w[lo], %w[q], %w[p], %w[lo]\n\t"
        "cmp     %w[lo], %w[p]               \n\t"
        "sub     %w[lo], %w[lo], %w[p], hs   \n\t"
        : [lo]"=&r"(lo), [q]"=&r"(q)
        : [a]"r"(a), [b]"r"(b),
          [p]"r"(p), [mu_b]"r"(mu_b)
        : "cc"
    );
    return lo;
#else
    uint64_t t  = (uint64_t)a * b;
    uint32_t q  = (uint32_t)(((uint64_t)a * mu_b) >> 32);
    uint32_t r  = (uint32_t)(t - (uint64_t)q * p);
    if (r >= p) r -= p;
    return r;
#endif
}
