#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include "utils.h"

static void poly_negate(poly_t *res, const poly_t *f);
static void poly_naive_negacyclic_convolution(poly_t *res, const poly_t *f, const poly_t *g);
static void poly_ntt_negacyclic_convolution(poly_t *res, const poly_t *f, const poly_t *g);
static void compute_norm(poly_t *res, const poly_t *f);

// Compute h(x) = f(-x) (negate odd-powered coefficients)
static void poly_negate(poly_t *res, const poly_t *f) {
    size_t n = f->n;
    for (size_t j = 0; j < n; ++j) {
        if (j % 2 == 1) {
            mpz_neg(res->coeff[j], f->coeff[j]);
        } else {
            mpz_set(res->coeff[j], f->coeff[j]);
        }
    }
}

// Naive O(n^2) negacyclic convolution: res = f * g mod (x^n+1)
static void poly_naive_negacyclic_convolution(poly_t *res, const poly_t *f, const poly_t *g) {
    size_t n = f->n;
    for (size_t k = 0; k < n; ++k) {
        mpz_set_ui(res->coeff[k], 0);
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            mpz_t temp;
            mpz_init(temp);
            mpz_mul(temp, f->coeff[i], g->coeff[j]);
            size_t sum_index = i + j;
            if (sum_index < n) {
                mpz_add(res->coeff[sum_index], res->coeff[sum_index], temp);
            } else {
                mpz_sub(res->coeff[sum_index - n], res->coeff[sum_index - n], temp);
            }
            mpz_clear(temp);
        }
    }
}

// NTT-based negacyclic convolution: res = f * g mod (x^n+1)
static void poly_ntt_negacyclic_convolution(poly_t *res, const poly_t *f, const poly_t *g) {
    size_t n = f->n;
    size_t N = 2 * n;
    // Prepare coefficient arrays for NTT (mod 1 and mod 2)
    int *a1 = malloc(N * sizeof(int)), *b1 = malloc(N * sizeof(int));
    int *a2 = malloc(N * sizeof(int)), *b2 = malloc(N * sizeof(int));
    for (size_t i = 0; i < n; ++i) {
        unsigned long u1 = mpz_fdiv_ui(f->coeff[i], MOD1);
        unsigned long v1 = mpz_fdiv_ui(g->coeff[i], MOD1);
        unsigned long u2 = mpz_fdiv_ui(f->coeff[i], MOD2);
        unsigned long v2 = mpz_fdiv_ui(g->coeff[i], MOD2);
        a1[i] = (int) u1;
        b1[i] = (int) v1;
        a2[i] = (int) u2;
        b2[i] = (int) v2;
    }
    for (size_t i = n; i < N; ++i) {
        a1[i] = b1[i] = 0;
        a2[i] = b2[i] = 0;
    }
    // Forward NTT on length-N arrays
    NTT(a1, N, root1, MOD1);
    NTT(b1, N, root1, MOD1);
    NTT(a2, N, root2, MOD2);
    NTT(b2, N, root2, MOD2);
    // Pointwise multiply in NTT domain
    for (size_t k = 0; k < N; ++k) {
        long long x1 = (long long)a1[k] * b1[k] % MOD1;
        long long x2 = (long long)a2[k] * b2[k] % MOD2;
        a1[k] = (int) x1;
        a2[k] = (int) x2;
    }
    // Inverse NTT for each modulus
    INTT(a1, N, inv_root1, MOD1);
    INTT(a2, N, inv_root2, MOD2);
    // Divide by N (since we performed unnormalized inverse NTT)
    int invN1 = invN_mod1[N];
    int invN2 = invN_mod2[N];
    for (size_t k = 0; k < N; ++k) {
        a1[k] = (int)(((long long)a1[k] * invN1) % MOD1);
        a2[k] = (int)(((long long)a2[k] * invN2) % MOD2);
    }
    // Combine results for res coefficients via CRT
    unsigned long inv_mod1_mod2 = inv_mod1;
    for (size_t r = 0; r < n; ++r) {
        int A1 = a1[r], A1n = a1[r+n];
        int A2 = a2[r], A2n = a2[r+n];
        int b1_val = A1 - A1n;
        if (b1_val < 0) b1_val += MOD1;
        int b2_val = A2 - A2n;
        if (b2_val < 0) b2_val += MOD2;
        long long diff = (long long)b2_val - b1_val;
        if (diff < 0) diff += MOD2;
        unsigned long k_val = (unsigned long)((diff % MOD2) * inv_mod1_mod2 % MOD2);
        unsigned long long x = (unsigned long long) b1_val + (unsigned long long) k_val * MOD1;
        // Set res->coeff[r] = x (combined result)
        mpz_set_ui(res->coeff[r], (unsigned long)(x & 0xFFFFFFFFUL));
        unsigned long high = (unsigned long)(x >> 32);
        if (high != 0) {
            mpz_add_ui(res->coeff[r], res->coeff[r], high * (1UL<<32));
        }
    }
    mpz_t combined; mpz_init(combined); mpz_clear(combined);  // (placeholder, not used further)
    free(a1); free(b1);
    free(a2); free(b2);
}

// Recursively compute Norm(f) and store in res->coeff[0]
static void compute_norm(poly_t *res, const poly_t *f) {
    size_t n = f->n;
    if (n == 1) {
        mpz_set(res->coeff[0], f->coeff[0]);
        return;
    }
    size_t m = n / 2;
    poly_t h = poly_init(n);
    poly_t g = poly_init(m);
    poly_negate(&h, f);
    if (n >= 64) {
        poly_ntt_negacyclic_convolution(&g, f, &h);
    } else {
        poly_t full = poly_init(n);
        poly_naive_negacyclic_convolution(&full, f, &h);
        for (size_t t = 0; t < m; ++t) {
            mpz_set(g.coeff[t], full.coeff[2*t]);
        }
        poly_clear(&full);
    }
    poly_clear(&h);
    poly_t rec_res = poly_init(m);
    compute_norm(&rec_res, &g);
    mpz_set(res->coeff[0], rec_res.coeff[0]);
    poly_clear(&rec_res);
    poly_clear(&g);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <n> <coeff0> <coeff1> ... <coeff{n-1}>\n", argv[0]);
        return 1;
    }
    size_t n = (size_t) atoi(argv[1]);
    if (n < 1) {
        fprintf(stderr, "Error: n must be positive and power of 2.\n");
        return 1;
    }
    size_t tmp = n;
    while (tmp % 2 == 0 && tmp > 1) tmp /= 2;
    if (tmp != 1) {
        fprintf(stderr, "Warning: n is not a power of 2.\n");
    }
    if (argc - 2 < (int)n) {
        fprintf(stderr, "Error: Expected %zu coefficients, got %d.\n", n, argc-2);
        return 1;
    }
    poly_t f = poly_init(n);
    for (size_t i = 0; i < n; ++i) {
        mpz_set_si(f.coeff[i], strtol(argv[2+i], NULL, 10));
    }
    poly_t result = poly_init(1);
    init_ntt_roots(2 * n);
    compute_norm(&result, &f);
    gmp_printf("%Zd\n", result.coeff[0]);
    poly_clear(&f);
    poly_clear(&result);
    return 0;
}

