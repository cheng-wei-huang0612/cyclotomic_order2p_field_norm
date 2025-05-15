#ifndef UTILS_H
#define UTILS_H

#include <gmp.h>
#include <stdlib.h>

typedef struct {
    size_t n;
    mpz_t *coeff;
} poly_t;

// Allocate a polynomial of length n (coeff initialized to 0)
static inline poly_t poly_init(size_t n) {
    poly_t p;
    p.n = n;
    p.coeff = (mpz_t*) malloc(n * sizeof(mpz_t));
    for (size_t i = 0; i < n; ++i) {
        mpz_init(p.coeff[i]);
    }
    return p;
}

// Free polynomial coefficients
static inline void poly_clear(poly_t *p) {
    for (size_t i = 0; i < p->n; ++i) {
        mpz_clear(p->coeff[i]);
    }
    free(p->coeff);
    p->coeff = NULL;
    p->n = 0;
}

// NTT parameters and functions
#define MOD1 998244353
#define MOD2 469762049
extern unsigned long inv_mod1;
extern int *root1;
extern int *inv_root1;
extern int *root2;
extern int *inv_root2;
extern int *invN_mod1;
extern int *invN_mod2;

void init_ntt_roots(size_t maxN);
void NTT(int *a, size_t N, const int *root, int mod);
void INTT(int *a, size_t N, const int *inv_root, int mod);

#endif

