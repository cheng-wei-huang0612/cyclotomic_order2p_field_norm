#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include "utils.h"

// Compute field norm via schoolbook method (matrix determinant)
static void build_matrix(const poly_t *f, mpz_t **M) {
    size_t n = f->n;
    // M is an n×n matrix (mpz_t entries) for multiplication-by-f
    for (size_t k = 0; k < n; ++k) {
        // Compute coefficients of f(α) * α^k mod (α^n + 1)
        mpz_t *res = (mpz_t*) malloc(n * sizeof(mpz_t));
        for (size_t i = 0; i < n; ++i) mpz_init(res[i]);
        for (size_t j = 0; j < n; ++j) {
            if (mpz_sgn(f->coeff[j]) == 0) continue;
            size_t exp = j + k;
            if (exp >= n) {
                exp -= n;
                mpz_sub(res[exp], res[exp], f->coeff[j]);  // wrap with - sign
            } else {
                mpz_add(res[exp], res[exp], f->coeff[j]);
            }
        }
        // Copy result into matrix column k
        for (size_t m = 0; m < n; ++m) {
            mpz_set(M[m][k], res[m]);
        }
        for (size_t i = 0; i < n; ++i) {
            mpz_clear(res[i]);
        }
        free(res);
    }
}

// Compute determinant of M (n×n) using Bareiss algorithm (fraction-free Gaussian elimination)
static void compute_determinant(mpz_t det, mpz_t **M, size_t n) {
    mpz_set_ui(det, 1);
    if (n == 0) return;
    if (n == 1) { mpz_set(det, M[0][0]); return; }
    mpz_t temp1, temp2, piv, prev_piv;
    mpz_inits(temp1, temp2, piv, prev_piv, NULL);
    mpz_set_ui(prev_piv, 1);
    for (size_t k = 0; k < n; ++k) {
        mpz_set(piv, M[k][k]);
        if (mpz_sgn(piv) == 0) {
            // Pivot 0: swap with a non-zero row
            size_t swap_row = k+1;
            while (swap_row < n && mpz_sgn(M[swap_row][k]) == 0) swap_row++;
            if (swap_row == n) { mpz_set_ui(det, 0); goto cleanup; }
            for (size_t j = k; j < n; ++j) {
                mpz_swap(M[k][j], M[swap_row][j]);
            }
            mpz_set(piv, M[k][k]);
        }
        if (k == n-1) break;
        for (size_t i = k+1; i < n; ++i) {
            for (size_t j = k+1; j < n; ++j) {
                // M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) / prev_piv
                mpz_mul(temp1, M[i][j], piv);
                mpz_mul(temp2, M[i][k], M[k][j]);
                mpz_sub(temp1, temp1, temp2);
                if (k > 0) mpz_divexact(temp1, temp1, prev_piv);
                mpz_set(M[i][j], temp1);
            }
        }
        mpz_set(prev_piv, piv);
    }
    mpz_set(det, M[n-1][n-1]);
cleanup:
    mpz_clears(temp1, temp2, piv, prev_piv, NULL);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <n> <coeff0> <coeff1> ... <coeff{n-1}>\n", argv[0]);
        return 1;
    }
    size_t n = (size_t) atoi(argv[1]);
    if (n < 1) {
        fprintf(stderr, "Error: n must be positive and a power of 2.\n");
        return 1;
    }
    // (Check if n is power of 2)
    size_t temp_n = n;
    while (temp_n % 2 == 0 && temp_n > 1) temp_n /= 2;
    if (temp_n != 1) {
        fprintf(stderr, "Warning: n is not a power of 2 (norm is defined for cyclotomic cases).\n");
    }
    if (argc - 2 < (int)n) {
        fprintf(stderr, "Error: Expected %zu coefficients, but got %d.\n", n, argc-2);
        return 1;
    }
    // Initialize polynomial f
    poly_t f;
    f.n = n;
    f.coeff = malloc(n * sizeof(mpz_t));
    for (size_t i = 0; i < n; ++i) {
        mpz_init(f.coeff[i]);
        mpz_set_si(f.coeff[i], strtol(argv[2+i], NULL, 10));
    }
    // Allocate and build matrix M
    mpz_t **M = (mpz_t**) malloc(n * sizeof(mpz_t*));
    for (size_t i = 0; i < n; ++i) {
        M[i] = (mpz_t*) malloc(n * sizeof(mpz_t));
        for (size_t j = 0; j < n; ++j) mpz_init(M[i][j]);
    }
    build_matrix(&f, M);
    // Compute determinant (field norm)
    mpz_t det;
    mpz_init(det);
    compute_determinant(det, M, n);
    gmp_printf("%Zd\n", det);
    // Cleanup
    mpz_clear(det);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) mpz_clear(M[i][j]);
        free(M[i]);
    }
    free(M);
    for (size_t i = 0; i < n; ++i) mpz_clear(f.coeff[i]);
    free(f.coeff);
    return 0;
}

