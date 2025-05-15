/* test_ref.c
 * Quick‐and‐dirty correctness test for compute_norm().
 * For each n = 2^k, k = 2..10 (4..1024), generate two random
 * Gaussian polynomials with sigma = 3 and print their norms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>
#include "flint/flint.h"

/* External functions implemented elsewhere */
void sample_poly(int8_t *poly, int n);
void compute_norm(mpz_t rop, const int8_t *poly, unsigned k);

int main(void)
{
    const unsigned k_vals[] = {2,3,4,5,6,7,8,9,10};
    const size_t   num_k    = sizeof(k_vals)/sizeof(k_vals[0]);

    for (size_t idx = 0; idx < num_k; idx++)
    {
        unsigned k = k_vals[idx];
        int n      = 1 << k;
        printf("================ n = %d (k = %u) ================\n", n, k);

        for (int s = 0; s < 2; s++)
        {
            int8_t *poly = malloc((size_t)n * sizeof(int8_t));
            if (!poly)
            {
                perror("malloc");
                return 1;
            }

            sample_poly(poly, n);

            mpz_t norm;
            mpz_init(norm);
            compute_norm(norm, poly, k);

            printf("Sample %d, first 16 coeffs: [", s + 1);
            int show = n < 16 ? n : 16;
            for (int i = 0; i < show; i++)
            {
                printf("%d", poly[i]);
                if (i + 1 < show) printf(", ");
            }
            if (n > 16) printf(", ...");
            printf("]\nNorm = ");
            gmp_printf("%Zd\n\n", norm);

            mpz_clear(norm);
            free(poly);
        }
    }

    return 0;
}
