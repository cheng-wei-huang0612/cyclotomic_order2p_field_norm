#include <stdio.h>
#include <gmp.h>
#include "utils.h"

unsigned long inv_mod1;
int *root1, *inv_root1;
int *root2, *inv_root2;
int *invN_mod1, *invN_mod2;

// Fast modular exponentiation
static int pow_mod(int base, unsigned long exp, int mod) {
    long long result = 1, b = base % mod;
    while (exp > 0) {
        if (exp & 1UL) result = (result * b) % mod;
        b = (b * b) % mod;
        exp >>= 1UL;
    }
    return (int) result;
}

// Compute modular inverse of a (mod m) using GMP (since m is prime, we could also use Fermat's little theorem)
static unsigned long modinv(unsigned long a, unsigned long m) {
    mpz_t A, M, Inv;
    mpz_init_set_ui(A, a);
    mpz_init_set_ui(M, m);
    mpz_init(Inv);
    if (mpz_invert(Inv, A, M) == 0) {  // no inverse
        mpz_clears(A, M, Inv, NULL);
        return 0;
    }
    unsigned long inv_val = mpz_get_ui(Inv);
    mpz_clears(A, M, Inv, NULL);
    return inv_val;
}

void init_ntt_roots(size_t maxN) {
    // Allocate arrays of length maxN+1 for roots and their inverses
    root1 = (int*) malloc((maxN+1) * sizeof(int));
    inv_root1 = (int*) malloc((maxN+1) * sizeof(int));
    root2 = (int*) malloc((maxN+1) * sizeof(int));
    inv_root2 = (int*) malloc((maxN+1) * sizeof(int));
    invN_mod1 = (int*) malloc((maxN+1) * sizeof(int));
    invN_mod2 = (int*) malloc((maxN+1) * sizeof(int));
    // Choose a primitive root (generator) for each field
    int g1 = 3, g2 = 3;
    // Compute primitive root of unity of order maxN for each prime
    unsigned long exp1 = (MOD1 - 1) / maxN;
    unsigned long exp2 = (MOD2 - 1) / maxN;
    int omega1 = pow_mod(g1, exp1, MOD1);
    int omega2 = pow_mod(g2, exp2, MOD2);
    // Precompute all powers: rootX[i] = omegaX^i mod MODX
    root1[0] = 1;
    root2[0] = 1;
    for (size_t i = 1; i <= maxN; ++i) {
        root1[i] = (int)((long long)root1[i-1] * omega1 % MOD1);
        root2[i] = (int)((long long)root2[i-1] * omega2 % MOD2);
    }
    // Inverse roots by symmetry: inv_root[i] = root[maxN - i]
    for (size_t i = 0; i <= maxN; ++i) {
        inv_root1[i] = root1[maxN - i];
        inv_root2[i] = root2[maxN - i];
    }
    // Precompute inverse of each 1 ≤ k ≤ maxN (for normalization factors)
    invN_mod1[1] = 1;
    invN_mod2[1] = 1;
    for (size_t j = 2; j <= maxN; ++j) {
        invN_mod1[j] = (int)(((long long)(MOD1 - MOD1/j) * invN_mod1[MOD1 % j]) % MOD1);
        invN_mod2[j] = (int)(((long long)(MOD2 - MOD2/j) * invN_mod2[MOD2 % j]) % MOD2);
    }
    // Compute inv_mod1 = (MOD1)^{-1} mod MOD2 for CRT combination
    inv_mod1 = modinv(MOD1, MOD2);
}

// In-place iterative NTT (Cooley-Tukey DIT)
void NTT(int *a, size_t N, const int *root, int mod) {
    // Bit-reverse reorder
    size_t j = 0;
    for (size_t i = 1; i < N; ++i) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            int tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }
    // Iterative FFT-like transform
    for (size_t len = 1; len < N; len <<= 1) {
        size_t step = N / (2 * len);
        for (size_t i = 0; i < N; i += 2 * len) {
            for (size_t k = 0; k < len; ++k) {
                int u = a[i+k];
                int v = (int)(((long long)a[i+k+len] * root[k * step]) % mod);
                int t = u + v;
                if (t >= mod) t -= mod;
                a[i+k] = t;
                t = u - v;
                if (t < 0) t += mod;
                a[i+k+len] = t;
            }
        }
    }
}

// In-place inverse NTT (using inv_root)
void INTT(int *a, size_t N, const int *inv_root, int mod) {
    // Same structure as NTT (DIT), but using inverse roots
    size_t j = 0;
    for (size_t i = 1; i < N; ++i) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            int tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }
    for (size_t len = 1; len < N; len <<= 1) {
        size_t step = N / (2 * len);
        for (size_t i = 0; i < N; i += 2 * len) {
            for (size_t k = 0; k < len; ++k) {
                int u = a[i+k];
                int v = (int)(((long long)a[i+k+len] * inv_root[k * step]) % mod);
                int t = u + v;
                if (t >= mod) t -= mod;
                a[i+k] = t;
                t = u - v;
                if (t < 0) t += mod;
                a[i+k+len] = t;
            }
        }
    }
    // (The calling function will perform division by N)
}

