/* compute_norm_fixedlen.c */

#include <stdint.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

/* 不再需要 length；長度隱含為 2^k */
void compute_norm(mpz_t rop,
                  const int8_t *poly,
                  unsigned k)
{
    const slong n = (slong)1 << k;          /* n = 2^k */

    /* 1) f(x) = Σ poly[i] x^i, i = 0..n-1 ---------------------------- */
    fmpz_poly_t f;
    fmpz_poly_init2(f, n);                  /* 預留 n 個係數位置 */
    for (slong i = 0; i < n; i++)
        fmpz_set_si(f->coeffs + i, (long)poly[i]);
    _fmpz_poly_set_length(f, n);
    _fmpz_poly_normalise(f);

    /* 2) g(x) = x^n + 1 --------------------------------------------- */
    fmpz_poly_t g;
    fmpz_poly_init2(g, n + 1);
    fmpz_one(g->coeffs + n);                /* x^n */
    fmpz_one(g->coeffs + 0);                /* +1  */
    _fmpz_poly_set_length(g, n + 1);

    /* 3) res = resultant(f, g) -------------------------------------- */
    fmpz_t res;
    fmpz_init(res);
    fmpz_poly_resultant(res, f, g);

    /* 4) Norm = (-1)^{n·deg f} · res -------------------------------- */
    const slong deg_f = fmpz_poly_degree(f);
    if ((n * deg_f) & 1)
        fmpz_neg(res, res);

    /* 5) 複製回 mpz_t                                                */
    fmpz_get_mpz(rop, res);

    /* 6) 清理 -------------------------------------------------------- */
    fmpz_clear(res);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
}