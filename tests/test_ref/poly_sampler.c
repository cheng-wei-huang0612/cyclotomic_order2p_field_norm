/* poly_sampler.c
 * Gaussian sampler for int8_t polynomials.
 * Each coefficient is drawn from N(0,3^2) and rounded to nearest int.
 * Any value outside [-128,127] is saturated to fit into int8_t.
 */

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Box–Muller transform helper: return a pair of independent standard
 * normal variates in z0,z1 (passed by pointer). */
static inline void box_muller(double *z0, double *z1)
{
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double mag = sqrt(-2.0 * log(u1));
    *z0 = mag * cos(2.0 * M_PI * u2);
    *z1 = mag * sin(2.0 * M_PI * u2);
}

/* Generate a polynomial of length n with i.i.d. Gaussian(0, sigma^2)
 * coefficients (sigma = 3).  The result is written into the int8_t array
 * `poly` supplied by the caller. */
void sample_poly(int8_t *poly, int n)
{
    const double sigma = 3.0;
    static int seeded = 0;
    if (!seeded)
    {
        srand((unsigned)time(NULL));
        seeded = 1;
    }

    for (int i = 0; i < n; i += 2)
    {
        double z0, z1;
        box_muller(&z0, &z1);

        int val0 = (int)llround(z0 * sigma);
        int val1 = (int)llround(z1 * sigma);

        /* Saturate to int8_t range */
        if (val0 < -128) val0 = -128;
        if (val0 > 127)  val0 = 127;
        if (val1 < -128) val1 = -128;
        if (val1 > 127)  val1 = 127;

        poly[i] = (int8_t)val0;
        if (i + 1 < n)
            poly[i + 1] = (int8_t)val1;
    }
}
