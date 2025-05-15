// red_once.c
#include "modmul.h"

/* If a does not lie in [p/2, p/2) then branchless plus or minus one */
int32_t red_once(int32_t a, int32_t p){
    const int32_t half = p >> 1;          /* floor(p/2) */

    /* mask_add  = 0xFFFFFFFF  if a < -half   (need a += p)
     * mask_sub  = 0xFFFFFFFF  if a >  half   (need a -= p)
     * otherwise 0x00000000 .
     *
     *  - (a + half) is negative  ⇔  a < -half
     *  - (a - half - 1) is >= 0  ⇔  a >  half
     *    (subtract 1 so that a==half will not trigger  -= p)
     */
    int32_t mask_add = (a + half) >> 31;              /* arithmetic shift */
    int32_t mask_sub = ~((a - half - 1) >> 31);

    a += mask_add & p;
    a -= mask_sub & p;
    return a;

}
    