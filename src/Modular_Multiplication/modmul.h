// modmul.h

#include <stdint.h>

#if !defined(BARRETT_USE_ARM64_ASM)

# if defined(__aarch64__)
#  define BARRETT_USE_ARM64_ASM 0
# else
#  define BARRETT_USE_ARM64_ASM 0
# endif
#endif

int32_t red_once(int32_t a, int32_t p);

int32_t barrett_mul(int32_t a, int32_t b, int32_t p, int32_t mu);
int32_t barrett_mul_const(int32_t a, int32_t b, int32_t p, int32_t mu_b);


