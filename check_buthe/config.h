#define OPEN_BSD
#define __USE_GNU
typedef unsigned u32_t;
typedef int i32_t;
typedef unsigned long long u64_t;
typedef long long int i64_t;
typedef unsigned short u16_t;
#define U32_MAX 0xffffffff
#define mpz_set_u64 mpz_set_ui
#define mpz_get_u64 mpz_get_ui
#define mpz_fdiv_u64 mpz_fdiv_ui
#define mpz_fdiv_q_u64 mpz_fdiv_q_ui
#define mpz_mul_u64 mpz_mul_ui
#define mpz_add_u64 mpz_add_ui
#define mpz_sub_u64 mpz_sub_ui
#define gmp_vfprintf __gmp_vfprintf

#if UINT_MAX < 0xffffffff
#error unsigned must be 32 bit
#endif
#define I64_MAX 0x7fffffffffffffff
#define I64_MIN 0x8000000000000000
#define U64_MAX 0xFFFFFFFFFFFFFFFF

#define mpz_get_u64 mpz_get_ui
#define mpz_fits_u64_p mpz_fits_ulong_p

/*
 * This limb of a mpz_t should start with the 64-th bit (the least
 * siginificant being the zeroth and should be at least 32 bits long.
 */

#define BIT64_LIMB 1

typedef unsigned char  u8_t;
