#include <mpfr.h>


int zs_init(char *filename);
int zs_close(void);
int zs_get_next_zero(mpfr_t zero);
