/* Number of 16 bit words in a q type number (12 or 24) */
#define NQ 24

/* Number of words in significand area */
#define OMG (NQ-2)

/* Byte offset to least significant word of significand */
#define OFFS (2*OMG+2)

/* Number of bits of precision */
#define NBITS ((OMG-1)*16)

/* Maximum number of decimal digits in conversion */
#define NDEC (NBITS*8/27)
