/*
------------------------------------------------------------------------------
isaac64.h: definitions for a random number generator
Bob Jenkins, 1996, Public Domain
------------------------------------------------------------------------------
*/
/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef rng_isaac_h_included
#define rng_isaac_h_included

/*
11/25/05: EWM -  typedefs to use standard int types defined in types.h :
*/
#include	"Mdata.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef  uint64  ub8;
#define UB8MAXVAL 0xffffffffffffffffLL
#define UB8BITS 64
typedef  sint64  sb8;
#define SB8MAXVAL 0x7fffffffffffffffLL
typedef  uint32  ub4;	/* unsigned 4-byte quantities */
#define UB4MAXVAL 0xffffffff
typedef  sint32  sb4;
#define UB4BITS 32
#define SB4MAXVAL 0x7fffffff
typedef  uint16  ub2;
#define UB2MAXVAL 0xffff
#define UB2BITS 16
typedef  sint16  sb2;
#define SB2MAXVAL 0x7fff
typedef uint8	 ub1;
#define UB1MAXVAL 0xff
#define UB1BITS 8
typedef sint8	 sb1;	/* signed 1-byte quantities */
#define SB1MAXVAL 0x7f
typedef int  	word;	/* fastest type available */


#ifndef ISAAC64
#define ISAAC64

#define RANDSIZL   (8)
#define RANDSIZ    (1<<RANDSIZL)

extern ub8 randrsl[RANDSIZ], randcnt;

/*
------------------------------------------------------------------------------
 If (flag==TRUE), then use the contents of randrsl[0..255] as the seed.
------------------------------------------------------------------------------
*/
void rng_isaac_init(word flag);

void isaac64();

/*
11/25/05: EWM - modified to add 2 types of double-precision floating rand() calls:

	- rng_isaac_rand_double() returns a random double via a 64-bit field
	which is (within the limits of the generator) a random 64-bit int;

	- rng_isaac_rand_double_norm_pos() returns a random double with
	probability uniformly distributed in [0, 1), insofar as IEEE64 doubles
	are capable of distributing such values, excluding underflows;

	- rng_isaac_rand_double_norm_pm1() returns a random double with
	probability uniformly distributed in (-1, 1), insofar as IEEE64 doubles
	are capable of distributing such values, excluding underflows;
*/
double	rng_isaac_rand_double();

double	rng_isaac_rand_double_norm_pos();

double	rng_isaac_rand_double_norm_pm1();

/*
------------------------------------------------------------------------------
 Call rand() to retrieve a single 64-bit random value
------------------------------------------------------------------------------
*/
#define rng_isaac_rand() \
   (!randcnt-- ? (isaac64(), randcnt=RANDSIZ-1, randrsl[randcnt]) : \
                 randrsl[randcnt])

#endif  /* RAND */

#ifdef __cplusplus
}
#endif

#endif	/* rng_isaac_h_included */

