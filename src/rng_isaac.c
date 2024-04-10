/*
------------------------------------------------------------------------------
isaac64.c: My random number generator for 64-bit machines.
By Bob Jenkins, 1996.  Public Domain.
------------------------------------------------------------------------------
*/

#include <stdio.h>
#include "rng_isaac.h"

/* externs declared in rng_isaac.h: */
ub8 randrsl[RANDSIZ], randcnt;

static    ub8 mm[RANDSIZ];
static    ub8 aa=0, bb=0, cc=0;

#define ind(mm,x)  (*(ub8 *)((ub1 *)(mm) + ((x) & ((RANDSIZ-1)<<3))))
#define rngstep(mix,a,b,mm,m,m2,r,x) \
{ \
  x = *m;  \
  a = (mix) + *(m2++); \
  *(m++) = y = ind(mm,x) + a + b; \
  *(r++) = b = ind(mm,y>>RANDSIZL) + x; \
}

void isaac64()
{
  register ub8 a,b,x,y,*m,*m2,*r,*mend;
  r = randrsl;	/* Need a variable address pointer to feed to rngstep */
  a = aa; b = bb + (++cc);
  for (m = mm, mend = m2 = m+(RANDSIZ/2); m<mend; )
  {
    rngstep(~(a^(a<<21)), a, b, mm, m, m2, r, x);
    rngstep(  a^(a>>5)  , a, b, mm, m, m2, r, x);
    rngstep(  a^(a<<12) , a, b, mm, m, m2, r, x);
    rngstep(  a^(a>>33) , a, b, mm, m, m2, r, x);
  }
  for (m2 = mm; m2<mend; )
  {
    rngstep(~(a^(a<<21)), a, b, mm, m, m2, r, x);
    rngstep(  a^(a>>5)  , a, b, mm, m, m2, r, x);
    rngstep(  a^(a<<12) , a, b, mm, m, m2, r, x);
    rngstep(  a^(a>>33) , a, b, mm, m, m2, r, x);
  }
  bb = b; aa = a;
}

#define mix(a,b,c,d,e,f,g,h) \
{ \
   a-=e; f^=h>>9;  h+=a; \
   b-=f; g^=a<<9;  a+=b; \
   c-=g; h^=b>>23; b+=c; \
   d-=h; a^=c<<15; c+=d; \
   e-=a; b^=d>>14; d+=e; \
   f-=b; c^=e<<20; e+=f; \
   g-=c; d^=f>>17; f+=g; \
   h-=d; e^=g<<14; g+=h; \
}

void rng_isaac_init(word flag)
{
   word i;
   ub8 a,b,c,d,e,f,g,h;
   aa=bb=cc=(ub8)0;
   a=b=c=d=e=f=g=h=0x9E3779B97F4A7C13ull;  /* the golden ratio */

   for (i=0; i<4; ++i)                    /* scramble it */
   {
     mix(a,b,c,d,e,f,g,h);
   }

   for (i=0; i<RANDSIZ; i+=8)   /* fill in mm[] with messy stuff */
   {
     if (flag)                  /* use all the information in the seed */
     {
       a+=randrsl[i  ]; b+=randrsl[i+1]; c+=randrsl[i+2]; d+=randrsl[i+3];
       e+=randrsl[i+4]; f+=randrsl[i+5]; g+=randrsl[i+6]; h+=randrsl[i+7];
     }
     mix(a,b,c,d,e,f,g,h);
     mm[i  ]=a; mm[i+1]=b; mm[i+2]=c; mm[i+3]=d;
     mm[i+4]=e; mm[i+5]=f; mm[i+6]=g; mm[i+7]=h;
   }

   if (flag)
   {        /* do a second pass to make all of the seed affect all of mm */
     for (i=0; i<RANDSIZ; i+=8)
     {
       a+=mm[i  ]; b+=mm[i+1]; c+=mm[i+2]; d+=mm[i+3];
       e+=mm[i+4]; f+=mm[i+5]; g+=mm[i+6]; h+=mm[i+7];
       mix(a,b,c,d,e,f,g,h);
       mm[i  ]=a; mm[i+1]=b; mm[i+2]=c; mm[i+3]=d;
       mm[i+4]=e; mm[i+5]=f; mm[i+6]=g; mm[i+7]=h;
     }
   }

   isaac64();          /* fill in the first set of results */
   randcnt=RANDSIZ;    /* prepare to use the first set of results */
}


#ifdef NEVER
/*int main() - Rename this for build purposes, since even with NEVER undefined, MSVC still views this as global main() */
int rng_isaac_main()
{
  word i,j;
  aa=bb=cc=(ub8)0;
  for (i=0; i<RANDSIZ; ++i) mm[i]=(ub8)0;
  rng_isaac_init(TRUE);
  for (i=0; i<2; ++i)
  {
    isaac64();
    for (j=0; j<RANDSIZ; ++j)
      printf("%.8lx%.8lx",(ub4)(randrsl[j]>>32),(ub4)randrsl[j]);
  }
}
#endif

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
double	rng_isaac_rand_double()
{
	uint64 iran64;
	uint32 fexp;

	/* Make sure resulting float will not be denormal: */
	for(;;)
	{
		iran64 = rng_isaac_rand();
		fexp = (uint32)(iran64 >> 52) & 0x7ff;
		if(fexp != 0 && fexp < 0x7f0) break;
	}
	return	*(double *)&iran64;
}

/* Assumes IEEE64-compliant: */
double	rng_isaac_rand_double_norm_pos()
{
	/*
	Obtain a result in [0, 1) by merging a sign/exponent field = 0x3ff with
	random 52-bit mantissa (52-bit because the hidden bit is assumed 1 via the
	choice of exponent - we only randomly generate the non-hidden 52 bits),
	yielding a result in [1, 2), and subtracting 1:
	*/
	uint64 iran64, itmp64;
	double retval;

	itmp64 = rng_isaac_rand();
	iran64 = 0x3FF0000000000000ull + (itmp64 & 0x000FFFFFFFFFFFFFull);
	retval=(*(double *)&iran64) - 1.0;
	/* GCC compiler bug: needed to insert the explicit range-check here, otherwise compiler 'optimized' the (*(double *)&iran64) to zero: */
	if(retval < 0.0 || retval > 1.0)
	{
		sprintf(cbuf, "rng_isaac_rand_double_norm_pos: itmp64 = %16llx, iran64 = %16llx, retval = %lf not in [0,1]!\n", itmp64, iran64, retval);
		ASSERT(0, cbuf);
	}
	return retval;
}


/* Assumes IEEE64-compliant: */
double	rng_isaac_rand_double_norm_pm1()
{
	/*
	Obtain a result in (-1, 1) by following the same procedure used in
	rng_isaac_rand_double_norm_pos to get a value in [0, 1) and multiplying
	the result by a random choice of -1 or +1. Note that this doubles the
	odds of getting a zero result, but we assume that won't be fatal -
	in essence one can consider that as though -0.0 and +0.0 were separate
	possible outputs, each occurring with probability equal to that of any
	of the discrete nonzero outputs.
	*/
	static double pm1[] = {-1.0, +1.0};
	double sign;
	uint64 itmp64, iran64;
	double retval;

	itmp64 = rng_isaac_rand();
	sign = pm1[itmp64 >> 63];	/* Use high bit of iran64 for sign */
	iran64 = 0x3FF0000000000000ull + (itmp64 & 0x000FFFFFFFFFFFFFull);
	retval=sign*((*(double *)&iran64) - 1.0);
	/* GCC compiler bug: needed to insert the explicit range-check here, otherwise compiler 'optimized' the (*(double *)&iran64) to zero: */
	if(retval < -1.0 || retval > 1.0)
	{
		sprintf(cbuf, "rng_isaac_rand_double_norm_pm1: itmp64 = %16llx, iran64 = %16llx, retval = %lf not in [0,1]!\n", itmp64, iran64, retval);
		ASSERT(0, cbuf);
	}
	return retval;
}

