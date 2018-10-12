#ifndef RAND_H
#define RAND_H

#define NK 624
#define	MK 397
#define MATRIX_A 0x9908b0dfUL // constant vector a //
#define UMASK 0x80000000UL // most significant w-r bits //
#define LMASK 0x7fffffffUL // least significant r bits //
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

#include <math.h>

static unsigned long state[NK]; // state vector //
static int lefto = 1;
static int initf = 0;
static unsigned long *next;

/* initializes state[NK] with a seed */
void ini_rand(unsigned long s)
{
	int j;
	state[0]= s & 0xffffffffUL;

	for (j=1; j<NK; j++) {
		state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j);

		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect */
		/* only MSBs of the array state[]. */
		/* 2002/01/09 modified by Makoto Matsumoto */

		state[j] &= 0xffffffffUL; /* for >32 bit machines */
	}

	lefto = 1; initf = 1;
}
static void next_state(void)
{
	unsigned long *p=state;
	int j;

	// if ini_rand() has not been called, //
	// a default initial seed is used //

	if (initf==0) ini_rand(5489UL);

	lefto = NK;
	next = state;

	for (j=NK-MK+1; --j; p++)
		*p = p[MK] ^ TWIST(p[0], p[1]);

	for (j=MK; --j; p++)
		*p = p[MK-NK] ^ TWIST(p[0], p[1]);

	*p = p[MK-NK] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,1)-real-interval */
double random(void)
{
	unsigned long y;
	if (--lefto == 0) next_state();
	y = *next++;

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return (double)y * (1.0/4294967296.0);
	/* divided by 2^32 */
}

double rndNormal()
{
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if(iset == 0)
	{
		do
		{
			v1 = 2.0*random() - 1.0;
			v2 = 2.0*random() - 1.0;
			rsq = v1*v1 + v2*v2;
		}
		while(rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);

		gset = v1*fac;
		iset = 1;

		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

#endif