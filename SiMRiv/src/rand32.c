//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File:   rand32.c
// Author: Bob Harris   rsharris at bx dot psu dot edu
//
// Support code for generic C applications
//
//----------

//----------
//
// rand32--
//	Generate 32-bit random numbers using the Mersenne Twister algorithm.
//
//----------
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Library General Public License as published by
// the Free Software Foundation (either version 2 of the License or, at your
// option, any later version).  This library is distributed in the hope that
// it will be useful, but WITHOUT ANY WARRANTY, without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
// the GNU Library General Public License for more details.  You should have
// received a copy of the GNU Library General Public License along with this
// library; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307, USA.
//
// The code as Shawn received it included the following notice:
//
//   Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When
//   you use this, send an e-mail to <matumoto@math.keio.ac.jp> with
//   an appropriate reference to your work.
//
// It would be nice to CC: <Cokus@math.washington.edu> when you write.
//
//----------
//
// The MT (Mersenne Twister) MT19937 generator was invented by Makoto Matsumoto
// and Takuji Nishimura.  It generates pseudorandom integers uniformly
// distributed in the range 0..(2^32 - 1), starting from any seed in the range
// 0..(2^31 - 1) (be sure to note that's a 31-bit seed).  For further details,
// see
//		www.math.keio.ac.jp/~matumoto/emt.html
//
// The Mersenne Twister is designed with consideration of the flaws of various
// existing generators, has a period of 2^19937 - 1, gives a sequence that is
// 623-dimensionally equidistributed, and has passed many stringent tests,
// including the die-hard test of G. Marsaglia and the load test of P.
// Hellekalek and S. Wegenkittl.  It is efficient in memory usage (typically
// using 2506 to 5012 bytes of static data, depending on data type sizes, and
// the code is quite short as well).  It generates random numbers in batches of
// 624 at a time (though the user need not be aware of this), so the caching
// and pipelining of modern systems is exploited.  It is also divide- and
// modulo-free.
//
// This implementation is based on one by Shawn Cokus (Cokus@math.washington.
// edu) on Mar/8/98.  That in turn was based on one by Takuji Nishimura (who
// also credits Topher Cooper and Marc Rieffel), Aug/97.  The Cokus code can
// be found at 
//		www.math.keio.ac.jp/~matumoto/cokus.c
// and the Nishimura code can be found at
//		www.math.keio.ac.jp/~matumoto/mt19937int.c
// [[ as of Jul/26/2005, that size no longer exists;  I was then able to ]]
// [[ locate what is probably the cokus code at                          ]]
// [[   http://www.physics.uc.edu/~pinskia/dla/src/cokus.c               ]]
//
// The code found at those two sites contains GPL copyright and free
// distribution notices, which are reproduced above.
//
// The limit to a 31-bit seed was introduced by Shawn Cokus.
//
//----------

#include <stdint.h>				// standard C sized integer stuff
#include "rand32.h"				// interface to this module

//----------
//
// definitions and global data
//
//----------

typedef uint32_t u32;
typedef uint64_t u64;

#define TwoTo32		4294967296.0

// period parameters

#define N			624					// length of state vector
#define M			397					// a period parameter
#define K			0x9908B0DFU			// a magic constant
#define hiBit(u)	((u) & 0x80000000U)	// mask all but highest   bit of u
#define loBit(u)	((u) & 0x00000001U)	// mask all but lowest    bit of u
#define loBits(u)	((u) & 0x7FFFFFFFU)	// mask     the highest   bit of u
#define mixBits(u,v) (hiBit(u)|loBits(v))//move hi bit of u to hi bit of v

// tempering parameters

#define TEMPERING_MASK_B		0x9D2C5680U
#define TEMPERING_MASK_C		0xEFC60000U
#define TEMPERING_SHIFT_U(y)	((y) >> 11)
#define TEMPERING_SHIFT_S(y)	((y) << 7)
#define TEMPERING_SHIFT_T(y)	((y) << 15)
#define TEMPERING_SHIFT_L(y)	((y) >> 18)

// the state vector

static u32	state[N];		// state vector
static u32*	next;			// pointer into state vector to access the next
							// .. random value
static int	left = -1;      // number of unused items remaining in the state
							// .. vector (-1 indicates that the state vector
							// has not been initialized)

//----------
//
// srand32--
//  Seed the generator, from a number.   Seeding from a number does not
//	provides as much seed variety as seeding from a string (see ssrand32).
//
//	Note that this should be called before calling any of the other routines.
//	But if not, the default seed will be used.
//
//----------
//
// Arguments:
//	unsigned long	seed:	The seed, which should be non-zero.
//
// Returns:
//	(nothing)
//
//----------
//
// The state vector is initialized using the generator
//		x <- 69069x mod 2^32
// which is from line 15 of table 1, p. 106, Sec. 3.3.4 of Knuth's "The Art of
// Computer Programming", Volume 2, 3rd ed.
//
// (much of the following is due to Shaw Cokus, with some changes by me)
//
// It seems this seeding generator could be better.  It achieves the maximum
// period for its modulus (2^32) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2,
// Knuth);  if x_initial can be even, you have sequences like 0, 0, 0, ...;
// 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31, 2^29, 2^29
// + 2^31, ..., etc.  So the seed is forced to be odd by shifting it left one
// bit and adding one.
//
// Even if x_initial is odd, if x_initial is 1 mod 4 then
//
//   the          lowest bit of x is always 1,
//   the  next-to-lowest bit of x is always 0,
//   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
//   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
//   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
//    ...
//
// and if x_initial is 3 mod 4 then
//
//   the          lowest bit of x is always 1,
//   the  next-to-lowest bit of x is always 1,
//   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
//   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
//   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
//    ...
//
// The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
// 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
// also does well in the dimension 2..5 spectral tests, but it could be
// better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
//
// Note that the random number user does not see the values generated
// here directly since reloadMT() will always munge them first, so maybe
// none of all of this matters.  In fact, the seed values made here could
// even be extra-special desirable if the Mersenne Twister theory says
// so-- that's why the only change I made is to restrict to odd seeds.
//
//----------

void srand32
   (unsigned long	seed)
	{
	u32		x = (seed << 1) + 1;
	u32*	s = state;
	int		j;

	*(s++) = x;
	for (j=N ; --j>0 ; )
		*(s++) = (x *= 69069U);

	left = 0;
	}

//----------
//
// ssrand32--
//  Seed the generator, from a string.   Seeding from a string provides more
//	seed variety than seeding from a number (see srand32).
//
//	Note that this should be called before calling any of the other routines.
//	But if not, the default seed will be used.
//
//----------
//
// Arguments:
//	char*	seed:	The seed, which should be a zero-terminated string.  The
//					string can be any non-zero length, but characters beyond
//					the Nth (624th) are ignored.  If it is empty, the default
//					numeric seed is used.
//
// Returns:
//	(nothing)
//
//----------
//
// The state vector is initialized using the generator
//		x <- ((69069x + seed[i]) | 1) mod 2^32
// When the end of the seed string is reached, the remaining seed[i] values
// are considered to be zero.
//
//----------

void ssrand32
   (char*	seed)
	{
	u32		x = 13013;
	u32*	s = state;
	int		j;

	if (*seed == 0)
		{ srand32 (0); return; }

	for (j=N ; j-->0 ; )
		{
		if (*seed != 0) seed++;
		x = (x * 69069U) + *seed;  *(s++) = x;  x |= 1;
		}

	left = 0;
	}

//----------
//
// rand32--
//  Generate one 32-bit random number.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A random 32-bit value.
//
//----------

unsigned long rand32
   (void)
	{
	u32		y;

	// if we have any more words available, use the next one to create a
	// 'random' value

	if (--left >= 0)
		{
		// get the next value from the state vector

	next_value:
	    y  = *(next++);

		// perform tempering transformation upon it, and return it to the user

	    y ^= (y >> 11);
	    y ^= (y <<  7) & 0x9D2C5680;
	    y ^= (y << 15) & 0xEFC60000;
	    return y ^ (y >> 18);
		}

	// we don't have any more words available, so we want to generate the next
	// batch, but if the user forgot to seed the darn thing, do it for him
	// (using the default seed that results from seed 0)

	if (left < -1)
		srand32 (0);

	// generate the next batch

	{
    u32*	p0 = state;
    u32*	p2 = state+2;
	u32*	pM = state + M;
	u32		s0, s1;
    int		j;

    for(s0=state[0], s1=state[1], j=N-M+1; --j; s0=s1, s1=*(p2++))
        *(p0++) = *(pM++) ^ (mixBits(s0,s1) >> 1) ^ (loBit(s1)? K : 0);

    for(pM=state, j=M; --j; s0=s1, s1=*(p2++))
        *(p0++) = *(pM++) ^ (mixBits(s0,s1) >> 1) ^ (loBit(s1)? K : 0);

	s1 = state[0];
    *p0 = *pM ^ (mixBits(s0,s1) >> 1) ^ (loBit(s1)? K : 0);
	}

	// go get the first value from this batch

	next = state;
    left = N-1;
	goto next_value;
	}

//----------
//
// rrand32--
//  Generate a uniformly distributed signed random number in a specified
//	range.
// urand32--
//  Generate a uniformly distributed unsigned random number in a specified
//	range.
//
//----------
//
// Arguments:
//	(unsigned) long	lo:	The low number of the desired range.
//	(unsigned) long	hi:	The high number of the desired range (inclusive).
//
// Returns:
//	A random value x such that lo <= x <= hi.
//
//----------

long rrand32
   (long	lo,
	long	hi)
	{
	u32		ulo = 0x80000000 + ((u32) lo);
	u32		uhi = 0x80000000 + ((u32) hi);

	return (long) (urand32 (ulo, uhi)) - 0x80000000;
 	}

unsigned long urand32
   (unsigned long	lo,
	unsigned long	hi)
	{
	static u32 x = 0;			// upon entry, x is a uniform random
	static u32 r = 1;			// .. integer such that 0 <= x < r
	u32 s = (hi - lo) + 1;		// s = our selection range
	u64 xx, rr;
	u32 leftovers, v;

	// fprintf (stderr, "%u-way choice (%lu..%lu) from 0<=%u<%u\n", s, lo, hi, x, r);

	// handle special cases where the selection range is the entire 32-bit
	// range;  in this case we can just keep the state we have and return a
	// newly generated random value

	if (s == 0) 				// (0 is the result of overflow)
		return rand32 ();		// selection range is the entire 32-bit range

	if (s == 1)					// lo == hi so no random choice needed
		return lo;

	// loop until we have a uniform value over a range that's a whole multiple
	// of the selection range;  note that the probability, on each pass, of
	// exiting the loop is always better than 1/2 (usually much better)

	xx = x;  rr = r;

	while (1==1)
		{
		// produce random xx such that  0 <= xx < rr  and  rr >= s

		// fprintf (stderr, "0<=%llu<%llu\n", xx, rr);

		if (rr < s)
			{
			xx = (xx << 32) + (u32) rand32 ();
			rr =  rr << 32;			// [[ s < 2^32 <= rr < s*2^32 ]]
			}

		// conceptually divide the [0,rr[ interval into packets of size s and
		// some leftovers;  if xx is in a leftover (there are rr % s of them),
		// we can't make a uniform choice, so we have to try again;  we can
		// use the value xx going forward, as it is uniform in [0,leftovers[

		leftovers = rr % s;			// [[ leftovers < s < 2^32 ]]

		// fprintf (stderr, "0<=%llu<%llu with %u leftovers\n", xx, rr, leftovers);

		if (xx >= leftovers)
			break;

		rr = leftovers;				// xx is a uniform random value in [0,rr[
									// [[ rr < 2^32 ]]
		}

	// xx is in one of the packets;  we can choose its position in the packet
	// to select our return value, and keep the packet number as our random
	// value moving forward, leaving us with x being uniform in [0,r[

	xx -= leftovers;
	rr -= leftovers;

	v = xx % s;				// v = position within packet
	x = xx / s;				// x = packet number
	r = rr / s;				// r = number of packets

	return lo + v;
 	}

//----------
//
// frand32--
//  Generate a real random number in the interval [0,1).
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A random value x such that 0 <= x < 1.
//
//----------

float frand32
   (void)
	{
	return rand32 () / TwoTo32;
	}

//----------
//
// crand32--
//  Generate a weighted coin flip.
//
//----------
//
// Arguments:
//	float	p:	The probability of a 'true' result.
//
// Returns:
//	true (with probability p) or false (with probability 1-p).
//
//----------

int crand32
   (float	p)
	{
	return (rand32 () < (u32) (p * TwoTo32));
 	}

//----------
//
// drand32--
//  Generate a random value according to a given distribution.
//
//----------
//
// Arguments:
//	int		piValues:	The number of different values to choose from.  This
//						should be postive, or strange things may result.
//	float*	piSum:		The 'integral' of the probability distribution.  This
//						must be monotonically increasing, so that
//						  probability(g) = piSum[g] - piSum[g-1].
//
// Returns:
//	A value in the range 0 through piValues-1.
//
//----------

int drand32
   (int		piValues,
	unsigned long *piSum)
	{
	long	lo, hi, mid;
	float	r;

	// choose a random value in the range of the probability function

	lo = 0;
	hi = piValues - 1;

	r = piSum[hi] * frand32 ();

	// use binary search to locate the appropriate cell

	while (lo < hi)
		{
		mid = (lo + hi) / 2;
		if (r < piSum[mid]) hi = mid; else lo = mid + 1;
		}
	return lo;
	}

