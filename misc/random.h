#ifndef SMALLSECRETLWE_RANDOM_H
#define SMALLSECRETLWE_RANDOM_H

#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

static uint64_t random_x=123456789u, random_y=362436069u, random_z=521288629u;

// Seeds the rng.
static void xorshf96_random_seed(uint64_t i){
	random_x += i;
	random_y = random_x*4095834;
	random_z = random_x + random_y*98798234;
}

// Main function that generates a random 64bit integer.
static uint64_t xorshf96() {          //period 2^96-1
	uint64_t t;
	random_x ^= random_x << 16u;
	random_x ^= random_x >> 5u;
	random_x ^= random_x << 1u;

	t = random_x;
	random_x = random_y;
	random_y = random_z;
	random_z = t ^ random_x ^ random_y;

	return random_z;
}

/* Fills a given buffer with random bytes. n = size of buffer in bytes, */
static int xorshf96_fastrandombytes(void *buf, size_t n){
	uint64_t *a = (uint64_t *)buf;

	const uint32_t rest = n%8;
	const size_t limit = n/8;
	size_t i = 0;

	for (; i < limit; ++i) {
		a[i] = xorshf96();
	}

	// last limb
	uint8_t *b = (uint8_t *)buf;
	b += n - rest;
	uint64_t limb = xorshf96();
	for (size_t j = 0; j < rest; ++j) {
		b[j] = (limb >> (j*8u)) & 0xFFu;
	}

	return 0;
}

// Wrapper function for xorshf96_fastrandombytes that takes an array (pointer)
// of type uint64_t instead of void as an input because the main function generates
// always random 64bit words.  n = size of buffer in bytes.
inline static int xorshf96_fastrandombytes_uint64_array(uint64_t *buf, size_t n){
	void *a = (void *)buf;
	xorshf96_fastrandombytes(a, n);
	return 0;
}

// Wrapper function for xorshf96.
static uint64_t xorshf96_fastrandombytes_uint64() {
	return xorshf96();
}

/*
pgc (permuted congruential generator?)
*/

struct pcg_state_setseq_64 {    // Internals are *Private*.
	uint64_t state;             // RNG state.  All values are possible.
	uint64_t inc;               // Controls which RNG sequence (stream) is
	// selected. Must *always* be odd.
};

typedef struct pcg_state_setseq_64 pcg32_random_t;

static pcg32_random_t pcg32_global = { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL };

static inline uint32_t pcg32_random_r(pcg32_random_t* rng) noexcept {
	uint64_t oldstate = rng->state;
	rng->state = oldstate * 6364136223846793005ULL + rng->inc;
	uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
	uint32_t rot = oldstate >> 59u;
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static inline uint32_t pcg32_random(void) noexcept {
	return pcg32_random_r(&pcg32_global);
}


static inline void pcg32_init_state(const uint32_t state) noexcept {
	pcg32_global.state = state;
}

static inline void pcg32_init_inc(const uint32_t inc) noexcept {
	pcg32_global.inc = inc | 1;
}

// bounded version.
template<const int range>
static int fastrandombytes_int() noexcept {
	uint64_t random32bit, multiresult;
	uint32_t leftover;
	uint32_t threshold;
	random32bit =  pcg32_random();
	multiresult = random32bit * range;
	leftover = (uint32_t) multiresult;
	if(leftover < range ) {
		threshold = -range % range ;
		while (leftover < threshold) {
			random32bit =  pcg32_random();
			multiresult = random32bit * range;
			leftover = (uint32_t) multiresult;
		}
	}
	return multiresult >> 32; // [0, range)
}

static int fastrandombytes_int(const uint32_t range) noexcept {
	uint64_t random32bit, multiresult;
	uint32_t leftover;
	uint32_t threshold;
	random32bit =  pcg32_random();
	multiresult = random32bit * range;
	leftover = (uint32_t) multiresult;
	if(leftover < range ) {
		threshold = -range % range ;
		while (leftover < threshold) {
			random32bit =  pcg32_random();
			multiresult = random32bit * range;
			leftover = (uint32_t) multiresult;
		}
	}
	return multiresult >> 32; // [0, range)
}



/*
	Written in 2016-2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
	To the extent possible under law, the author has dedicated all copyright and
	related and neighboring rights to this software to the public domain worldwide.
	This software is distributed without any warranty.

	See <http://creativecommons.org/publicdomain/zero/1.0/>.
    This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than xoroshiro128**. It
   passes all tests we are aware of except for the four lower bits, which might
   fail linearity tests (and just those), so if low linear complexity is not
   considered an issue (as it is usually the case) it can be used to generate
   64-bit outputs, too; moreover, this generator has a very mild Hamming-weight
   dependency making our test (http://prng.di.unimi.it/hwd.php) fail after 5 TB
   of output; we believe this slight bias cannot affect any application. If you
   are concerned, use xoroshiro128++, xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and right
   shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have a
   64-bit seed, we suggest to seed a splitmix64 generator and use its output to
   fill s.

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly better
   results in our test than the 2016 version (a=55, b=14, c=36).
*/

static inline uint64_t rotl(const uint64_t x, const int k) noexcept {
	return (x << k) | (x >> (64 - k));
}

uint64_t xoroshiro128_random_uint64_t(uint64_t *S0, uint64_t *S1) {
	const uint64_t s0 = *S0;
	uint64_t s1 = *S1;
	const uint64_t result = s0 + s1;

	s1 ^= s0;
	*S0 = rotl(s0, 24) ^ s1 ^ (s1 << 16);  // a, b
	*S1 = rotl(s1, 37);                    // c

	return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64 non-overlapping
   subsequences for parallel computations. */
void xoroshiro128_jump(uint64_t *S0, uint64_t *S1) {
	static const uint64_t JUMP[] = {0xdf900294d8f554a5, 0x170865df4b3201fc};

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for (unsigned i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
		for (int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= *S0;
				s1 ^= *S1;
			}
			xoroshiro128_random_uint64_t(S0, S1);
		}

	*S0 = s0;
	*S1 = s1;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^96 calls to next(); it can be used to generate 2^32 starting points, from
   each of which jump() will generate 2^32 non-overlapping subsequences for
   parallel distributed computations. */
void xoroshiro128_long_jump(uint64_t *S0, uint64_t *S1) noexcept {
	static const uint64_t LONG_JUMP[] = {0xd2a98b26625eee7b, 0xdddf9b1090aa7ac1};

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for (unsigned i = 0; i < sizeof(LONG_JUMP) / sizeof(*LONG_JUMP); i++)
		for (int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= *S0;
				s1 ^= *S1;
			}

			xoroshiro128_random_uint64_t(S0, S1);
		}

	*S0 = s0;
	*S1 = s1;
}

int xoroshiro128_seed_random(uint64_t *S0, uint64_t *S1) noexcept {
	uint64_t new_s[2];
	FILE *urandom_fp;

	urandom_fp = fopen("/dev/urandom", "r");
	if (urandom_fp == NULL) return 0;
	if (fread(&new_s, 8, 2, urandom_fp) != 2) return 0;
	fclose(urandom_fp);

	*S0 = new_s[0];
	*S1 = new_s[1];

	return 1;
}

uint64_t xoroshiro128_random_lim(uint64_t limit, uint64_t *S0, uint64_t *S1) noexcept {
	uint64_t divisor = 0xffffffffffffffffUL / (limit + 1);
	uint64_t retval;

	do {
		retval = xoroshiro128_random_uint64_t(S0, S1) / divisor;
	} while (retval > limit);

	return retval;
}

/**
 * FINAL WRAPPER FUNCTIONS FOR SEEDING AND NUMBER GENERATION. HERE ONE CAN CHOOSE BETWEEN
 * xorshf96, pcg, and xoroshiro128.
*/

/*
 *  Final wrapper function. Can easily be replaced with whatever function you want.
 */
static int fastrandombytes(void *buf, size_t n) noexcept {
	return xorshf96_fastrandombytes(buf, n);
}

/// Wrapper function that seeds the chosen rng (Probably always xorshf96).
/// \param i
static void random_seed(uint64_t i) noexcept {
	xorshf96_random_seed(i);
}

/// Final wrapper function that generates a random 64bit integer with the chosen rng
/// (prabably always xorshf96).
/// \return
static uint64_t fastrandombytes_uint64() noexcept {
	return xorshf96_fastrandombytes_uint64();
}

/// Creates a weight w uint64_t.
/// \param w
/// \return
template<typename T>
static T fastrandombytes_weighted(const uint32_t w) noexcept {
	T ret = (1u << w) - 1u;
	for (uint32_t i = 0; i < w; ++i) {
		const size_t to_pos = fastrandombytes_uint64() % ((sizeof(T)*8) - i);
		const size_t from_pos = i;

		const T from_mask = 1u << from_pos;
		const T to_mask = 1u << to_pos;

		const T from_read = (ret & from_mask) >> from_pos;
		const T to_read = (ret & to_mask) >> to_pos;

		ret ^= (-from_read ^ ret) & (1ul << to_pos);
		ret ^= (-to_read ^ ret) & (1ul << from_pos);
	}

	return ret;
}

/// Generates random bits and returns them in the specified type T.
/// \tparam T - The type in which the bits are returned in.
/// \tparam bits - The number of bits to generate.
/// \return The integer-type in which the bits are stored.
template<typename T, const uint32_t bits>
static inline T fastrandombits() noexcept {
	constexpr uint32_t Tbits = sizeof(T)*8;
	constexpr T mask = (T(1u) << bits) - T(1u);
	static_assert(Tbits >= bits);

	return T(xorshf96_fastrandombytes_uint64()) & mask;
}
#endif //SMALLSECRETLWE_RANDOM_H
