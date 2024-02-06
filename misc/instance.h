#ifndef INCLUDE_INSTANCE
#define INCLUDE_INSTANCE

#include <iostream>
#include <cstdint>
#include "helper.h"
#include "m4ri/m4ri.h"
#include "custom_matrix.h"


// problem instance for the algorithm.
struct DecodingInstance {
public:
	uint32_t n, k;

	// CPU Parameters
	mzd_t *AT;          	// parity check matrix transposed: n \times n-k
	mzd_t *A;           	// parity check matrix: n-k \times n
	mzd_t *syndrome;     	// syndrome 1 \times n-k
	mzd_t *syndromeT;    	// syndrome: n-k \times 1
	mzd_t *error;    		// error 1 \times n
	const uint32_t iters = 2;

	/// print some information
	void print() {
		std::cout << "n: " << n
		          << ", k: " << k
				  << std::endl;
	}

	// checks if the set error vector results in the syndrome. Returns true if He=s, false otherwise.
	bool check() {
		mzd_t *ss_tmp = mzd_init(syndrome->ncols, syndrome->nrows);
		mzd_t *ss_tmp_T = mzd_init(syndrome->nrows, syndrome->ncols);
		mzd_t *ee_T = mzd_init(n, 1);
		ASSERT(ss_tmp);
		ASSERT(ss_tmp_T);
		ASSERT(ee_T);

		bool ret = true;
		mzd_transpose(ee_T, error);
		mzd_mul_naive(ss_tmp, A, ee_T);
		mzd_transpose(ss_tmp_T, ss_tmp);

		print_matrix("FINAL  e: ", error);
		print_matrix("SHOULD s: ", syndrome);
		print_matrix("IS     s:", ss_tmp_T);

		for (uint32_t i = 0; i < n-k; ++i) {
			if (mzd_read_bit(syndrome, 0, i) != mzd_read_bit(ss_tmp, i, 0)) {
				ret = false;
				break;
			}
		}

		mzd_free(ss_tmp);
		mzd_free(ss_tmp_T);
		mzd_free(ee_T);
		return ret;
	}

	// Constructor.
	DecodingInstance(
			const char *H,
			const char *S,
			const uint32_t n,
	        const uint32_t k) : n(n), k(k) {
		AT = mzd_from_str(n, n - k, H),
		syndrome = mzd_from_str(1, n - k, S),
	
		A = mzd_transpose(nullptr, AT);
		syndromeT = mzd_transpose(nullptr, syndrome);

		error = mzd_init(1, n);
	}

	// Destructor.
	~DecodingInstance() {
		mzd_free(syndrome);
		mzd_free(syndromeT);
		mzd_free(A);
		mzd_free(AT);
		mzd_free(error);
	}
};

#endif
