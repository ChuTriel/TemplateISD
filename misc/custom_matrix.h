#ifndef SMALLSECRETLWE_CUSTOM_MATRIX_H
#define SMALLSECRETLWE_CUSTOM_MATRIX_H

//#ifdef __APPLE__
//#include <boost/align/aligned_alloc.hpp>
//using boost::alignment::aligned_alloc;
//#endif

#include "m4ri/m4ri.h"
#include "m4ri/mzd.h"
#include "helper.h"
#include "random.h"

#include <type_traits>
#include <utility>
#include <array>
#include <memory>

//#include <omp.h>

#ifdef USE_AVX2
#include <immintrin.h>
#include <emmintrin.h>
#endif

#define MATRIX_AVX_PADDING(len) (((len + 255) / 256) * 256)

#define MAX_K 8
#define WORD_SIZE (8 * sizeof(word))



uint32_t hamming_weight_column(const mzd_t *A, const uint32_t col) {
	// Allow to calc the weight of the syndrome?
	ASSERT(uint32_t(A->ncols) >= col);

	uint32_t ret = 0;
	for (uint32_t i = 0; i < uint32_t(A->nrows); ++i) { ret += mzd_read_bit(A, i, col); }
	return ret;
}

// TODO, create a class out of this
typedef struct customMatrixData {
	uint32_t **rev;
	uint32_t **diff;

	// Unpadded number of columns
	uint32_t real_nr_cols;

	uint32_t working_nr_cols;

	alignas(64) uint64_t *lookup_table;
} customMatrixData;

typedef struct matrix_perm_t {
	// The swap operations in LAPACK format.
	uint16_t *values;

   	// The length of the swap array.
	size_t length;
} matrix_perm;

matrix_perm* matrix_perm_init(const size_t length) {
	matrix_perm *P = (matrix_perm *)aligned_alloc(64, sizeof(matrix_perm));
	P->values = (uint16_t *)aligned_alloc(64, sizeof(uint16_t) * length);
	P->length = length;

	if (P == nullptr || P->values == nullptr) {
		std::cout << "Alloc Error\n";
		exit(1);
	}

	for (uint32_t i = 0; i < length; ++i) {
		P->values[i] = i;
	}

	return P;
}

void inline matrix_random(matrix_perm *P) {
	for (uint32_t i = 0; i < uint32_t(P->length-1); ++i) {
		word pos = fastrandombytes_uint64() % (P->length - i);

		ASSERT(i+pos < uint64_t(P->length));
		std::swap(P->values[i], P->values[i+pos]);
	}
}

void inline mzp_random(mzp_t *P) {
	for (uint32_t i = 0; i < uint32_t(P->length-1); ++i) {
		word pos = fastrandombytes_uint64() % (P->length - i);

		ASSERT(i+pos < uint64_t(P->length));
		std::swap(P->values[i], P->values[i+pos]);
	}
}



// implementation taken from gray_code.c
// a = #rows, b = #cols
constexpr int matrix_opt_k(const int a, const int b) noexcept {
	return  MIN(__M4RI_MAXKAY, MAX(1, (int)(0.75 * (1 + const_log(MIN(a, b))))));
}

inline void xor_avx1_new(const uint8_t *x,
                         const uint8_t *y,
                         uint8_t *z,
                         unsigned nn) noexcept {
	for (uint32_t i = 0; i < nn; i += 1) {
#ifdef USE_AVX2
#ifdef USE_AVX2_SPECIAL_ALIGNMENT
		__m256i x_avx = _mm256_load_si256((__m256i *)x + i);
		__m256i y_avx = _mm256_load_si256((__m256i *)y + i);
		__m256i z_avx = _mm256_xor_si256(x_avx, y_avx);
		_mm256_store_si256((__m256i *)z + i, z_avx);
#else
		__m256 x_avx = _mm256_loadu_ps((float*)x + 8*i);
		__m256 y_avx = _mm256_loadu_ps((float*)y + 8*i);
		__m256 z_avx = _mm256_xor_ps(x_avx,y_avx);
		_mm256_storeu_ps((float*)z + 8*i, z_avx);
#endif
#else
		((uint64_t *)z)[4*i] = ((uint64_t *)x)[4*i] ^ ((uint64_t *)y)[4*i];
		((uint64_t *)z)[4*i+1] = ((uint64_t *)x)[4*i+1] ^ ((uint64_t *)y)[4*i+1];
		((uint64_t *)z)[4*i+2] = ((uint64_t *)x)[4*i+2] ^ ((uint64_t *)y)[4*i+2];
		((uint64_t *)z)[4*i+3] = ((uint64_t *)x)[4*i+3] ^ ((uint64_t *)y)[4*i+3];
#endif
	}
}
// TODO AVX2
void mzd_row_xor(mzd_t *out,
                 const rci_t i,
                 const rci_t j) noexcept {
	ASSERT(out->nrows > i && out->nrows > j);
	for (uint32_t l = 0; l < uint32_t(out->width); ++l) {
		out->rows[i][l] ^= out->rows[j][l];
	}
}

// Is this really smart?
// this changes only the pointer not the content
inline void matrix_swap_rows_new(mzd_t *__restrict__ M,
                                 const size_t i,
                                 const size_t j) noexcept {
	std::swap(M->rows[i], M->rows[j]);
}

// in comparison to `matrix_swap_rows_new` this changes the content and not just the pointesr
inline void matrix_swap_rows_new2(mzd_t *__restrict__ M,
                                  const size_t i,
                                  const size_t j) noexcept {
	const size_t cols = M->ncols;
	const size_t cols_padded_ymm = MATRIX_AVX_PADDING(cols) / 256;

	xor_avx1_new((uint8_t *) M->rows[i], (uint8_t *) M->rows[j], (uint8_t *) M->rows[i], cols_padded_ymm);
	xor_avx1_new((uint8_t *) M->rows[j], (uint8_t *) M->rows[i], (uint8_t *) M->rows[j], cols_padded_ymm);
	xor_avx1_new((uint8_t *) M->rows[i], (uint8_t *) M->rows[j], (uint8_t *) M->rows[i], cols_padded_ymm);
}

/// custom function to append in_matrix to out_matrix.
/// \param out_matrix
/// \param in_matrix
/// \param start_r
/// \param start_c
void mzd_append(mzd_t *out_matrix, const mzd_t *const in_matrix, const int start_r, const int start_c) noexcept {
	assert(out_matrix->nrows > start_r);
	assert(out_matrix->ncols > start_c);
	assert(out_matrix->nrows - start_r == in_matrix->nrows);
	for (int row = 0; row < in_matrix->nrows; ++row) {
		for (int col = 0; col < in_matrix->ncols; ++col) {
			auto bit = mzd_read_bit(in_matrix, row, col);
			mzd_write_bit(out_matrix, row + start_r, col + start_c, bit);
		}
	}
}


mzd_t *_mzd_transpose(mzd_t *DST, mzd_t const *A);

mzd_t *matrix_transpose(mzd_t *DST, mzd_t const *A) noexcept {
	if (DST == NULL) {
		DST = mzd_init(A->ncols, A->nrows);
	} else if (__M4RI_UNLIKELY(DST->nrows < A->ncols || DST->ncols < A->nrows)) {
		m4ri_die("mzd_transpose: Wrong size for return matrix.\n");
	} else {
		/** it seems this is taken care of in the subroutines, re-enable if running into problems **/
		// mzd_set_ui(DST,0);
	}

	if (A->nrows == 0 || A->ncols == 0) return mzd_copy(DST, A);

	if (__M4RI_LIKELY(!mzd_is_windowed(DST) && !mzd_is_windowed(A)))
		return _mzd_transpose(DST, A);
	int A_windowed = mzd_is_windowed(A);
	if (A_windowed) A = mzd_copy(NULL, A);
	if (__M4RI_LIKELY(!mzd_is_windowed(DST)))
		_mzd_transpose(DST, A);
	else {
		mzd_t *D = mzd_init(DST->nrows, DST->ncols);
		_mzd_transpose(D, A);
		mzd_copy(DST, D);
		mzd_free(D);
	}
	if (A_windowed) mzd_free((mzd_t *)A);
	return DST;
}

// The same as `mzd_concat` but with the oddity that a new matrix with avx padding is generated if necessary.
// This is sometimes called augment
mzd_t *matrix_concat(mzd_t *C, mzd_t const *A, mzd_t const *B) noexcept {
	if (A->nrows != B->nrows) { m4ri_die("mzd_concat: Bad arguments to concat!\n"); }

	if (C == NULL) {
		C = mzd_init(A->nrows, MATRIX_AVX_PADDING(A->ncols + B->ncols));
		// This is the important point, we `resize` the ouput matrix to its normal expected size.
		C->ncols = A->ncols + B->ncols;
	} else if (C->nrows != A->nrows || C->ncols < (A->ncols + B->ncols)) {
		m4ri_die("mzd_concat: C has wrong dimension!\n");
	}

	for (rci_t i = 0; i < A->nrows; ++i) {
		word *dst_truerow = C->rows[i];
		word *src_truerow = A->rows[i];
		for (wi_t j = 0; j < A->width; ++j) { dst_truerow[j] = src_truerow[j]; }
	}

	for (rci_t i = 0; i < B->nrows; ++i) {
		for (rci_t j = 0; j < B->ncols; ++j) {
			mzd_write_bit(C, i, j + A->ncols, mzd_read_bit(B, i, j));
		}
	}

	__M4RI_DD_MZD(C);
	return C;
}

// generates the DOOM shits
// in must be the syndrom in column form
// writes the shifts directly into the out matrix.
mzd_t *matrix_down_shift_into_matrix(mzd_t *out, const mzd_t *in, const uint32_t col, const uint32_t i) noexcept {
	if ((in->ncols != 1) || (in->nrows != out->nrows)) m4ri_die("matrix_shift: Bad argument!\n");

	if (out == nullptr) {
		out = mzd_init(in->nrows, in->nrows);
	}

	for (int j = 0; j < in->nrows; ++j) {
		mzd_write_bit(out, (j + i) % in->nrows, col, mzd_read_bit(in, j, 0));
	}

	return out;
}

// generates the DOOM shits
// in must be the syndrom in column form
// writes the shifts directly into the out matrix.
mzd_t *matrix_up_shift_into_matrix(mzd_t *out, const mzd_t *in, const uint32_t col, const uint32_t i) noexcept {
	if (in->ncols != 1) m4ri_die("matrix_shift: Bad argument!\n");

	if (out == nullptr) {
		out = mzd_init(in->nrows, in->nrows);
	}

	for (int j = 0; j < in->nrows; ++j) {
		mzd_write_bit(out, j, col, mzd_read_bit(in, (j + i) % in->nrows, 0));
	}

	return out;
}

// generates the DOOM shits
// in must be the syndrom in column form
mzd_t *matrix_down_shift(mzd_t *out, const mzd_t *in, uint32_t i) noexcept {
	if (in->ncols != 1) m4ri_die("matrix_shift: Bad argument!\n");

	if (out == nullptr) {
		out = mzd_init(in->nrows, in->ncols);
	}

	for (int j = 0; j < in->nrows; ++j) {
		mzd_write_bit(out, (j + i) % in->nrows, 0, mzd_read_bit(in, j, 0));
	}

	return out;
}

// generates the DOOM shits
// in must be the syndrom in column form
mzd_t *matrix_up_shift(mzd_t *out, const mzd_t *in, uint32_t i) noexcept {
	if (in->ncols != 1) m4ri_die("matrix_shift: Bad argument!\n");

	if (out == nullptr) {
		out = mzd_init(in->nrows, in->ncols);
	}

	for (int j = 0; j < in->nrows; ++j) {
		mzd_write_bit(out, j, 0, mzd_read_bit(in, (j + i) % in->nrows, 0));
	}

	return out;
}

// create A with `mzp_t *A = mzp_init(length)
// input permutation `A` should have initilaised with: `for(int i = 0; i < A->length; i++) A->value[i] = i` or something.
void matrix_create_random_permutation(mzd_t *A, mzp_t *P) noexcept {
	mzd_t * AT= mzd_init(A->ncols,A->nrows);
	mzd_transpose(AT, A);
	// dont permute the last column since it is the syndrome
	for (uint32_t i = 0; i < uint32_t(P->length-1); ++i) {
		word pos = fastrandombytes_uint64() % (P->length - i);

		ASSERT(i+pos < uint64_t(P->length));
		std::swap(P->values[i], P->values[i+pos]);

		mzd_row_swap(AT, i, i+pos);
	}
	mzd_transpose(A, AT);
	mzd_free(AT);
}

// optimized version to which you have additionally pass the transposed, So its not created/freed every time
void matrix_create_random_permutation(mzd_t *__restrict__ A,
                                      mzd_t *__restrict__ AT,
                                      mzp_t *__restrict__ P) noexcept {
	matrix_transpose(AT, A);

	// dont permute the last column since it is the syndrome
	for (uint32_t i = 0; i < uint32_t(P->length-1); ++i) {
		word pos = fastrandombytes_uint64() % (P->length - i);

		ASSERT(i+pos < uint32_t(P->length));
		std::swap(P->values[i], P->values[i+pos]);
		mzd_row_swap(AT, i, i+pos);
		//matrix_swap_rows_new2(AT, i, i+pos);
	}
	matrix_transpose(A, AT);
}


static size_t matrix_gauss_submatrix(mzd_t *__restrict__ M,
                              const size_t r,
                              const size_t c,
                              const size_t rows,
                              const size_t cols,
                              const size_t k) noexcept {
	size_t start_row = r, j;
	const size_t cols_padded_ymm = MATRIX_AVX_PADDING(cols) / 256;

	for (j = c; j < c + k; ++j) {
		int found = 0;
		for (size_t i = start_row; i < rows; ++i) {
			for (size_t l = 0; l < j - c; ++l) {
				if ((M->rows[i][(c + l) / WORD_SIZE] >> ((c + l) % WORD_SIZE)) & 1) {
					xor_avx1_new((uint8_t *) M->rows[r + l],
					             (uint8_t *) M->rows[i],
					             (uint8_t *) M->rows[i],
					             cols_padded_ymm);
				}
			}

			if ((M->rows[i][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
				matrix_swap_rows_new(M, i, start_row);
				for (size_t l = r; l < start_row; ++l) {
					if ((M->rows[l][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
						xor_avx1_new((uint8_t *) M->rows[start_row],
						             (uint8_t *) M->rows[l],
						             (uint8_t *) M->rows[l],
						             cols_padded_ymm);
					}
				}
				++start_row;
				found = 1;
				break;
			}
		}

		if (found == 0) {
			break;
		}
	}
	return j - c;
}

// specialised gaus which looks ahead in the l window
inline size_t matrix_gauss_submatrix_opt(mzd_t *__restrict__ M,
                                         const size_t r,
                                         const size_t c,
                                         const size_t rows,
                                         const size_t cols,
                                         const size_t rstop,
                                         const size_t lwin_start,
                                         const size_t k) noexcept {
	size_t start_row = r, j;
	const size_t cols_padded_ymm = MATRIX_AVX_PADDING(cols) / 256;

	// for each column
	for (j = c; j < c + k; ++j) {
		int found = 0;
		for (size_t i = start_row; i < rstop; ++i) {
			for (size_t l = 0; l < j - c; ++l) {
				if ((M->rows[i][(c + l) / WORD_SIZE] >> ((c + l) % WORD_SIZE)) & 1) {
					xor_avx1_new((uint8_t *) M->rows[r + l],
					             (uint8_t *) M->rows[i],
					             (uint8_t *) M->rows[i],
					             cols_padded_ymm);
				}
			}

			if ((M->rows[i][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
				matrix_swap_rows_new(M, i, start_row);
				for (size_t l = r; l < start_row; ++l) {
					if ((M->rows[l][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
						xor_avx1_new((uint8_t *) M->rows[start_row],
						             (uint8_t *) M->rows[l],
						             (uint8_t *) M->rows[l],
						             cols_padded_ymm);
					}
				}

				++start_row;
				found = 1;
				break;
			}
		}

		// second try: find the pivot element in the l window
		if (found == 0) {
			for (size_t i = lwin_start; i < rows; ++i) {
				for (size_t l = 0; l < j - c; ++l) {
					if ((M->rows[i][(c + l) / WORD_SIZE] >> ((c + l) % WORD_SIZE)) & 1) {
						xor_avx1_new((uint8_t *) M->rows[r + l],
						             (uint8_t *) M->rows[i],
						             (uint8_t *) M->rows[i],
						             cols_padded_ymm);
					}
				}

				if ((M->rows[i][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
					matrix_swap_rows_new(M, i, start_row);
					for (size_t l = r; l < start_row; ++l) {
						if ((M->rows[l][j / WORD_SIZE] >> (j % WORD_SIZE)) & 1) {
							xor_avx1_new((uint8_t *) M->rows[start_row],
							             (uint8_t *) M->rows[l],
							             (uint8_t *) M->rows[l],
							             cols_padded_ymm);
						}
					}

					++start_row;
					found = 1;
					break;
				}
			}
		}

		// okk we really couldnt find a pivot element
		if (found == 0) {
			break;
		}
	}

	return j - c;
}

void matrix_make_table(mzd_t *__restrict__ M,
                       const size_t r,
                       const size_t cols,
                       const size_t k,
                       uint64_t *T,
                       uint32_t **diff) noexcept {
	const size_t cols_padded = MATRIX_AVX_PADDING(cols);
	const size_t cols_padded_word = cols_padded / 64;
	const size_t cols_padded_ymm = cols_padded / 256;

	for (size_t i = 0; i < cols_padded_word; ++i) {
		T[i] = 0L;
	}

	for (size_t i = 0; i + 1 < 1UL << k; ++i) {
		xor_avx1_new((uint8_t *)M->rows[r + diff[k][i]],
		             (uint8_t *)T,
		             (uint8_t *)(T + cols_padded_word),
		             cols_padded_ymm);
		T += cols_padded_word;
	}
}


uint64_t matrix_read_bits(mzd_t *M,
                      const size_t x,
                      const size_t y,
                      const size_t nn) noexcept {
	const uint32_t spot  = y % WORD_SIZE;
	const uint64_t block = y / WORD_SIZE;

	// this must be negative...
	const int32_t spill = spot + nn - WORD_SIZE;
	uint64_t temp = (spill <= 0) ? M->rows[x][block] << -spill
	                         : (M->rows[x][block + 1] << (WORD_SIZE - spill)) |
	                           (M->rows[x][block] >> spill);
	return temp >> (WORD_SIZE - nn);
}


void matrix_process_rows(mzd_t *__restrict__ M,
                         const size_t rstart,
                         const size_t cstart,
                         const size_t rstop,
                         const size_t k,
                         const size_t cols,
                         uint64_t *T,
                         uint32_t **rev) noexcept {
	const size_t cols_padded = MATRIX_AVX_PADDING(cols);
	const size_t cols_padded_ymm = cols_padded / 256;
	const size_t cols_padded_word = cols_padded / 64;

	for (size_t r = rstart; r < rstop; ++r) {
		size_t x0 = rev[k][matrix_read_bits(M, r, cstart, k)];
		if (x0) {
			xor_avx1_new((uint8_t *) (T + x0 * cols_padded_word),   // in
			             (uint8_t *) M->rows[r],                    // in
			             (uint8_t *) M->rows[r],                    // out
			             cols_padded_ymm);
		}
	}
}

///
/// \param M 			matrix
/// \param k 			m4ri k
/// \param rstop 		row stop
/// \param matrix_data 	helper data
/// \param cstart 		column start.
/// \return
size_t matrix_echelonize_partial(mzd_t *M,
                                 const size_t k,
                                 const size_t rstop,
                                 customMatrixData *matrix_data,
                                 const size_t cstart) noexcept {
	alignas(64) uint32_t **rev      = matrix_data->rev;
	alignas(64) uint32_t **diff     = matrix_data->diff;
	alignas(64) uint64_t *xor_rows  = matrix_data->lookup_table;

	const size_t rows = M->nrows;
	const size_t cols = matrix_data->working_nr_cols;

	size_t kk = k;

	// current row
	size_t r = 0;
	// current column
	size_t c = cstart;

	while (c < rstop) {
		if (c + kk > rstop) {
			kk = rstop - c;
		}

		size_t kbar = matrix_gauss_submatrix(M, r, c, rows, cols, kk);
		if (kk != kbar)
			break;

		if (kbar > 0) {
			matrix_make_table  (M, r, cols, kbar, xor_rows, diff);
			// fix everything below
			matrix_process_rows(M, r + kbar, c, rows, kbar, cols, xor_rows, rev);
			// fix everything over it
			matrix_process_rows(M, 0, c, r, kbar, cols, xor_rows, rev);
		}

		r += kbar;
		c += kbar;
	}

	return r;
}

// specialised m4ri which is allowd to perform special look aheads
// in the l window
size_t matrix_echelonize_partial_opt(mzd_t *M,
                                 const size_t k,
                                 const size_t rstop,
                                 const size_t lwin_start,
                                 customMatrixData *matrix_data) noexcept {
	alignas(64) uint32_t **rev      = matrix_data->rev;
	alignas(64) uint32_t **diff     = matrix_data->diff;
	alignas(64) uint64_t *xor_rows  = matrix_data->lookup_table;

	const size_t rows = M->nrows;
	const size_t cols = matrix_data->working_nr_cols;

	ASSERT(lwin_start <= rows);
	ASSERT(rstop <= rows);

	size_t kk = k;

	// current row
	size_t r = 0;
	// current column
	size_t c = 0;

	while (c < rstop) {
		if (c + kk > rstop) {
			kk = rstop - c;
		}

		size_t kbar = matrix_gauss_submatrix_opt(M, r, c, rows, cols, rstop, lwin_start, kk);
		if (kk != kbar)
			break;

		if (kbar > 0) {
			matrix_make_table  (M, r, cols, kbar, xor_rows, diff);
			// fix everything below
			matrix_process_rows(M, r + kbar, c, rows, kbar, cols, xor_rows, rev);
			// fix everything over it
			matrix_process_rows(M, 0, c, r, kbar, cols, xor_rows, rev);
		}

		r += kbar;
		c += kbar;
	}

	return r;
}

// additionally, to the m4ri algorithm a fix is applied if m4ri fails
/// \param M
/// \param k			m4ri k
/// \param rstop
/// \param matrix_data
/// \param cstart
/// \param fix_col 		how many columns must be solved at the end of this function. Must be != 0
/// \param look_ahead 	how many the algorithm must start ahead of `fix_col` to search for a pivo element. Can be zero
/// \param permutation
/// \return
size_t matrix_echelonize_partial_plusfix(mzd_t *__restrict__ M,
                                         const size_t k,
                                         const size_t rstop,
										 customMatrixData *__restrict__ matrix_data,
                                         const size_t cstart,
										 const size_t fix_col,
                                         const size_t look_ahead,
                                         mzp_t *__restrict__ permutation) noexcept {
	size_t rang = matrix_echelonize_partial(M, k, rstop, matrix_data, cstart);
	// printf("%d %d\n", rang, rstop);

	for (size_t b = rang; b < rstop; ++b) {
		bool found = false;
		// find a column where in the last row is a one
		for (size_t i = fix_col+look_ahead; i < size_t(M->ncols); ++i) {
			if (mzd_read_bit(M, b, i)) {
				found = true;
				// if (i == b)
				// 	break;

				std::swap(permutation->values[i-look_ahead], permutation->values[b]);
				mzd_col_swap(M, b, i);
				break;
			}
		}

		// if something found, fix this row by adding it to each row where a 1 one.
		if (found) {
			for (size_t i = 0; i < b; ++i) {
				if (mzd_read_bit(M, i, b)) {
					mzd_row_xor(M, i, b);
				}
			}

			// and solve it below
			for (size_t i = b+1; i < M->nrows; ++i) {
				if (mzd_read_bit(M, i, b)) {
					mzd_row_xor(M, i, b);
				}
			}
		} else {
			// if we were not able to fix the gaussian elimination, we return the original rang which was solved
			return rang;
		}
	}

	// return the full rang, if we archived it.
	// return fix_col;
	return rstop;
}




size_t matrix_echelonize_partial(mzd_t *M,
                                 const size_t k,
                                 const size_t rstop,
                                 customMatrixData *matrix_data) noexcept {
	return matrix_echelonize_partial(M, k, rstop, matrix_data, 0);
}

static int gray_new(const uint32_t i,
                    const uint32_t k) noexcept {
	int lastbit = 0;
	int res = 0;
	for (int j = k; j-- > 0;) {
		int bit = i & (1 << j);
		res |= (lastbit >> 1) ^ bit;
		lastbit = bit;
	}
	return res;
}


///
void matrix_alloc_gray_code(uint32_t ***__restrict__ rev, uint32_t ***__restrict__ diff) noexcept {
	constexpr uint32_t alignment = 1024;
	*rev  = (uint32_t **) aligned_alloc(alignment, MAX(alignment, (MAX_K + 1) * sizeof(uint32_t *)));
	*diff = (uint32_t **) aligned_alloc(alignment, MAX(alignment, (MAX_K + 1) * sizeof(uint32_t *)));

	for (size_t k = 0; k <= MAX_K; ++k) {
		// (*rev)[k]  = (uint32_t *) aligned_alloc(alignment, (1 << k) * sizeof(uint32_t));
		// (*diff)[k] = (uint32_t *) aligned_alloc(alignment, (1 << k) * sizeof(uint32_t));

		(*rev)[k]  = (uint32_t *) malloc((1 << k) * sizeof(uint32_t));
		(*diff)[k] = (uint32_t *) malloc((1 << k) * sizeof(uint32_t));
	}
}


void matrix_free_gray_code(uint32_t **rev, uint32_t **diff) noexcept {
	for (size_t k = 0; k <= MAX_K; ++k) {
		free(rev[k]);
		free(diff[k]);
	}

	free(rev);
	free(diff);
}


void matrix_build_gray_code(uint32_t **rev, uint32_t **diff) noexcept {
	for (size_t k = 0; k <= MAX_K; ++k) {
		for (size_t i = 0; i < 1UL << k; ++i) {
			rev[k][gray_new(i, k)] = i;
		}

		for (size_t i = k + 1; i-- > 0;) {
			for (size_t j = 1; j < (1UL << i) + 1; ++j) {
				diff[k][j * (1 << (k - i)) - 1] = k - i;
			}
		}
	}
}

// special copy, which can be useful if you working with avx alignment
mzd_t *matrix_copy(mzd_t *N, mzd_t const *P) noexcept {
	if (N == P) {
		return N;
	}

	if (N == NULL) {
		N = mzd_init(P->nrows, P->ncols);
	}
	word *p_truerow, *n_truerow;
	wi_t const wide = P->width - 1;
	word mask_end   = P->high_bitmask;
	for (rci_t i = 0; i < P->nrows; ++i) {
		p_truerow = P->rows[i];
		n_truerow = N->rows[i];
		for (wi_t j = 0; j < wide; ++j) n_truerow[j] = p_truerow[j];
		n_truerow[wide] = (n_truerow[wide] & ~mask_end) | (p_truerow[wide] & mask_end);
	}
	__M4RI_DD_MZD(N);
	return N;
}

customMatrixData* init_matrix_data(int nr_columns) noexcept {
	// this should be alignend
	customMatrixData *matrix_data = (customMatrixData *) malloc(sizeof(customMatrixData));

	matrix_data->real_nr_cols = nr_columns;
	matrix_data->working_nr_cols = nr_columns;  // this can be overwritten

	matrix_alloc_gray_code(&matrix_data->rev, &matrix_data->diff);
	matrix_build_gray_code(matrix_data->rev, matrix_data->diff);

	matrix_data->lookup_table = (uint64_t *)aligned_alloc(4096, (MATRIX_AVX_PADDING(nr_columns) / 8) * (1<<MAX_K));
	return matrix_data;
}

void free_matrix_data(customMatrixData* matrix_data) noexcept {
	matrix_free_gray_code(matrix_data->rev, matrix_data->diff);
	free(matrix_data->lookup_table);
	free(matrix_data);
}

mzd_t *matrix_init(rci_t r, rci_t c) noexcept {
	auto padding = MATRIX_AVX_PADDING(c);
	return mzd_init(r, padding);
}

/// ein versuch die matrix so zu splitten, dass die H matrix 256 alignent ist.
mzd_t *matrix_init_split(const mzd_t *A, const mzd_t *s, const uint32_t nkl, const uint32_t c) noexcept {
	assert(s->nrows==1);

	auto padding = MATRIX_AVX_PADDING(A->ncols + c);
	mzd_t *r = mzd_init(A->nrows, padding);

	for (uint32_t row = 0; row < uint32_t(A->nrows); row++) {
		for (uint col = 0; col < nkl; ++col) {
			mzd_write_bit(r, row, col, mzd_read_bit(A, row, col));
		}

		for (int col = nkl; col < A->ncols; ++col) {
			mzd_write_bit(r, row, col+c, mzd_read_bit(A, row, col));
		}

		mzd_write_bit(r, row, A->ncols+c, mzd_read_bit(s, 0, row));

	}
	return r;
}

// the same as above.
static void print_matrix(const char *s, const mzd_t *A, int start_row=-1, int end_row=-1) noexcept {
	const int ss = start_row == -1 ? 0 : start_row;
	const int ee = end_row == -1 ? A->nrows : end_row;

	if (s != nullptr)
		printf("%s\n", s);
	for (int i = ss; i < ee; ++i) {
		mzd_fprint_row(stdout, A, i);
	}
	//printf("\n");
}

static void print_matrix(const char *s, const mzd_t *A, int start_row, int end_row, int start_col, int end_col) noexcept {
	const int sstart_row = start_row == -1 ? 0 : start_row;
	const int eend_row = end_row == -1 ? A->nrows : end_row;

	const int sstart_col = start_col == -1 ? 0 : start_col;
	const int eend_col = end_col == -1 ? A->ncols : end_col;

	mzd_t *B = mzd_init(eend_row-sstart_row, eend_col-sstart_col);
	mzd_submatrix(B, A, sstart_row, sstart_col, eend_row, eend_col);

	printf("%s\n", s);
	for (int i = 0; i < eend_row-sstart_row; ++i) {
		mzd_fprint_row(stdout, B, i);
	}
	printf("\n");
	mzd_free(B);
}


template<typename T, std::size_t N, std::size_t... I>
constexpr auto create_array_impl(std::index_sequence<I...>) {
	return std::array<T, N>{ {I...} };
}

template<const uint32_t n, const uint32_t k, const uint32_t l>
size_t matrix_echelonize_partial_plusfix_opt (
		mzd_t *__restrict__ A,
		mzd_t *__restrict__ AT,
		const size_t m4ri_k,
        const size_t max_column,
		customMatrixData *__restrict__ matrix_data,
		mzp_t *__restrict__ P) noexcept {
	constexpr uint32_t nkl  = n-k-l;
	constexpr uint16_t zero = uint16_t(-1);

#if 1
	// create the transposed
	mzd_transpose(AT, A);

	// constexpr auto unitpos = create_array_impl<uint16_t, n-k-l>(std::make_index_sequence<n-k-l>{});
	// std::vector<uint16_t> unitpos(nkl, zero);
	// std::vector<uint16_t> posfix(n, zero);
	//auto posfix = create_array_impl<uint16_t, n>(std::make_index_sequence<n>{});

	//uint16_t unitpos[nkl] = {[0 ... nkl-1] = zero};
	//uint16_t posfix [n] = {[0 ... n-1] = zero};

	static uint16_t unitpos[nkl];
	static uint16_t posfix [n];
	for (uint32_t i = 0; i < nkl; ++i) {
		unitpos[i] = zero;
	}

	for (uint32_t i = 0; i < nkl; ++i) {
		posfix[i] = i;
	}

	for (uint32_t i = nkl; i < n; ++i) {
		posfix[i] = zero;
	}
	uint32_t unitctr = 0, switchctr = 0;

	// dont permute the last column since it is the syndrome
	for (uint32_t i = 0; i < uint32_t(P->length-1); ++i) {
		uint64_t pos = fastrandombytes_uint64() % (P->length - i);
		ASSERT(i+pos < uint64_t(P->length));

		std::swap(posfix[i], posfix[i+pos]);
		std::swap(P->values[i], P->values[pos+i]);
		// matrix_swap_rows_new2(AT, i, i+pos);
		mzd_row_swap(AT, i, i+pos);

		if (posfix[i] != zero && (i < n-k-l)) {
			unitpos[i] = posfix[i];
			unitctr += 1;
		}
	}

	// move in the matrix everything to the right
	for (uint32_t i = nkl; i > 0; i--) {
		// move the right pointer to the next free place from the right
		if (unitpos[i-1] == zero) {
			// move the left pointer to the next element
			while (unitpos[switchctr] == zero)
				switchctr += 1;

			if (i-1 <= switchctr)
				break;

			std::swap(unitpos[i-1], unitpos[switchctr]);
			std::swap(P->values[i-1], P->values[switchctr]);

			// matrix_swap_rows_new2(AT, i-1, switchctr);
			mzd_row_swap(AT, i-1, switchctr);
			switchctr += 1;
		}
	}

	// for (int i = 0; i < n-k-l; ++i) {
	// 	std::cout << i << ":" << unitpos[n-k-l-1-i] << ":" << n-k-l-1-i << "\n";
	// }
	// std::cout << "\n" << switchctr << ":" << unitctr << "\n";
	// std::cout << std::flush;

	// transpose it back
	mzd_transpose(A, AT);

	// mzd_print(A);
	// std::cout << "\n";

	// reset the array
	for (uint32_t i = 0; i < n-k-l; ++i) {
		posfix[i] = i;
	}

	// swap everything down
	for (uint32_t i = 0; i < unitctr; i++){
		// matrix_swap_rows_new2(A, n-k-l-1-i, posfix[unitpos[n-k-l-1-i]]);
		mzd_row_swap(A, nkl-1-i, posfix[unitpos[nkl-1-i]]);
		std::swap(posfix[nkl-1-i], posfix[unitpos[nkl-1-i]]);
	}

	// mzd_print(A);
	// std::cout << "\n";

	// final fix
	for (uint32_t i = nkl-unitctr; i < nkl; ++i) {
		if (mzd_read_bit(A, i, i) == 0) {
			for (uint32_t j = 0; j < nkl; ++j) {
				if (mzd_read_bit(A, j, i)) {
					// matrix_swap_rows_new2(A, i, j);
					mzd_row_swap(A, i, j);
					break;
				}
			}
		}
	}

	// mzd_print(A);
	// std::cout << "\n";

	// apply the final gaussian elimination on the first coordinates
	matrix_echelonize_partial_opt(A, m4ri_k, nkl-unitctr, nkl, matrix_data);
	// std::cout << unitctr << " " << nkl-unitctr << "\n";
	return nkl;
#else
	return matrix_echelonize_partial_plusfix(A, m4ri_k, n - k - l, matrix_data, 0, n - k - l, 0, P);
#endif
}

template<const uint32_t n,
	const uint32_t k,
	const uint32_t l,
	const uint32_t c>
size_t matrix_echelonize_partial_plusfix_opt_onlyc(
		mzd_t *__restrict__ A,
		mzd_t *__restrict__ AT,
		const size_t m4ri_k,
        const size_t max_column,
		customMatrixData *__restrict__ matrix_data,
		mzp_t *__restrict__ P) noexcept {
	constexpr uint32_t nkl  = n-k-l;
	constexpr uint16_t zero = uint16_t(-1);
	// constexpr uint32_t mod = k+l;

	// some security checks
	static_assert(nkl >= c);
	ASSERT(c > m4ri_k);

	// create the transposed
	mzd_transpose(AT, A);

	//std::vector<uint16_t> unitpos1(c, zero);
	//std::vector<uint16_t> unitpos2(c, zero);
	//std::vector<uint16_t> posfix  (nkl, zero);
	//std::vector<uint16_t> perm    (n, zero);
	static uint16_t unitpos1[c];
	static uint16_t unitpos2[c];
	static uint16_t posfix[nkl];
	static uint16_t perm[n];

	#pragma unroll
	for (uint32_t i = 0; i < c; ++i) {
		unitpos1[i] = zero;
	}

	#pragma unroll
	for (uint32_t i = 0; i < c; ++i) {
		unitpos2[i] = zero;
	}

	for (uint32_t i = 0; i < nkl; ++i) {
		posfix[i] = i;
	}

	for (uint32_t i = 0; i < n; ++i) {
		perm[i] = i;
	}

	for (uint32_t i = 0; i < c; ++i) {
		const uint32_t pos1 = fastrandombytes_uint64() % (k+l-i);
		const uint32_t pos2 = fastrandombytes_uint64() % (n-k-l-i);
		std::swap(perm[n-k-l+i], perm[n-k-l+i+pos1]);
		std::swap(perm[i], perm[i+pos2]);

	}

	// dont permute the last column since it is the syndrome
	for (uint32_t i = 0; i < c; ++i) {
		const uint32_t pos1 = perm[i];
		const uint32_t pos2 = perm[n-k-l+i];
		ASSERT(pos1 < P->length);
		ASSERT(pos2 < P->length);

		std::swap(P->values[pos1], P->values[pos2]);
		mzd_row_swap(AT, pos1, pos2);

		unitpos1[i] = pos1;
		if (pos1 < c) {
			unitpos2[pos1] = 0;
		}
	}

	// move in the matrix everything to the left
	uint32_t left_ctr, right_ctr = 0;//, switch_ctr = 0;
	for (left_ctr = 0; left_ctr < c && right_ctr < c; left_ctr++) {
		if (unitpos2[left_ctr] == 0)
			continue;

		while(right_ctr < (c-1) && unitpos1[right_ctr] < c) {
			right_ctr += 1;
		}

		//if (right_ctr >= c || unitpos1[right_ctr] >= (n-k-l)) {
		//	std::cout << "kekw: ";
		//	std::cout << omp_get_thread_num() << ":";
		//	std::cout << right_ctr << "\n";
		//	for (uint32_t i = 0; i < c; i++){
		//		std::cout << unitpos1[i] << " ";
		//	}
		//}
		//ASSERT(right_ctr < c);
		//ASSERT(unitpos1[right_ctr] < n-k-l);

		mzd_row_swap(AT, left_ctr, unitpos1[right_ctr]);
		std::swap(P->values[left_ctr], P->values[unitpos1[right_ctr]]);
		posfix[unitpos1[right_ctr]] = left_ctr;

		// step to the next element
		right_ctr += 1;
	}

	// transpose it back
	mzd_transpose(A, AT);

	for (uint32_t i = nkl; i > c; i--){
		mzd_row_swap(A, i-1, posfix[i-1]);
		std::swap(posfix[i-1], posfix[posfix[i-1]]);
	}

	// apply the final gaussian elimination on the first coordinates
	const uint32_t rang = matrix_echelonize_partial_opt(A, m4ri_k, c, nkl, matrix_data);

	// fix the first
	for (size_t b = rang; b < c; ++b) {
		bool found = false;
		// find a column where in the last row is a one
		for (size_t i = n-k-l; i < n; ++i) {
			if (mzd_read_bit(A, b, i)) {
				found = true;

				std::swap(P->values[i], P->values[b]);
				mzd_col_swap(A, b, i);
				break;
			}
		}

		ASSERT(found);

		// if something found, fix this row by adding it to each row where a 1 one.
		if (found) {
			for (size_t i = 0; i < b; ++i) {
				if (mzd_read_bit(A, i, b)) {
					mzd_row_xor(A, i, b);
				}
			}

			// and solve it below
			for (size_t i = b+1; i < n-k; ++i) {
				if (mzd_read_bit(A, i, b)) {
					mzd_row_xor(A, i, b);
				}
			}
		}
	}


	// std::cout << unitctr << " " << nkl-unitctr << "\n";

	// mzd_print(A);
	// std::cout << "\n";
	return nkl;
}
#endif //SMALLSECRETLWE_CUSTOM_MATRIX_H
