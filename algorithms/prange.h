#ifndef DECODING_PRANGE_H
#define DECODING_PRANGE_H

#include "custom_matrix.h"
#include "instance.h"
#include <iostream>
#include <stdint.h>
#include <fstream>
#include <string>
#include <atomic>
#include <chrono>


class ConfigPrange 
{
public:

    const uint32_t n, k, w;
    const uint32_t m4ri_k;


	constexpr ConfigPrange(const uint32_t n,
	                              const uint32_t k,
	                              const uint32_t w) :
                                  n(n),
                                  k(k),
                                  w(w),
                                  m4ri_k(matrix_opt_k(n - k, MATRIX_AVX_PADDING(n)))
	{}

	void print() const
	{
		std::cout << "n: " << n
                  << ", k: " << k
                  << ", w: " << w
                  << ", m4ri_k: " << m4ri_k
		          << std::endl;
	}
};

template<const ConfigPrange& config>
class Prange
{
    public:
    constexpr static uint32_t n = config.n;
	constexpr static uint32_t k = config.k;
	constexpr static uint32_t w = config.w;
	constexpr static uint32_t m4ri_k = config.m4ri_k;

    mzd_t* H;
    mzd_t* S;
    mzd_t* e;

    mzp_t* P;
    mzd_t* wH;
	mzd_t* wHT;
	customMatrixData* cm;

    static std::atomic_bool not_found;

    Prange(DecodingInstance& I) : H(I.A), S(I.syndrome), e(I.error)
    {
        auto ST = mzd_transpose(NULL, S);
        wH = matrix_init(n-k, n+1);
        wHT = mzd_init(wH->ncols, wH->nrows);
        P = mzp_init(n);
        cm = init_matrix_data(wH->ncols);

        mzd_t* tmp = matrix_concat(nullptr, H, ST);
        mzd_copy(wH, tmp);

        mzd_free(tmp);
        mzd_free(ST);
    }

    ~Prange()
    {
        mzp_free(P);
        mzd_free(wH);
        mzd_free(wHT);
        free_matrix_data(cm);
    }

    uint64_t run()
    {
        uint64_t loops = 0;

        while(not_found)
        {
            loops++;
            matrix_create_random_permutation(wH, wHT, P);
            matrix_echelonize_partial_plusfix(wH, m4ri_k, n-k, cm, 0, n-k, 0, P);
            if(unlikely(hamming_weight_column(wH, n) <= w)){
                bool expect = true;
                if(!not_found.compare_exchange_strong(expect, false))
                    continue;

                for(rci_t i = 0; i < n-k; i++)
                    mzd_write_bit(e, 0, P->values[i], mzd_read_bit(wH, i, n));
            }
        }
        return loops;
    }

    void bench_time(uint32_t iterations = 10000, std::string file_name = "Prange.txt")
    {
        uint32_t loops = 0;
        auto start = std::chrono::high_resolution_clock::now();
        while(loops < iterations)
        {
            loops++;
            matrix_create_random_permutation(wH, wHT, P);
            matrix_echelonize_partial_plusfix(wH, m4ri_k, n-k, cm, 0, n-k, 0, P);
            if(unlikely(hamming_weight_column(wH, n) <= w)){
                bool expect = true;
                if(!not_found.compare_exchange_strong(expect, false))
                    continue;

                for(rci_t i = 0; i < n-k; i++)
                    mzd_write_bit(e, 0, P->values[i], mzd_read_bit(wH, i, n));
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        
        if ( access( file_name.c_str(), F_OK ) == -1 )
		{
			std::ofstream tmpFile(file_name);
			tmpFile << "# n loops microsec\n";
			tmpFile.close();
		}
		std::ofstream file;
		file.open(file_name, std::ios_base::app);
		if(file.good())
			file << n << " " << iterations << " " << duration << "\n";

    }

};
template<const ConfigPrange& config>
std::atomic_bool Prange<config>::not_found{true};

#endif