#ifndef DECODING_TEMPLATE_PRANGE_H
#define DECODING_TEMPLATE_PRANGE_H

#include "template_alg.h"
#include "custom_matrix.h"
#include <fstream>
#include <atomic>

// Basic configuration, no additional variables are needed
class ConfigTemplatePrange : public TemplateConfig
{
public:
	// switches between standard and advanced permutation
	const bool use_adv_perm;

	constexpr ConfigTemplatePrange(const uint32_t n,
	                              const uint32_t k,
	                              const uint32_t w,
	                              const uint32_t additional_rows,
	                              const uint32_t new_n,
								  const bool use_adv_perm = true) : TemplateConfig(n, k, w, additional_rows, new_n), use_adv_perm(use_adv_perm)
	{}

	void print() const override
	{
		TemplateConfig::print();
		std::cout << "Prange specific config: "
				  << "use_adv_perm: " << use_adv_perm
		          << std::endl;
	}
};

template<const ConfigTemplatePrange& config>
class TemplatePrange : TemplateBaseAlg
{
private:
	constexpr static uint32_t n = config.n;
	constexpr static uint32_t k = config.k;
	constexpr static uint32_t w = config.w;
	constexpr static uint32_t new_n = config.new_n;
	constexpr static uint32_t add_rows = config.additional_rows;
	constexpr static uint32_t m4ri_k = config.m4ri_k;
	constexpr static bool use_adv_perm = config.use_adv_perm;

	mzp_t* P;
	mzd_t* wHT;
	mzd_t* wHTTemp;
	customMatrixData* c_m;
public:
	static std::atomic<bool> not_found;
public:
	TemplatePrange(DecodingInstance& I, const std::vector<ColumnsBlock>& B) : TemplateBaseAlg(config, B, I)
	{
		//setup own variables, also critical part for omp
		#pragma omp critical
		{
			P = mzp_init(new_n);
			wHT = mzd_init(wH->ncols, wH->nrows);
			c_m = init_matrix_data(wH->ncols);

			//used as a tmp matrix to copy the columns that are chosen as the information set to.
			wHTTemp = mzd_init(wHT->nrows, wHT->ncols);
			mzd_transpose(wHT, wH);
			mzd_copy_row(wHTTemp, new_n, wHT, new_n);
		}
	}
	~TemplatePrange()
	{
		#pragma omp critical
		{
			mzp_free(P);
			mzd_free(wHT);
			mzd_free(wHTTemp);
			free_matrix_data(c_m);
		}
	}

	uint64_t __attribute__ ((noinline)) run() noexcept
	{
		uint64_t loops = 0;
		size_t rang = 0;

		while(not_found)
		{
			if constexpr (use_adv_perm)
			{
				permutation_special_prange(P, wH, wHT, wHTTemp);
				rang = matrix_echelonize_partial(wH, m4ri_k, n-k+add_rows, c_m, 0);
				if(rang != (n-k+add_rows))
					continue;
			} else
			{
				matrix_create_random_permutation(wH, wHT, P);
				matrix_echelonize_partial_plusfix(wH, m4ri_k, n-k+add_rows, c_m, 0, n-k+add_rows, 0, P);
			}
			loops++;
			if(unlikely(hamming_weight_column(wH, new_n) <= w)) {
				bool expect = true; //whyever i cant just write true in the function...
				if(!not_found.compare_exchange_strong(expect, false))
					break;

				if constexpr (use_adv_perm)
					construct_error_vector(P, wH);
				else
					construct_error_vector2();
			}
		}
		return loops;
	}

	void bench_time(uint32_t iterations = 10000, std::string file_name = "TemplatePrange.txt")
	{
		uint32_t loops = 0;
		uint32_t rang = 0;
		auto start = std::chrono::high_resolution_clock::now();
		while(loops < iterations)
		{
			if constexpr (use_adv_perm)
			{
				permutation_special_prange(P, wH, wHT, wHTTemp);
				rang = matrix_echelonize_partial(wH, m4ri_k, n-k+add_rows, c_m, 0);
				if(rang != (n-k+add_rows))
					continue;
			} else
			{
				matrix_create_random_permutation(wH, wHT, P);
				matrix_echelonize_partial_plusfix(wH, m4ri_k, n-k+add_rows, c_m, 0, n-k+add_rows, 0, P);
			}
			loops++;
			if(unlikely(hamming_weight_column(wH, new_n) <= w)) {
				bool expect = true;
				if(!not_found.compare_exchange_strong(expect, false))
					continue;

				if constexpr (use_adv_perm)
					construct_error_vector(P, wH);
				else
					construct_error_vector2();

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

private:
	// applies a special permutation where a precomputed number of columns of each block
	// are chosen as the information set
	void permutation_special_prange(mzp_t* perm, mzd_t* wH, mzd_t* wHT, mzd_t* tmp){
		//1: permutation
		int block_offset = 0;
		for(const auto &block : blocks){
			for(int i = 0; i < block.length; i++){
				word pos =  fastrandombytes_uint64() % (block.length-i);
				ASSERT(block_offset+i+pos < uint32_t(perm->length));
				std::swap(perm->values[block_offset+i], perm->values[block_offset+i+pos]);
			}
			block_offset+=block.length;
		}

		//2: copy first columns of each block to tmp matrix and then transpose to target matrix
		int row_dst_ctr = 0, offset = 0;
		for(const auto &block : blocks){
			for(int j = 0; j < block.nrColsToChoose; j++){
				mzd_copy_row(tmp, row_dst_ctr, wHT,perm->values[offset + j]);
				row_dst_ctr++;
			}
			offset+= block.length;
		}
		mzd_transpose(wH, tmp);
	}

	void construct_error_vector2()
	{
		mzd_t* tmp = mzd_init(1, new_n);
		for (rci_t j = 0; j < n-k+add_rows; ++j)
			mzd_write_bit(tmp, 0, P->values[j], mzd_read_bit(wH, j, new_n));

		int map[new_n], b = 0;
		for(const auto &block : blocks){
			for(int i = 0; i <  block.length; i++) {
				map[b] = block.startIdx + i;
				b++;
			}
		}
		//mzd_print(tmp);
		// actual map to e
		for(rci_t c = 0; c < tmp->ncols; c++){
			if(mzd_read_bit(tmp, 0, c) == 1){
				mzd_write_bit(e, 0, map[c], 1);
			}
		}
		mzd_free(tmp);


	}

	// reconstructs the error vector when using advanced permutation
	void construct_error_vector(mzp_t* perm, mzd_t* wH){
		mzd_t* tmp = mzd_init(1, new_n);
		int offset_error = 0, offset_syndrome = 0;
		for(const auto &block : blocks){
			for(int i = 0; i < block.nrColsToChoose; i++){
				mzd_write_bit(tmp, 0,perm->values[offset_error + i], mzd_read_bit(wH, offset_syndrome+i, new_n));
			}
			offset_error += block.length;
			offset_syndrome += block.nrColsToChoose;
		}

		int map[new_n], b = 0;
		for(const auto &block : blocks){
			for(int i = 0; i <  block.length; i++) {
				map[b] = block.startIdx + i;
				b++;
			}
		}
		//mzd_print(tmp);
		// actual map to e
		for(rci_t c = 0; c < tmp->ncols; c++){
			if(mzd_read_bit(tmp, 0, c) == 1){
				mzd_write_bit(e, 0, map[c], 1);
			}
		}
		mzd_free(tmp);
	}

};
// oh boi, are there better ways?
template<const ConfigTemplatePrange& config>
std::atomic<bool> TemplatePrange<config>::not_found{true};

#endif//DECODING_TEMPLATE_PRANGE_H
