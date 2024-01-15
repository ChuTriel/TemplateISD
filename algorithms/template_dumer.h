#ifndef DECODING_TEMPLATE_DUMER_H
#define DECODING_TEMPLATE_DUMER_H

#include "template_alg.h"
#include "chase_manager.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <bits/stdc++.h> //only gcc?
#include <atomic>
#include <omp.h>

class ConfigTemplateDumer : public TemplateConfig
{
public:
	const uint32_t p;
	const uint32_t l;

	const uint32_t bucket_factor = 5;
	const uint32_t I_size = n-k+additional_rows-l; // was nkl
	const uint32_t right_size = new_n-I_size; // mod. H without identity, was kl
	const uint32_t right_size_half = right_size/2; // was kl_half
	const uint32_t comb_half = cceil(float(right_size)/2.); // was kl_comb_half
	const uint32_t l_table_size = 1 << l;
	const uint32_t nr_combs = bc(comb_half, p);
	const uint32_t bucket_size = cceil(double(nr_combs) / double(l_table_size))*bucket_factor;
	const uint32_t nr_top_limbs = (I_size + sizeof(uint64_t)*8 -1) / (sizeof(uint64_t)*8); // mzd words are fix uint64_t
	const uint32_t top_weight = w - 2*p;

	const uint32_t nr_threads = 1;

	// switches between standard and advanced permutation
	const bool use_adv_perm;

	constexpr ConfigTemplateDumer(const uint32_t n,
	                              const uint32_t k,
	                              const uint32_t w,
	                              const uint32_t additional_rows,
	                              const uint32_t new_n,
	                              const uint32_t p,
	                              const uint32_t l,
								  const bool use_adv_perm = false,
	                              const uint32_t nr_threads = 1) : TemplateConfig(n, k, w, additional_rows, new_n),
	                                                  p(p),
	                                                  l(l),
													  use_adv_perm(use_adv_perm),
	                                                  nr_threads(nr_threads)
	{

	}
	void print() const override
	{
		TemplateConfig::print();
		std::cout << "Dumer specific config: "
		          << "l: " << l
		          << ", p: " << p
		          << ", I_size: "<< I_size
		          << ", right_size: " << right_size
		          << ", right_half: " << right_size_half
		          << ", comb_half: " << comb_half
		          << ", l_table_size: " << l_table_size
		          << ", nr_combs: " << nr_combs
		          << ", bucket_size: " << bucket_size
		          << ", nr_top_limbs: " << nr_top_limbs
		          << ", top_weight: " << top_weight
				  << ", use_adv_perm: " << use_adv_perm
		          << ", nr_threads: " << nr_threads
		          << std::endl;
	}

	// TODO: adjust all types to be dynamically uint32_t or uint64_t in config and alg
	void print_mem_consumption() const
	{
		double factor = 1.0 / 1000000; // MB for now
		std::string unit = "MB";
		double buck_mem = double(l_table_size) * double(bucket_size) * double(sizeof(uint32_t));
		double ctr_mem = double(l_table_size) * double(sizeof(uint32_t));
		double comb_mem = nr_combs * p * sizeof(uint16_t);
		double mat_syn_mem = (right_size+1) * (nr_top_limbs+1) * sizeof(uint64_t); // one more limb for l part

		buck_mem*=factor; ctr_mem*=factor; comb_mem*=factor; mat_syn_mem*=factor;
		std::cout << "Approx. memory consumption: "
		          << "bucket: " << buck_mem << unit
		          << ", counter: " << ctr_mem << unit
		          << ", comb_mem: " << comb_mem << unit
		          << ", mat+syn: " << mat_syn_mem << unit
		          << std::endl;

	}

	[[nodiscard]] std::vector<ColumnsBlock> parse_weight_string(const std::string &eW) const override
	{
		std::vector<ColumnsBlock> blocks;
		if (use_adv_perm){
			parse_basic_blocks(blocks, eW);
			compute_nr_nols_chosen_per_block(blocks);
		} else {
			blocks = TemplateConfig::parse_weight_string(eW);
		}
		return blocks;
	}

	protected:

	void compute_nr_nols_chosen_per_block(std::vector<ColumnsBlock> &blocks) const override
	{
		// compute number of cols to chose from each block
		std::multimap<float, ColumnsBlock&> map;
		float wD = float(int(top_weight)), nkm = float(int(I_size));
		int sum = 0;

		// we must also scale the weight of each block down, use double for more accuracy?
		float percentage_to_schale = wD / float(int(w));
		for(auto &block : blocks){
			float scaled_bw = float(block.weight) * percentage_to_schale;
			float resF = ((scaled_bw / wD)*nkm);
			int res = (int)resF;
			block.nrColsToChoose = (res < block.weight) ? block.weight : (res > block.length ? block.length : res);
			sum += block.nrColsToChoose;
			map.insert(std::make_pair(resF-float(res), std::ref(block)));
		}

		while(sum < (int)I_size) // shouldnt be necessary, but to be safe if multiple iterations are needed
		{
			for(auto it = map.rbegin(); it != map.rend(); it++){
				if(it->second.nrColsToChoose < it->second.length) {
					it->second.nrColsToChoose = it->second.nrColsToChoose + 1;
					sum++;
				}
				if(sum == (int)I_size)
					break;
			}
		}
		// if(sum != (int)I_size)
		// 	std::cout << "Sum does not match I_size..why??\n";

	}

};
// const important to be able to put reference "&". Otherwise config cannot be passed to class
template<const class ConfigTemplateDumer& config>
class TemplateDumer : public TemplateBaseAlg
{
	constexpr static uint32_t n = config.n;
	constexpr static uint32_t k = config.k;
	constexpr static uint32_t w = config.w;
	constexpr static uint32_t new_n = config.new_n;
	constexpr static uint32_t add_rows = config.additional_rows;
	constexpr static uint32_t m4ri_k = config.m4ri_k;
	constexpr static uint32_t l = config.l;
	constexpr static uint32_t p = config.p;
	constexpr static uint32_t I_size = config.I_size; // was nkl
	constexpr static uint64_t kl_half = config.right_size_half;
	constexpr static uint64_t kl_combinations_half = config.comb_half;
	constexpr static uint32_t l_table_size = config.l_table_size;
	constexpr static uint32_t nr_combs = config.nr_combs;
	constexpr static uint32_t bucket_size = config.bucket_size;
	constexpr static uint32_t nr_top_limbs = config.nr_top_limbs;
	constexpr static uint32_t top_weight = config.top_weight;
	constexpr static uint32_t threads = config.nr_threads;
	constexpr static bool use_adv_perm = config.use_adv_perm;

	mzp_t* P;
	mzp_t* P2;
	mzd_t* wHT;
	mzd_t* wHTTemp;
	customMatrixData* c_m;

	// old-style master thesis structures
	mzd_t* H; // extracted wH (real size)
	mzd_t* HT; //transposed so we can access column as rows
	mzd_t* botPart;
	mzd_t* topPart;

	uint32_t bucket[l_table_size*bucket_size] = {0};
	//std::atomic_uint32_t counter[l_table_size];
	uint32_t counter[l_table_size];
	uint16_t* combs;

	static std::atomic_bool not_found;

public:
	TemplateDumer(DecodingInstance& inst, const std::vector<ColumnsBlock>& blocks, uint16_t* combis =  nullptr) : TemplateBaseAlg(config, blocks, inst)
	{
		P = mzp_init(new_n);
		P2 = mzp_init(new_n - I_size);
		wHT = mzd_init(wH->ncols, wH->nrows);
		c_m = init_matrix_data(wH->ncols);

		wHTTemp = mzd_init(wHT->nrows, wHT->ncols);
		mzd_transpose(wHT, wH);
		mzd_copy_row(wHTTemp, new_n, wHT, new_n);

		// old-style master thesis structures
		H = mzd_init(n-k+add_rows, new_n - I_size +1);
		HT = mzd_init(H->ncols, H->nrows);
		botPart = mzd_init(HT->nrows, l);
		topPart = mzd_init(HT->nrows, I_size);

		combs = combis ? combis : ChaseManager::getInstance().get_chase_sequence(kl_combinations_half,p);
	}
	~TemplateDumer()
	{
		mzp_free(P);
		mzp_free(P2);
		mzd_free(wHT);
		mzd_free(wHTTemp);
		free_matrix_data(c_m);

		mzd_free(H);
		mzd_free(HT);
		mzd_free(botPart);
		mzd_free(topPart);
	}

	uint64_t __attribute__ ((noinline)) run() noexcept
	{
		uint64_t loops = 0;
		size_t rang = 0;

		while(not_found)
		{
			
			if constexpr (use_adv_perm)
			{
				adv_perm();
				rang = matrix_echelonize_partial(wH, m4ri_k, I_size, c_m, 0);
				if(rang != (I_size))
					continue ;
			} else
			{
				matrix_create_random_permutation(wH, wHT, P);
				matrix_echelonize_partial_plusfix(wH, m4ri_k, I_size, c_m, 0, I_size, 0, P);
			}
			

			loops++;

			mzd_submatrix(H, wH, 0, I_size, n-k+add_rows, new_n+1);
			mzd_transpose(HT, H);
			mzd_submatrix(botPart, HT, 0, I_size, HT->nrows, n-k+add_rows);
			mzd_submatrix(topPart, HT, 0, 0, HT->nrows, I_size);

			birthday_decoding();

			// #pragma omp parallel default(none) num_threads(threads)
			// {
			// 	birthday_decoding();
			// }
		}
		return loops;
	}

	void bench_time(uint32_t iterations = 10000, std::string file_name = "TemplateDumer.txt")
	{
		uint32_t loops = 0;
		size_t rang = 0;
		auto start = std::chrono::high_resolution_clock::now();
		while(loops < iterations)
		{
			if constexpr (use_adv_perm)
			{
				adv_perm();
				rang = matrix_echelonize_partial(wH, m4ri_k, I_size, c_m, 0);
				if(rang != (I_size))
					continue ;
			} else
			{
				matrix_create_random_permutation(wH, wHT, P);
				matrix_echelonize_partial_plusfix(wH, m4ri_k, I_size, c_m, 0, I_size, 0, P);
			}

			loops++;

			mzd_submatrix(H, wH, 0, I_size, n-k+add_rows, new_n+1);
			mzd_transpose(HT, H);
			mzd_submatrix(botPart, HT, 0, I_size, HT->nrows, n-k+add_rows);
			mzd_submatrix(topPart, HT, 0, 0, HT->nrows, I_size);
			birthday_decoding();
			// #pragma omp parallel default(none) num_threads(threads)
			// {
			// 	birthday_decoding();
			// }
		}
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
		std::cout << "Single threaded time: " << duration << std::endl;
		std::cout << "Iterations: " << loops << std::endl;

		if ( access( file_name.c_str(), F_OK ) == -1 )
		{
			std::ofstream tmpFile(file_name);
			tmpFile << "# n loops microsec l p\n";
			tmpFile.close();
		}
		std::ofstream file;
		file.open(file_name, std::ios_base::app);
		if(file.good())
			file << n << " " << loops << " " << duration << " " << l << " " << p << "\n";
	}

	uint32_t bench_perform_fixed_loops(uint32_t iterations = 10000)
	{
		uint32_t loops = 0;
		size_t rang = 0;
		while(loops < iterations)
		{
			if constexpr (use_adv_perm)
			{
				adv_perm();
				rang = matrix_echelonize_partial(wH, m4ri_k, I_size, c_m, 0);
				if(rang != (I_size))
					continue ;
			} else
			{
				matrix_create_random_permutation(wH, wHT, P);
				matrix_echelonize_partial_plusfix(wH, m4ri_k, I_size, c_m, 0, I_size, 0, P);
			}

			loops++;

			mzd_submatrix(H, wH, 0, I_size, n-k+add_rows, new_n+1);
			mzd_transpose(HT, H);
			mzd_submatrix(botPart, HT, 0, I_size, HT->nrows, n-k+add_rows);
			mzd_submatrix(topPart, HT, 0, 0, HT->nrows, I_size);
			birthday_decoding();
			// #pragma omp parallel default(none) num_threads(threads)
			// {
			// 	birthday_decoding();
			// }
		}
		return loops;
	}

private:
	void birthday_decoding()
	{
		reset_ctr();
		//#pragma omp barrier
		fill_buckets();
		//#pragma omp barrier
		join_buckets();
	}

	void reset_ctr(){
		//for(uint32_t idx = omp_get_thread_num(); idx < nr_combs; idx+=threads)
		for(uint32_t i = 0;  i < l_table_size; i += threads)
			counter[i] = 0;
	}

	void fill_buckets()
	{
		uint64_t syndrome = botPart->rows[new_n-I_size][0];
		//print_bin(syndrome);
		//for(uint32_t idx = omp_get_thread_num(); idx < nr_combs; idx+=threads)
		for(uint32_t idx = 0; idx < nr_combs; idx += threads){
			uint64_t tmp = syndrome;
			for(int i = 0; i < p; i++)
				tmp ^= botPart->rows[ combs[p*idx + i] ][0];
			ASSERT(tmp < l_table_size);

			auto pos = counter[tmp]++; //usually atomic
			//uint32_t val = pos.load(std::memory_order_relaxed);
			//pos = std::min(bucket_size-1, pos);
			pos = (pos < bucket_size-1) ? pos : bucket_size-1;
			ASSERT(pos < bucket_size);

			bucket[tmp*bucket_size + pos] = idx;
		}
	}

	void join_buckets()
	{
		//check each combination
		//for(uint32_t idx = omp_get_thread_num(); idx < nr_combs; idx+=threads)
		for(uint32_t idx = 0; idx < nr_combs; idx+=threads)
		{
			uint32_t xLow = 0;
			#pragma unroll
			for(uint32_t i = 0; i < p; i++)
				xLow ^= botPart->rows[ kl_half + combs[idx*p + i] ][0];
			ASSERT(xLow < l_table_size);

			//uint32_t numEleInBuckets = std::min(bucket_size, counter[xLow].load(std::memory_order_relaxed));
			//uint32_t numEleInBuckets = (bucket_size < counter[xLow]) ? bucket_size : counter[xLow];
			//uint32_t numEleInBuckets = std::min(bucket_size, counter[xLow].load(std::memory_order_relaxed));
			uint32_t numEleInBuckets = std::min(bucket_size, counter[xLow]);
			ASSERT(numEleInBuckets <= bucket_size);

			// go through each match in a bucket
			for(uint32_t j = 0; j < numEleInBuckets; j++){
				uint32_t idx_left = bucket[xLow * bucket_size + j];

				uint64_t xTop[nr_top_limbs] = {0};
				#pragma unroll
				for(uint32_t s = 0; s < p; s++)
					#pragma unroll
					for(uint32_t i = 0; i < nr_top_limbs; i++)
						xTop[i] ^= topPart->rows[kl_half + combs[idx*p + s]][i];

				uint64_t STop[nr_top_limbs] = {0};
				#pragma unroll
				for(uint32_t i = 0; i < nr_top_limbs; i++)
					STop[i] = topPart->rows[new_n-I_size][i];

				#pragma unroll
				for(uint32_t s = 0; s < p; s++)
					#pragma unroll
					for(uint32_t i = 0; i < nr_top_limbs; i++)
						STop[i] ^= topPart->rows[ combs[idx_left*p + s] ][i];

				uint32_t weight = 0;
				#pragma unroll
				for(uint32_t i = 0; i < nr_top_limbs; i++)
					weight+= __builtin_popcount(xTop[i]^STop[i]); //only g++/gcc

				// reconstruct solution when found
				if(unlikely(weight <= top_weight)){
					bool expect = true; //whyever i cant just write true in the function...
					if(!not_found.compare_exchange_strong(expect, false))
						break;
					//reconstruct error before permutation
					if constexpr (use_adv_perm)
						construct_solution_adv_perm(STop, xTop, idx_left, idx);
					else
						construct_solution(STop, xTop, idx_left, idx);
				}
			}
		}
	}

	void construct_solution(const uint64_t* STop, const uint64_t* xTop, const uint64_t idx_left, const uint64_t idx_right)
	{
		mzd_t* tmp_pre_perm = mzd_init(1, new_n);
		mzd_t* tmp_after_perm = mzd_init(1, new_n);
		for(uint32_t i = 0; i < nr_top_limbs; i++)
			tmp_pre_perm->rows[0][i] = xTop[i]^STop[i];

		for(uint32_t s = 0; s < p; s++)
		{
			uint32_t intPosLeft     = (I_size + combs[p*idx_left + s]) / 64;
			uint32_t shiftPosLeft   = (I_size + combs[p*idx_left + s]) % 64;
			uint32_t intPosRight    = (I_size + kl_half + combs[p*idx_right + s]) / 64;
			uint32_t shiftPosRight  = (I_size + kl_half + combs[p*idx_right + s]) % 64;

			tmp_pre_perm->rows[0][intPosLeft] ^= 1ull << shiftPosLeft;
			tmp_pre_perm->rows[0][intPosRight] ^= 1ull << shiftPosRight;
		}

		//reverse permutation
		for(uint32_t i = 0; i < P->length; i++)
			mzd_write_bit(tmp_after_perm, 0, P->values[i], mzd_read_bit(tmp_pre_perm, 0, i));

		//map solution for small instance (new_n) to the real instance (n) using the blocks
		map_small_e_to_big_e(new_n, tmp_after_perm);
		mzd_free(tmp_pre_perm);
		mzd_free(tmp_after_perm);
	}

	void construct_solution_adv_perm(const uint64_t* STop, const uint64_t* xTop, const uint64_t idx_left, const uint64_t idx_right)
	{
		// 1. Set "1" according to the indices
		mzd_t* tmp_pre_perm = mzd_init(1, new_n);
		for(uint32_t i = 0; i < nr_top_limbs; i++)
			tmp_pre_perm->rows[0][i] = xTop[i]^STop[i];

		for(uint32_t s = 0; s < p; s++)
		{
			uint32_t intPosLeft     = (I_size + combs[p*idx_left + s]) / 64;
			uint32_t shiftPosLeft   = (I_size + combs[p*idx_left + s]) % 64;
			uint32_t intPosRight    = (I_size + kl_half + combs[p*idx_right + s]) / 64;
			uint32_t shiftPosRight  = (I_size + kl_half + combs[p*idx_right + s]) % 64;

			tmp_pre_perm->rows[0][intPosLeft] ^= 1ull << shiftPosLeft;
			tmp_pre_perm->rows[0][intPosRight] ^= 1ull << shiftPosRight;
		}

		// 2. Reverse the second permutation for the right side (which is a normal permutation)
		mzd_t* tmp_rev_sec_perm = mzd_init(1, new_n);
		for(int i = 0; i < I_size; i++)
			mzd_write_bit(tmp_rev_sec_perm, 0, i, mzd_read_bit(tmp_pre_perm, 0, i));
		for(int i = 0; i < P2->length; i++)
			mzd_write_bit(tmp_rev_sec_perm, 0, I_size + P2->values[i], mzd_read_bit(tmp_pre_perm,0,I_size+i));

		// 3. Reverse first permutation ("choose x columns per block" permutation)
		mzd_t* tmp_final = mzd_init(1, new_n);
		int block_offset = 0, offset_id = 0, offset_rs = 0;
		for(const auto &block : blocks){
			for(int i = 0; i < block.nrColsToChoose; i++)
				mzd_write_bit(tmp_final, 0, P->values[block_offset + i], mzd_read_bit(tmp_rev_sec_perm,0,offset_id + i));
			for(int i = 0; i < block.length - block.nrColsToChoose; i++)
				mzd_write_bit(tmp_final, 0, P->values[block_offset + block.nrColsToChoose + i], mzd_read_bit(tmp_rev_sec_perm,0,I_size + offset_rs + i));

			block_offset += block.length;
			offset_id += block.nrColsToChoose;
			offset_rs += block.length-block.nrColsToChoose;
		}

		map_small_e_to_big_e(new_n, tmp_final);
		mzd_free(tmp_pre_perm);
		mzd_free(tmp_rev_sec_perm);
		mzd_free(tmp_final);
	}

	void adv_perm()
	{
		// 1: permutation within each block to choose first nrCols from each block for the IS
		int block_offset = 0;
		for(const auto &block : blocks){
			for(int i = 0; i < block.length; i++){
				word pos =  fastrandombytes_uint64() % (block.length-i);
				ASSERT(block_offset+i+pos < uint32_t(P->length));
				std::swap(P->values[block_offset+i], P->values[block_offset+i+pos]);
			}
			block_offset+=block.length;
		}

		// 2: copy first nrCols from each block to the IS, copy the rest to the cols after identity part
		int id_ctr = 0, after_id_ctr = 0, offset = 0;
		for(const auto &block : blocks){
			for(int j = 0; j < block.nrColsToChoose; j++){ // take first columns for IS
				mzd_copy_row(wHTTemp, id_ctr, wHT, P->values[offset+j]);
				id_ctr++;
			}
			for(int j = block.nrColsToChoose; j < block.length; j++){ // take the rest for right side
				mzd_copy_row(wHTTemp, I_size + after_id_ctr, wHT, P->values[offset+j]);
				after_id_ctr++;
			}
			offset += block.length;
		}

		// 3: permute columns after identity part
		for(int i = 0; i < P2->length; i++) P2->values[i] = i; // reset important
		for (uint32_t i = 0; i < uint32_t(P2->length-1); ++i) { // perhaps not -1?
			word pos = fastrandombytes_uint64() % (P2->length - i);
			ASSERT(i+pos < uint32_t(P2->length));
			std::swap(P2->values[i], P2->values[i+pos]);
			mzd_row_swap(wHTTemp, I_size + i, I_size+i+pos);
		}
		mzd_transpose(wH, wHTTemp);
	}


	template<typename T = uint64_t>
	    requires std::is_integral_v<T>
	void print_bin(T val){
		for(uint32_t i = 0; i < sizeof(T)*8; i++){
			std::cout << (val & 1);
			val >>= 1;
		}
	}
};
template<const ConfigTemplateDumer& config>
std::atomic_bool TemplateDumer<config>::not_found{true};

#endif