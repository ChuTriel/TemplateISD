#ifndef DECODING_DUMER_H
#define DECODING_DUMER_H

#include "custom_matrix.h"
#include "chase_manager.h"
#include "instance.h"
#include <cstdint>
#include <string>
#include <atomic>

// Config for the basic (non-template) dumer. Similar to ConfigTemplateDumer, refer to that
// config for more information about the parameter and variables.
class ConfigDumer
{
    public:

    // instance and alg. parameter
    const uint32_t n, k, w;
    const uint64_t p, l;
    const uint32_t m4ri_k;

    // helper variable and variables for hashmap (bucket sort)
    const uint32_t bucket_factor = 3;
    const uint32_t I_size = n-k-l;
    const uint32_t right_size = n - I_size;
    const uint32_t right_size_half = right_size / 2;
    const uint32_t comb_half = cceil(float(right_size)/2.);
    const uint64_t l_table_size = 1ull << l;
    const uint64_t nr_combs = bc(comb_half, p);
    const uint32_t bucket_size = cceil(double(nr_combs) / double(l_table_size))*bucket_factor;
    const uint32_t nr_top_limbs = (I_size + sizeof(uint64_t)*8 -1) / (sizeof(uint64_t)*8);
    const uint32_t top_weight = w - 2*p;
    const uint32_t nr_threads;

    // Constructor
    constexpr ConfigDumer(const uint32_t n,
                        const uint32_t k,
                        const uint32_t w,
                        const uint64_t l,
                        const uint64_t p,
                        const uint32_t nr_threads = 1) :
                        n(n),
                        k(k),
                        w(w),
                        l(l),
                        p(p),
                        nr_threads(nr_threads),
                        m4ri_k(matrix_opt_k(n-k, MATRIX_AVX_PADDING(n)))
    {}

    // Prints dumer config.
    void print() const
    {
        std::cout << "Dumer config: "
            << "n: " << n
            << ", k: " << k
            << ", w: " << w
            << ", l: " << l
            << ", p: " << p
            << ", m4ri_k: " << m4ri_k
            << ", I_size: "<< I_size
            << ", right_size: " << right_size
            << ", right_half: " << right_size_half
            << ", comb_half: " << comb_half
            << ", l_table_size: " << l_table_size
            << ", nr_combs: " << nr_combs
            << ", bucket_size: " << bucket_size
            << ", nr_top_limbs: " << nr_top_limbs
            << ", top_weight: " << top_weight
            << std::endl;
    }

};

// Basic (non-template) dumer class. Similar to TemplateDumer, for more information refer to
// that class.
template <const class ConfigDumer &config>
class Dumer
{
    public:

    constexpr static uint32_t n = config.n;
	constexpr static uint32_t k = config.k;
	constexpr static uint32_t w = config.w;
	constexpr static uint64_t l = config.l;
	constexpr static uint64_t p = config.p;
    constexpr static uint32_t m4ri_k = config.m4ri_k;
	constexpr static uint32_t I_size = config.I_size; // was nkl
	constexpr static uint64_t kl_half = config.right_size_half;
	constexpr static uint64_t kl_combinations_half = config.comb_half;
	constexpr static uint64_t l_table_size = config.l_table_size;
	constexpr static uint64_t nr_combs = config.nr_combs;
	constexpr static uint32_t bucket_size = config.bucket_size;
	constexpr static uint32_t nr_top_limbs = config.nr_top_limbs;
	constexpr static uint32_t top_weight = config.top_weight;
    constexpr static uint32_t nr_threads = config.nr_threads;

    // instance structures
    mzd_t* A, *s, *e;
    
    // permutation and working matrix
    mzp_t* P;
    mzd_t* wH;
    mzd_t* wHT;
    customMatrixData* cm;

    // helper structures for birthday decoding
    mzd_t* H; 
	mzd_t* HT; 
	mzd_t* botPart;
	mzd_t* topPart;

    // hashmap (bucket sort) structures
    uint32_t bucket[l_table_size*bucket_size] = {0};
    uint32_t counter[l_table_size];

    // case array and static found variable (multi-threaded purposes)
    uint16_t* combs;
    static std::atomic_bool not_found;

    // Constructor.
    Dumer(DecodingInstance& I, uint16_t* combis = nullptr) : A(I.A), s(I.syndrome), e(I.error)
    {
        wH = matrix_init(n-k, n+1);
        wHT = mzd_init(wH->ncols, wH->nrows);
        P = mzp_init(n);
        cm = init_matrix_data(wH->ncols);

        mzd_t* tmp = matrix_concat(nullptr, I.A, I.syndromeT);
        mzd_copy(wH, tmp);
        mzd_free(tmp);

        H = mzd_init(n-k, k+l+1);
        HT = mzd_init(H->ncols, H->nrows);
        botPart = mzd_init(HT->nrows, l);
        topPart = mzd_init(HT->nrows, I_size); //nkl

        combs = combis ? combis : ChaseManager::getInstance().get_chase_sequence(kl_combinations_half, p);
    }

    // Destuctor
    ~Dumer()
    {
        mzp_free(P);
        mzd_free(wHT);
        mzd_free(wH);
        free_matrix_data(cm);

        mzd_free(H);
        mzd_free(HT);
        mzd_free(botPart);
        mzd_free(topPart);
    }

    // Main method that executes dumer. Returns the number of executed loops.
    uint64_t __attribute__ ((noinline)) run() noexcept
    {
        uint64_t loops = 0;

        while(not_found)
        {
            loops++;
            // permutation + gauss
            matrix_create_random_permutation(wH, wHT, P);
            matrix_echelonize_partial_plusfix(wH, m4ri_k, I_size, cm, 0, I_size, 0, P);

            // extraction of neccessary structures for the birthday decoding
            mzd_submatrix(H, wH, 0, I_size, n-k, n+1);
            mzd_transpose(HT, H);
            mzd_submatrix(botPart, HT, 0, I_size, HT->nrows, n-k);
            mzd_submatrix(topPart, HT, 0, 0, HT->nrows, I_size);
            // pragma omp
            birthday_decoding();
        }

        return loops;
    }

    // Helper method to bench the runtume. Executes "iterations" many loops and writes the
    // runtime in the specified file.
    void bench_time(uint32_t iterations = 10000, std::string file_name = "Dumer.txt")
    {
        uint32_t loops = 0;
        auto start = std::chrono::high_resolution_clock::now();
        while(loops < iterations)
        {
            loops++;
            matrix_create_random_permutation(wH, wHT, P);
            matrix_echelonize_partial_plusfix(wH, m4ri_k, I_size, cm, 0, I_size, 0, P);

            mzd_submatrix(H, wH, 0, I_size, n-k, n+1);
            mzd_transpose(HT, H);
            mzd_submatrix(botPart, HT, 0, I_size, HT->nrows, n-k);
            mzd_submatrix(topPart, HT, 0, 0, HT->nrows, I_size);
            // pragma omp
            birthday_decoding();
        }
        auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

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

    private:

    // Main birthday decoding method.
    void birthday_decoding()
    {
        reset_ctr();
        fill_buckets();
        join_buckets();
    }

    // Step 1: resets the counter.
    void reset_ctr()
    {
        for(int i = 0; i < l_table_size; i+=nr_threads)
            counter[i] = 0;
    }

    // Step 2: builds L2=H1e1+s on l bit and sorts it implicitely by using bucket sort.
    void fill_buckets()
    {
        uint64_t syndrome = botPart->rows[n-I_size][0];
		//print_bin(syndrome);
		//for(uint32_t idx = omp_get_thread_num(); idx < nr_combs; idx+=threads)
		for(uint32_t idx = 0; idx < nr_combs; idx += 1){
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

    // Step 3: builds L2=H2e2 and matches on l bit. If a collision is found, the remaining n-k-l bits
    // are computed and checked for the target weight of w-2p.
    void join_buckets()
    {
        for(uint32_t idx = 0; idx < nr_combs; idx+=nr_threads)
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
					STop[i] = topPart->rows[n-I_size][i]; //or k+l

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

                    construct_solution(STop, xTop, idx_left, idx);
				}
			}
		}

    }

    // Reconstructs the solution. It is called in join_buckets.
    void construct_solution(const uint64_t* STop, const uint64_t* xTop, const uint64_t idx_left, const uint64_t idx_right)
    {
		mzd_t* tmp_pre_perm = mzd_init(1, n);
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
			mzd_write_bit(e, 0, P->values[i], mzd_read_bit(tmp_pre_perm, 0, i));

        mzd_free(tmp_pre_perm);
    }
};
template<const class ConfigDumer& config>
std::atomic_bool Dumer<config>::not_found{true};

#endif