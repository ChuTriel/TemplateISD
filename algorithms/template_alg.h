#ifndef DECODING_TEMPLATE_ALG_H
#define DECODING_TEMPLATE_ALG_H

#include "custom_matrix.h"
#include "instance.h"
#include <vector>
#include <map>

struct ColumnsBlock{
public:
	const int startIdx, length, weight;
	int nrColsToChoose = 0; // restructure later to make this also const, need to adapt parseWeightString in Config

	ColumnsBlock(const int startIdx,
	             const int length,
	             const int weight) :
	             startIdx(startIdx),
	             length(length),
	             weight(weight)
	{}

	// convenience method for debugging
	void print() const
	{
		std::cout << "startIdx: " << startIdx
		          << ", length: " << length
		          << ", weight: " << weight
		          << ", nrColsChosen: " << nrColsToChoose << std::endl << std::flush;
	}
};

// should also hold the weightstring in addition to n,k,w,m4ri_k?
class TemplateConfig
{
public:
	// enables parity-check equations as additional rows (1 row for each non-zero block)
	const bool enable_pce = true;

	const uint32_t n, k, w;
	const uint32_t additional_rows, new_n;
	const uint32_t m4ri_k;

	constexpr TemplateConfig(const uint32_t n,
	                         const uint32_t k,
	                         const uint32_t w,
	                         const uint32_t additional_rows,
	                         const uint32_t new_n) :
 							n(n),
	                        k(k),
	                        w(w),
	                        additional_rows(enable_pce ? additional_rows : 0),
	                        new_n(new_n),
	                        m4ri_k(matrix_opt_k(n - k, MATRIX_AVX_PADDING(n))) // new_n and n-k+addrows?
	{

	}


	virtual void print() const
	{
		std::cout << "Default template config: "
				<< "n: " << n
				<< ", k: " << k
				<< ", w: " << w
		        << ", n-k: " << n-k
		        << ", m4ri_k: " << m4ri_k
				<< ", enable_pce: " << enable_pce
				<< ", add_rows: " << additional_rows
				<< ", new_n: " << new_n
				<< std::endl;
	}

	// default arg should also be explicitely specified when overriding this function
	// or: make compute_nr_cols... function protected and virtual? -> the actual parsing should be the same
	// for each template alg. and only the number of cols chosen in each block might be different
	[[nodiscard]] virtual std::vector<ColumnsBlock> parse_weight_string(const std::string &eW) const
	{
		std::vector<ColumnsBlock> blocks;
		TemplateConfig::parse_basic_blocks(blocks, eW);
		TemplateConfig::compute_nr_nols_chosen_per_block(blocks);
		return blocks;
	}

protected:

	// should be the same in every parse_weight_string function, derived configs can use this function as
	// the first parsing step.
	void parse_basic_blocks(std::vector<ColumnsBlock> &blocks, const std::string &eW) const
	{
		for(int i = 0; i < eW.length()-1; i++)
			if(eW.at(i) != '0')
				blocks.emplace_back(i*32, 32, eW.at(i) - '0');

		if(eW.at(eW.length()-1) != '0')
			blocks.emplace_back(((int)eW.length()-1)*32, ((int)n%32 == 0) ? 32 : (int)n%32, eW.at(eW.length()-1) - '0');
	}

	virtual void compute_nr_nols_chosen_per_block(std::vector<ColumnsBlock> &blocks) const
	{
		// implicitely ordered, cannot sort a vector / array using std::sort because the ColumnsBlocks use const member...
		std::multimap<float, ColumnsBlock&> map;
		// calculate the number of columns to chose from each block (first round down, then correct result if necc)
		float wD = float(int(w)), nkm = float(int(n-k+additional_rows));
		int sum = 0;
		for(auto &block : blocks){
			float resF = ((float(block.weight) / wD)*nkm);
			int res = (int)resF;
			block.nrColsToChoose = (res < block.weight) ? block.weight : (res > block.length ? block.length : res);
			sum += block.nrColsToChoose;
			map.insert(std::make_pair(resF-float(res), std::ref(block)));
		}

		while(sum < (int)nkm) // shouldnt be necessary, but to be safe if multiple iterations are needed
		{
			for(auto it = map.rbegin(); it != map.rend(); it++){
				if(it->second.nrColsToChoose < it->second.length) {
					it->second.nrColsToChoose = it->second.nrColsToChoose + 1;
					sum++;
				}
				if(sum == (int)nkm)
					break;
			}
		}
		ASSERT(sum == (int)nkm);
	}
};

//template<const class TemplateConfig& t>
class TemplateBaseAlg
{
	//Add steps to modify matrix and have variables (wH,...) as member vars
public:
	mzd_t* H; //parity-check matrix
	mzd_t* S; //syndrome (row vector)
	mzd_t* e; //error
	mzd_t* wH; // working matrix (column-reduced, syndrome attached, parity-check rows added, padded for avx)
	std::vector<ColumnsBlock> blocks;


	TemplateBaseAlg(const TemplateConfig& config,
	                const std::vector<ColumnsBlock>& blocks,
	                DecodingInstance& I) : H(I.A), S(I.syndrome), e(I.error), blocks(blocks)
	{
		//potentionally critical code section for omp
		#pragma omp critical
		{
			wH = prepare_template_matrix(config);
		}
	}

	~TemplateBaseAlg()
	{
		#pragma omp critical
		{
			mzd_free(wH);
		}
	}

private:
	// absolute final (reduced columns, appended rows, avx padding for gauss)
	mzd_t* prepare_template_matrix(const TemplateConfig& c)
	{
		auto H_col_reduced = reduce_columns(c);
		auto ST = mzd_transpose(NULL, S);
		mzd_t* H_col_reduced_s = mzd_concat(NULL, H_col_reduced, ST);
		auto ret = append_rows_and_padding(c, H_col_reduced_s);

		mzd_free(H_col_reduced);
		mzd_free(ST);
		mzd_free(H_col_reduced_s);
		return ret;
	}

	mzd_t* reduce_columns(const TemplateConfig& c)
	{
		mzd_t* HT = mzd_transpose(NULL, H);
		mzd_t* HWT = mzd_init(c.new_n, c.n-c.k); // so we can copy the columns as rows! :D
		mzd_t* HW = mzd_init(HWT->ncols, HWT->nrows);
		int currentRowHWT = 0;
		for(const auto &block : blocks){
			for(rci_t r = 0; r < block.length; r++){ // use transposed matrices so copy_row can be used
				mzd_copy_row(HWT, currentRowHWT, HT, block.startIdx + r);
				currentRowHWT++;
			}
		}
		mzd_transpose(HW, HWT);
		mzd_free(HT);
		mzd_free(HWT);
		return HW;
	}

	mzd_t* append_rows_and_padding(const TemplateConfig& c, mzd_t* in)
	{
		ASSERT(in->ncols == (c.new_n+1)); // in = HRedCol | s

		mzd_t* ret = matrix_init(c.n-c.k + c.additional_rows, in->ncols);
		mzd_copy(ret, in);
		rci_t offset = 0;
		for(rci_t rowCtr = 0; rowCtr < c.additional_rows; rowCtr++){
			mzd_t* row = mzd_init(1, in->ncols); // all-zero init
			auto block = blocks.at(rowCtr);
			for(rci_t c = 0; c < block.length; c++)
				mzd_write_bit(row, 0, offset + c, 1);
			offset += block.length;
			mzd_write_bit(row, 0, in->ncols-1, block.weight % 2);
			mzd_copy_row(ret, c.n-c.k + rowCtr, row, 0);
			mzd_free(row);
		}
		return ret;
	}
protected:
	void map_small_e_to_big_e(uint32_t newN, mzd_t* small_e)
	{
		int map[newN], b = 0;
		for(const auto &block : blocks){
			for(int i = 0; i <  block.length; i++) {
				map[b] = block.startIdx + i;
				b++;
			}
		}
		// actual map to e
		for(rci_t c = 0; c < small_e->ncols; c++){
			if(mzd_read_bit(small_e, 0, c) == 1){
				mzd_write_bit(e, 0, map[c], 1);
			}
		}
	}

};


#endif//DECODING_TEMPLATE_ALG_H
