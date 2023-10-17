#pragma once

#include <iostream>
#include <map>
#include "m4ri/m4ri.h"
#include "custom_matrix.h"

// represents a block of columns and stores the startIdx, the length, and the weight.
struct ColumnsBlock{
public:
	const int startIdx, length, weight;
	int nrColsToChoose = 0; // restructure later to make this also const, need to adapt parseWeightString in Config
	uint64_t nrCombs = 0;
	uint8_t* combs = nullptr;

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

struct ConfigSpecialPrange{
public:
	// instance parameter
	const uint32_t n, k, w;

	// info about new instance matrix
	// newN = new number of columns, additionalRows = additional parity check rows resulting from the blocks
	// Also save the weightDistribution string by adding it to the constructor and then immediately parse it?
	const uint32_t newN;
	const uint32_t additionalRows;

	// used for the method of the 4 russians algorithm (gauss)
	const uint32_t m4ri_k;

	//const char* weight_distro;

	// configuration for normal prange
	constexpr ConfigSpecialPrange(const uint32_t n,
	                              const uint32_t k,
	                              const uint32_t w) :
	        n(n),
	        k(k),
	        w(w),
	        newN(0),
	        additionalRows(0),
	        m4ri_k(matrix_opt_k(n - k, MATRIX_AVX_PADDING(n)))
	{};

	// configuration for special prange
	constexpr ConfigSpecialPrange(const uint32_t n,
	                              const uint32_t k,
	                              const uint32_t w,
	                              const uint32_t newN,
	                              const uint32_t addRows,
	                              const char* w_distro = ""):
	                              n(n),
	                              k(k),
	                              w(w),
	                              newN(newN),
	                              additionalRows(addRows),
	                              m4ri_k(matrix_opt_k(n - k + addRows, MATRIX_AVX_PADDING(newN)))
	{};

	void print() const
	{
		std::cout << "n: " << n
				  << ", k: " << k
				  << ", w: " << w
		          << ", m4ri_k: " << m4ri_k
		          << ", newN: " << newN
		          << ", addRows: " << additionalRows
				  << std::endl << std::flush;

//		int pos = 0;
//		while(weight_distro[pos] != '\0'){
//			std::cout << weight_distro[pos] << "\n";
//			pos++;
//		}
	}

	/// Parses a weightString (e.g. 0011020110) and returns a vector of the resulting columnsblocks.
	/// For now assume each digit represents the weight of a block of 32 columns. This format is also
	/// temporary as the maximum weight per block is limited to 9 and also does not allow for variable
	/// block length, but it is the easiest to implement so that the main algorithm can be tested.
	/// \param eW The weight string (Perhaps instead of a string use a constexpr const char* ?).
	/// \return A vector of ColumnBlocks that only contains the blocks in which the weight != 0.
	[[nodiscard]] std::vector<ColumnsBlock> parse_weight_string(const std::string eW = "") const
	{
		std::vector<ColumnsBlock> blocks;
		for(int i = 0; i < eW.length()-1; i++)
			if(eW.at(i) != '0')
				blocks.emplace_back(i*32, 32, eW.at(i) - '0');

		if(eW.at(eW.length()-1) != '0')
			blocks.emplace_back(((int)eW.length()-1)*32, ((int)n%32 == 0) ? 32 : (int)n%32, eW.at(eW.length()-1) - '0');

		compute_nr_nols_chosen_per_block(blocks);

		return blocks;
	}
private:
	/// Calculates and sets the number of columns to choose from each block when performing the permutation
	/// in prange. Perhaps just put this in parse function and adapt it (compute and cache in temporary structures
	/// so that nrColsToChose can be passed in the constructor so that it can be const?)git
	/// \param blocks The previously calculated blocks (without nrCols)
	void compute_nr_nols_chosen_per_block(std::vector<ColumnsBlock> &blocks) const {

		// implicitely ordered, cannot sort a vector / array using std::sort because the ColumnsBlocks use const member...
		std::multimap<float, ColumnsBlock&> map;
		// calculate the number of columns to chose from each block (first round down, then correct result if necc)
		float wD = float(int(w)), nkm = float(int(n-k+additionalRows));
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
}; // end config

template<const ConfigSpecialPrange &config>
class SpecialPrange
{
public:
    constexpr static uint32_t n = config.n;
	constexpr static uint32_t k = config.k;
	constexpr static uint32_t w = config.w;
	constexpr static uint32_t newN = config.newN;
	constexpr static uint32_t additionalRows = config.additionalRows;
	constexpr static uint32_t m4ri_k = config.m4ri_k;

	mzd_t* e;
	const mzd_t* s;
	const mzd_t* H;
	const std::vector<ColumnsBlock> blocks;


    SpecialPrange(){}

    	SpecialPrange(mzd_t* __restrict__ e,
	              const mzd_t* __restrict__ s,
	              const mzd_t* __restrict__ H,
	              const std::vector<ColumnsBlock> &blocks) noexcept :
	        e(e), s(s), H(H), blocks(blocks)
	{}

    uint64_t run()
    {
        ASSERT(blocks.size() == additionalRows);
		auto HReducedColumns = reduce_columns(H);
		auto sT = mzd_transpose(NULL, s);
		mzd_t* HReducedColumnsS = mzd_concat(NULL, HReducedColumns, sT);
		mzd_t* wH = append_rows_and_padding(HReducedColumnsS); //final matrix

		mzp_t* perm = mzp_init(newN);
		mzd_t* wHT = mzd_init(wH->ncols, wH->nrows);
		customMatrixData *matrix_data = init_matrix_data(wH->ncols);

		// used as target matrix to copy the rows from wHT into, this matrix gets transposed to wH
		// copy syndrome (once) to the newN-th+1 (idx newN) position of wHTTemp and leave the rest 0
		mzd_t* wHTTemp = mzd_init(wHT->nrows, wHT->ncols);
		mzd_transpose(wHT, wH);
		mzd_copy_row(wHTTemp, newN, wHT, newN);

		bool found = false;
		uint64_t loops = 0;
		size_t rang = 0;

		while(!found)
		{
			permutation_special_prange(perm, wH, wHT, wHTTemp);
			rang = matrix_echelonize_partial(wH, m4ri_k, n-k+additionalRows, matrix_data, 0);
			if(rang != (n-k+additionalRows))
				continue;
			loops++;
			if(unlikely(hamming_weight_column(wH, newN) <= w)) {
				found = true;
				construct_error_vector(perm, wH);
			}
		}

		// dont forget the free of structures!
		mzd_free(wHTTemp);
		mzd_free(HReducedColumns);
		mzd_free(sT);
		mzd_free(HReducedColumnsS);
		mzd_free(wH);
		mzd_free(wHT);
		mzp_free(perm);
		free_matrix_data(matrix_data);
		return loops;
    }

private:
    /// Removes the 0-weight blocks and returns a new, modified matrix
	/// \param H The original (n-k) x n parity-check matrix.
	/// \return The new matrix with dimensions (n-k) x newN.
	mzd_t* reduce_columns(const mzd_t* H){
		mzd_t* HT = mzd_transpose(NULL, H);
		mzd_t* HWT = mzd_init(newN, n-k); // so we can copy the columns as rows! :D
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

	/// Appends addRows-many additional parity check equations that can be constructed using the known weight of
	/// each columnsblock. It furthermore adds column padding (for the avx register in gauss).
	/// The input matrix is the column-reduced parity-check matrix to which the syndrome was appended.
	/// \param in The reduced parity-check matrix with appended syndrome (dimensions: (n-k) x (newN+1) ).
	/// \return The padded final working matrix (dimension: (n-k+addRows) x avx_padding(newN+1) )
	mzd_t* append_rows_and_padding(mzd_t* in){
		ASSERT(in->ncols == (newN+1)); // in = HRedCol | s

		mzd_t* ret = matrix_init(n-k + blocks.size(), in->ncols);
		mzd_copy(ret, in);
		rci_t offset = 0;
		for(rci_t rowCtr = 0; rowCtr < blocks.size(); rowCtr++){
			mzd_t* row = mzd_init(1, in->ncols); // all-zero init
			auto block = blocks.at(rowCtr);
			for(rci_t c = 0; c < block.length; c++)
				mzd_write_bit(row, 0, offset + c, 1);
			offset += block.length;
			mzd_write_bit(row, 0, in->ncols-1, block.weight % 2);
			mzd_copy_row(ret, n-k + rowCtr, row, 0);
			mzd_free(row);
		}

		return ret;
	}

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

    void construct_error_vector(mzp_t* perm, mzd_t* wH){
		mzd_t* tmp = mzd_init(1, newN);
		int offset_error = 0, offset_syndrome = 0;
		for(const auto &block : blocks){
			//const auto &block = blocks.at(j);
			//const auto perm = perms[j];
			for(int i = 0; i < block.nrColsToChoose; i++){
				mzd_write_bit(tmp, 0,perm->values[offset_error + i], mzd_read_bit(wH, offset_syndrome+i, newN));
			}
			offset_error += block.length;
			offset_syndrome += block.nrColsToChoose;
		}

		int map[newN], b = 0;
		for(const auto &block : blocks){
			for(int i = 0; i <  block.length; i++) {
				map[b] = block.startIdx + i;
				b++;
			}
		}
		// actual map to e
		for(rci_t c = 0; c < tmp->ncols; c++){
			if(mzd_read_bit(tmp, 0, c) == 1){
				mzd_write_bit(e, 0, map[c], 1);
			}
		}
		mzd_free(tmp);
	}
};