#ifndef CHASE_MANAGER_H
#define CHASE_MANAGER_H

#include <map>

class ChaseManager
{
	using T = uint16_t;
	typedef std::pair<int, int> Key;
	typedef std::map<Key, T*> Map;
public:
	ChaseManager(const ChaseManager&) = delete; //copy
	ChaseManager(ChaseManager&&) = delete; //move
	ChaseManager& operator=(const ChaseManager&) = delete;
	ChaseManager& operator=(ChaseManager&&) = delete;

	static ChaseManager& getInstance(){
		//created when used the first time, otherwise, if member var, it needs to be declared somewhere and instantiated
		// this object lives until the program's termination
		static ChaseManager instance;
		return instance;
	}

	//perhaps a pair with size and pointer to array? Or funky C-style by adding an additional output parameter?
	T* get_chase_sequence(int nn, int pp)
	{
		//std::cout << "Doing stuff for key: (" << nn << ", " << pp << ")" << std::endl;
		Key key(nn, pp);
		if(map.contains(key))
			return map[key];

		auto size = bin_coeff(nn, pp);
		auto combs = (T*)malloc(size * pp * sizeof(T));
		this->calc_p_combinations<T>(combs, nn, pp);
		map.insert(std::make_pair(key, combs));

		return combs;
	}

	uint64_t bin_coeff(uint64_t nn, uint64_t kk) noexcept
	{
		return
		        (kk > nn  ) ? 0 :       // out of range
		                (kk == 0 || kk == nn  ) ? 1 :       // edge
		                (kk == 1 || kk == nn - 1) ? nn :       // first
		                (kk + kk < nn  ) ?           // recursive:
		                (bin_coeff(nn - 1, kk - 1) * nn) / kk :       //  path to k=1   is faster
		                (bin_coeff(nn - 1, kk) * nn) / (nn - kk);      //  path to k=n-1 is faster
	}

private:
	ChaseManager() = default;
	~ChaseManager()
	{
		for(auto ele : map)
			free(ele.second);
		map.clear();
	}

	template <typename Iterator>
	bool next_combination(const Iterator first, Iterator k, const Iterator last) {
		// src: https://stackoverflow.com/questions/127704/algorithm-to-return-all-combinations-of-k-elements-from-n
		// Credits: Mark Nelson http://marknelson.us
		if ((first == last) || (first == k) || (last == k))
			return false;
		Iterator i1 = first;
		Iterator i2 = last;
		++i1;
		if (last == i1)
			return false;
		i1 = last;
		--i1;
		i1 = k;
		--i2;
		while (first != i1) {
			if (*--i1 < *i2) {
				Iterator j = k;
				while (!(*i1 < *j)) ++j;
				std::iter_swap(i1,j);
				++i1;
				++j;
				i2 = k;
				std::rotate(i1,j,last);
				while (last != j) {
					++j;
					++i2;
				}
				std::rotate(k,i2,last);
				return true;
			}
		}

		std::rotate(first,k,last);
		return false;
	}

	/// generate n over p combinations and save the result in combinations
	/// output of the form: (10, 2): [0, 1, 0, 2, 0, 3, ..., 0, 9, 1, 2, 1, 3, ....],
	/// get the i-th element easily via combinations[i*p + j]
	template<typename T>
	uint64_t calc_p_combinations(T *combinations, const uint64_t nn, const uint64_t pp) {
		uint64_t ctr = 0;

		std::vector<T> s(nn);
		for (uint64_t i = 0; i < nn; i++){
			s[i] = i;
		}

		do {
			// copy it back
			for (uint64_t i = 0; i < pp; i++) {
				combinations[ctr*pp + i] = s[i];
			}

			ctr += 1;
		} while (next_combination(s.begin(), s.begin() + pp, s.end()) && ctr < bc(nn, pp));

		return ctr;
	}

	Map map;
};

#endif // CHASE_MANAGER_H