#include <iostream>
#include "challenges/mceliece/mce240sp.h"
#include "special_prange.h"

int main()
{
    srand(time(NULL));
    random_seed(rand());

    mzd_t* HT = mzd_from_str(n, n-k, h);
	mzd_t* H = mzd_init(n-k, n);
	mzd_transpose(H, HT);
	mzd_t* sRV = mzd_from_str(1, n-k, s);
	mzd_t* e = mzd_init(1, n);

	static constexpr ConfigSpecialPrange config(n, k, w, newN, addRows);
	auto blocks = config.parse_weight_string(eW);
	config.print();
    std::cout << "Parsed blocks: \n";
    for(const auto &block : blocks)
        block.print();
    
    auto sp = SpecialPrange<config>(e, sRV, H, blocks);
    auto loops = sp.run();
    std::cout << "Found e after " << loops << " loops.\n";
    mzd_print(e);

    return 0;
}