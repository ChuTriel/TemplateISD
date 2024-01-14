#include <iostream>
#include <vector>
#include <omp.h>
#include "instance.h"
#include "template_prange.h"
#include "template_dumer.h"
#include "prange.h"
#include "chase_manager.h"
#include "challenges/mceliece/mce640sp.h"
//#include "challenges/mceliece/mce808sp.h"
//#include "challenges/mceliece/mce982sp.h"
//#include "challenges/mceliece/mce1101sp.h"
//#include "challenges/mceliece/mce1223sp.h"
//#include "challenges/mceliece/mce1473sp.h"
//#include "challenges/mceliece/mce1665sp.h"
//#include "challenges/mceliece/mce1995sp.h"
//#include "challenges/mceliece/mce2197sp.h"
//#include "challenges/mceliece/mce3488sp.h"

// #include "challenges/mceliece/mce2403sp.h"
// #include "challenges/mceliece/mce2681sp.h"
// #include "challenges/mceliece/mce2893sp.h"
// #include "challenges/mceliece/mce3108sp.h"
// #include "challenges/mceliece/mce3180sp.h"

// global vars for easier adjustment
#define ITERATIONS 10000
constexpr uint32_t l = 14;
constexpr uint32_t p = 2;


void BENCH_STANDARD_PRANGE()
{
    std::cout << "Benching Standard Prange.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigPrange config(n, k, w);
    config.print();
    auto prange = new Prange<config>(I);
    prange->bench_time(ITERATIONS, "StandardPrange.txt");
    delete prange;
}

void BENCH_PRANGE_STANDARD_PERM()
{
    std::cout << "Benching Prange using standard permutation.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplatePrange config(n, k, w, addRows, newN, false);
    const auto B = config.parse_weight_string(eW);
    config.print();
    auto prange = new TemplatePrange<config>(I, B);
    prange->bench_time(ITERATIONS, "TemplatePrangeStandardPerm.txt");
    delete prange;
}

void BENCH_PRANGE_ADVANCED_PERM()
{
    std::cout << "Benching Prange using advanced permutation.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplatePrange config(n, k, w, addRows, newN, true);
    const auto B = config.parse_weight_string(eW);
    config.print();
    auto prange = new TemplatePrange<config>(I, B);
    prange->bench_time(ITERATIONS, "TemplatePrangeAdvancedPerm.txt");
    delete prange;
}


void BENCH_DUMER_STANDARD_PERM()
{
    std::cout << "Benching Dumer using standard permutation.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplateDumer config(n, k, w, addRows, newN, p, l, false);
    const auto B = config.parse_weight_string(eW);
    config.print();
    auto dumer = new TemplateDumer<config>(I, B);
    dumer->bench_time(ITERATIONS, "TemplatDumerStandardPerm.txt");
    delete dumer;
}

void BENCH_DUMER_ADVANCED_PERM()
{
    std::cout << "Benching Dumer using advanced permutation.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplateDumer config(n, k, w, addRows, newN, p, l, true);
    const auto B = config.parse_weight_string(eW);
    config.print();
    auto dumer = new TemplateDumer<config>(I, B);
    dumer->bench_time(ITERATIONS, "TemplatDumerAdvancedPerm.txt");
    delete dumer;
}

int main(int argc, char* argv[])
{
    srand(time(NULL));
    random_seed(rand());

    std::cout << "Instance to bench: " << n << "\n";
    BENCH_STANDARD_PRANGE();
    BENCH_PRANGE_STANDARD_PERM();
    BENCH_PRANGE_ADVANCED_PERM();
    BENCH_DUMER_STANDARD_PERM();
    BENCH_DUMER_ADVANCED_PERM();
    
    return 0;
}