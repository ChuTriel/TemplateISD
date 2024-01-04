#include <iostream>
#include <chrono>
#include <stdint.h>
#include "template_dumer.h"
#include "template_prange.h"
#include "challenges/mceliece/mce2197sp.h"

// global vars for easier adjustment
#define ITERATIONS 10000
constexpr uint32_t l = 17;
constexpr uint32_t p = 2;

void BENCH_PRANGE()
{
    std::cout << "Benching Prange.\n";
    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplatePrange config(n, k, w, addRows, newN);
    const auto B = config.parse_weight_string(eW);
    TemplatePrange<config> prange(I, B);
    prange.bench_time(ITERATIONS);
}

void BENCH_DUMER_STANDARD_PERM()
{
    std::cout << "Benching Dumer using standard permutation.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplateDumer config(n, k, w, addRows, newN, p, l, false);
    const auto B = config.parse_weight_string(eW);
    TemplateDumer<config> dumer(I, B);
    dumer.bench_time(ITERATIONS, "TemplatDumerStandardPerm.txt");
}

void BENCH_DUMER_ADVANCED_PERM()
{
    std::cout << "Benching Dumer using advanced permutation.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplateDumer config(n, k, w, addRows, newN, p, l, true);
    const auto B = config.parse_weight_string(eW);
    TemplateDumer<config> dumer(I, B);
    dumer.bench_time(ITERATIONS, "TemplatDumerAdvancedPerm.txt");
}

int main(int argc, char* argv[])
{
    std::cout << "Instance to bench: " << n << "\n";
    BENCH_PRANGE();
    BENCH_DUMER_STANDARD_PERM();
    BENCH_DUMER_ADVANCED_PERM();
    return 0;
}