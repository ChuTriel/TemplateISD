#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>
#include "instance.h"
#include "template_prange.h"
#include "template_dumer.h"
#include "prange.h"
#include "chase_manager.h"
//#include "challenges/mceliece/mce640sp.h"
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
#include "challenges/mceliece/mce2062sp.h"

// global vars for easier adjustment
#define ITERATIONS 10000
#define NR_THREADS 16
constexpr uint32_t l = 16;
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

void BENCH_PRANGE_ADVANCED_PERM_MULTI()
{
    std::cout << "Benching Prange (adv perm) multithreaded.\n";

    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplatePrange config(n, k, w, addRows, newN, true);
    const auto B = config.parse_weight_string(eW);
    config.print();

    std::vector<TemplatePrange<config>*> pranges(NR_THREADS);
    for(int i = 0; i < NR_THREADS; i++)
        pranges[i] = new TemplatePrange<config>(I, B);

    uint32_t overall_loops = 0;
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel default(none) shared(pranges, overall_loops) num_threads(NR_THREADS)
    {
        auto id = omp_get_thread_num();
        auto loops = pranges[id]->bench_perform_fixed_loops(ITERATIONS);
        #pragma omp atomic
        overall_loops+=loops;
    }
    auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    std::cout << "Nr_threads: " << NR_THREADS << std::endl;
    std::cout << "Iteration/thrd: " << ITERATIONS << std::endl;
    std::cout << "All iterations: " << overall_loops << std::endl;
    std::cout << "Time: " << duration << std::endl;

    for(int i = 0; i < NR_THREADS; i++)
        delete pranges[i];

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

void BENCH_DUMER_ADVANCED_PERM_MULTI()
{
    std::cout << "Benching Dumer (adv perm) multithreaded.";
    DecodingInstance I(h, s, n, k);
    static constexpr ConfigTemplateDumer config(n, k, w, addRows, newN, p, l, true);
    const auto B = config.parse_weight_string(eW);
    config.print();
    auto chase = ChaseManager::getInstance().get_chase_sequence(config.comb_half, p);
    std::vector<TemplateDumer<config>*> dumers(NR_THREADS);
    for(size_t i = 0; i < NR_THREADS; i++)
        dumers[i] = new TemplateDumer<config>(I, B, chase);
    
    uint32_t overall_loops = 0;
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel default(none) shared(dumers, overall_loops) num_threads(NR_THREADS)
    {
        auto t_id = omp_get_thread_num();
        auto loops = dumers[t_id]->bench_perform_fixed_loops(ITERATIONS);
        #pragma omp atomic
        overall_loops+=loops;
    }
    auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    std::cout << "Nr_threads: " << NR_THREADS << std::endl;
    std::cout << "Iteration/thrd: " << ITERATIONS << std::endl;
    std::cout << "All iterations: " << overall_loops << std::endl;
    std::cout << "Time: " << duration << std::endl;

    for(size_t i = 0; i < NR_THREADS; i++)
        delete dumers[i];
}

int main(int argc, char* argv[])
{
    srand(time(NULL));
    random_seed(rand());

    std::cout << "Instance to bench: " << n << "\n";
    BENCH_STANDARD_PRANGE();
    BENCH_PRANGE_STANDARD_PERM();
    BENCH_PRANGE_ADVANCED_PERM();
    BENCH_PRANGE_ADVANCED_PERM_MULTI();
    BENCH_DUMER_STANDARD_PERM();
    BENCH_DUMER_ADVANCED_PERM();
    BENCH_DUMER_ADVANCED_PERM_MULTI();
    
    return 0;
}