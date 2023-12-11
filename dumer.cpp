#include <iostream>
#include <vector>
#include <omp.h>
#include "instance.h"
#include "template_dumer.h"
#include "chase_manager.h"
#include "challenges/mceliece/mce923sp.h"

int main(int argc, char* argv[])
{
    if(argc > 2){
        std::cerr << "Wrong number of input parameters! Usage: ./dumer <nr_threads>. nr_threads is optional.\n";
        return 1;
    }

    const uint32_t nr_threads = (argc == 2) ? atoi(argv[1]) : omp_get_max_threads();
    const uint32_t seed0 = time(NULL);
    srand(seed0);
    const uint32_t seed1 = rand();
    random_seed(seed1);

    DecodingInstance I(h, s, n, k);
    constexpr const uint32_t p = 2;
    constexpr const uint32_t l = 14;
    static constexpr ConfigTemplateDumer config(n, k, w, addRows, newN, p, l);
    auto B = config.parse_weight_string(eW);
    config.print();
    config.print_mem_consumption();

    auto chase = ChaseManager::getInstance().get_chase_sequence(config.comb_half, p);
    std::vector<TemplateDumer<config>*> dumers(nr_threads);
    for(size_t i = 0; i < nr_threads; i++)
        dumers[i] = new TemplateDumer<config>(I, B, chase);

    // some infos
    std::cout << "Seed0: " << seed0 << "\n";
    std::cout << "Seed1: " << seed1 << "\n";
    std::cout << "Threads: " << nr_threads << "\n"; 

    uint64_t overall_loops = 0;
    #pragma omp parallel default(none) shared(dumers, overall_loops) num_threads(nr_threads)
    {
        auto t_id = omp_get_thread_num();
        auto loops = dumers[t_id]->run();
        #pragma omp atomic
        overall_loops += loops;
    }

    if(I.check())
        std::cout << "Found solution after " << overall_loops << " loops.\n";
    else
        std::cout << "Something went wrong! Solution not correct!\n";

    mzd_print(I.error);

    return 0;
}