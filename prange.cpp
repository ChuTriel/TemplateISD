#include <iostream>
#include <stdlib.h>
//#include "challenges/mceliece/mce240sp.h"
#include "template_prange.h"
#include "challenges/mceliece/mce923sp.h"
#include "template_prange.h"
#include "instance.h"
#include "helper.h"
#include <omp.h>
#include <time.h>

int main(int argc, char* argv[])
{
    if(argc > 2){
        std::cerr << "Wrong number of input parameters! Usage: ./prange <nr_threads>. nr_threads is optional.\n";
        return 1;
    }
     std::cout << "Starting multithreaded prange" << std::endl;


    const uint32_t nr_threads = (argc == 2) ? atoi(argv[1]) : omp_get_max_threads();
    const uint32_t seed0 = time(NULL);
    srand(seed0);
    const uint32_t seed1 = rand();
    random_seed(seed1);

    DecodingInstance I(h, s, n, k);
	static constexpr ConfigTemplatePrange config(n, k, w, addRows, newN);
	config.print();
    const auto B = config.parse_weight_string(eW);

    // instantiate sequentially prange instances outside the parallel region because the omp critical fix
    // does not seem to fix the not-threadsafety of m4ri inits/frees...
    std::vector<TemplatePrange<config>*> pranges(nr_threads);
    for(int i = 0; i < nr_threads; i++)
        pranges[i] = new TemplatePrange<config>(I, B);

    // some info
    std::cout << "seed0 = " << seed0 << std::endl;
    std::cout << "seed1 = " << seed1 << std::endl;
    std::cout << "Threads: " << nr_threads << std::endl;

    // some timings etc
    uint64_t overall_loops = 0;
    auto start_real = timer_start(CLOCK_REALTIME);
    auto start_CPU = timer_start(CLOCK_PROCESS_CPUTIME_ID);

    #pragma omp parallel default(none) shared(pranges, overall_loops) num_threads(nr_threads)
    {
        auto id = omp_get_thread_num();
        auto loops = pranges[id]->run();
        #pragma omp atomic
        overall_loops+=loops;
    }

    auto elapsed_real_ns = timer_end(start_real, CLOCK_REALTIME, true);
    auto elapsed_CPU_ns = timer_end(start_CPU, CLOCK_PROCESS_CPUTIME_ID, true);
    double elapsed_real_s = double(elapsed_real_ns) / double(1e9);
    double elapsed_CPU_s = double(elapsed_CPU_ns) / double(1e9);
    
    if(I.check())
        std::cout << "Found solution after " << overall_loops << " loops" << std::endl;
    else
        std::cout << "Something went wrong, solution not correct!\n";
    
    // in a multithreaded scenario cpu time can exceed wall time (each thread/core counts individually)
    std::cout << std::fixed << "Real time: " << elapsed_real_s << "s" << std::endl;
    std::cout << std::fixed << "CPU time: " << elapsed_CPU_s << "s" << std::endl;
    mzd_print(I.error);

    // sequentially delete instances so mzd_free does not throw errors right and left
    for(int i = 0; i < nr_threads; i++)
        delete pranges[i];

    return 0;
}
