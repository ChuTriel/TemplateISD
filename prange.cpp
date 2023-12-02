#include <iostream>
//#include "challenges/mceliece/mce240sp.h"
#include "challenges/mceliece/mce923sp.h"
#include "template_prange.h"
#include "instance.h"
#include <omp.h>

int main()
{
    std::cout << "Starting multithreaded prange" << std::endl;
    srand(time(NULL));
    random_seed(rand());

    DecodingInstance I(h, s, n, k);
	static constexpr ConfigTemplatePrange config(n, k, w, addRows, newN);
	config.print();
    auto B = config.parse_weight_string(eW);

    // instantiate sequentially prange instances outside the parallel region because the omp critical fix
    // does not seem to fix the not-threadsafety of m4ri inits/frees...
    const auto nr_threads = omp_get_max_threads();
    std::vector<TemplatePrange<config>*> pranges(nr_threads);
    for(int i = 0; i < nr_threads; i++)
        pranges[i] = new TemplatePrange<config>(I, B);

    uint64_t overall_loops = 0;

    std::cout << "Executing program with " << nr_threads << " prange instance(s)" << std::endl;
    #pragma omp parallel default(none) shared(pranges, overall_loops) num_threads(nr_threads)
    {
        auto id = omp_get_thread_num();
        auto loops = pranges[id]->run();
        #pragma omp atomic
        overall_loops+=loops;
    }

    if(I.check()) {
        std::cout << "Found solution after " << overall_loops << " loops" << std::endl;
        mzd_print(I.error);
    }
    else
        std::cout << "Something went wrong, solution not correct!\n";

    // sequentially delete instances so mzd_free does not throw errors right and left
    for(int i = 0; i < nr_threads; i++)
        delete pranges[i];

    return 0;
}
