#include "mcmc.h"

int main(int argc, char *argv[])
{

    std::cout << "\n\nTESTING NON-THREADED SAMPLER...\n\n"
              << std::endl;

    // in general I want these two strings to be different!
    std::string fastwig_tables_folder = argv[1];
    std::string store_path_assigned = argv[1];
    int dspin_assigned = 2;
    int length_assigned = 100;
    double sigma_assigned = 0.3;
    double burnin_assigned = 10;
    int verbosity = 2;

    Chain test_chain(dspin_assigned, length_assigned, sigma_assigned, burnin_assigned, store_path_assigned, verbosity);

    dmc_run(test_chain);
}
