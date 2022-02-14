#include "mcmc.h"

#include <getopt.h>

int main(int argc, char **argv)

{

    std::cout << "\n\nTESTING NON-THREADED SAMPLER...\n\n"
              << std::endl;

    std::string store_path_assigned = "./data_folder/collected_draws";
    std::string hashed_tables_path_assigned = "./data_folder/hashed_21j";
    int dspin_assigned = 2;
    int length_assigned = 10;
    double sigma_assigned = 0.3;
    double burnin_assigned = 10;
    int verbosity = 2;

    Chain test_chain(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned, sigma_assigned, burnin_assigned, verbosity);

    Metropolis_Hastings_run(test_chain);

    return EXIT_SUCCESS;
}
