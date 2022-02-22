#include "mcmc.h"
#include "/home/frisus95/Scrivania/Final_project/src/python_mirror.cpp"

int main(int argc, char **argv)

{
    std::cout << "\n\nTesting threaded Metropolis-Hastings algorithm...\n\n"
              << std::endl;

    std::string store_path_assigned = "./data_folder/collected_draws";
    std::string hashed_tables_path_assigned = "./data_folder/hashed_21j";
    int dspin_assigned = 2;
    int length_assigned = 100000;
    double sigma_assigned = 0.40;
    double burnin_assigned = 10;
    int verbosity = 0;
    int number_of_threads = 1;

    MH_parallel_run(&store_path_assigned[0], &hashed_tables_path_assigned[0], dspin_assigned, length_assigned,
                    sigma_assigned, burnin_assigned, verbosity, number_of_threads);

    return EXIT_SUCCESS;
}
