#include "mcmc.h"
#include "/home/frisus95/Scrivania/Final_project/src/python_mirror.cpp"

int main(int argc, char **argv)

{
    std::cout << "\n\nTesting threaded Metropolis-Hastings algorithm...\n\n"
              << std::endl;

    std::string store_path_assigned = "./data_folder/collected_draws";
    std::string hashed_tables_path_assigned = "./data_folder/hashed_21j";
    int dspin_assigned = 5;
    int length_assigned = 100;
    double sigma_assigned = 0.40;
    double burnin_assigned = 10;
    int verbosity = 0;
    int number_of_threads = 6;

    MH_parallel_run(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned,
                    sigma_assigned, burnin_assigned, verbosity, number_of_threads);

    /*
       Chain test_chain(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned,
                        sigma_assigned, burnin_assigned, verbosity, thread);

       Metropolis_Hastings_run(test_chain);


   char tmp[1024];

    sprintf(tmp, "j_%.8g", ((double)(dspin_assigned) / 2.0));

       auto path = hashed_tables_path_assigned + "/hashed_21j_symbols_" + std::string(tmp);

      

    MH_run("./data_folder/collected_draws", "path", dspin_assigned, length_assigned,
           sigma_assigned, burnin_assigned, verbosity, 6);


            */

    return EXIT_SUCCESS;
}
