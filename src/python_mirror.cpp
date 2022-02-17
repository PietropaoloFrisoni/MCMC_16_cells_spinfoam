#include "mcmc.h"
#include "hash_21j_symbols.h"

int MH_run(char *store_path_from_python, char *hashed_tables_path_from_python, int dspin_assigned, int length_assigned, double sigma_assigned, int burnin_assigned, int verbosity)
{
    std::string store_path_assigned(store_path_from_python);
    std::string hashed_tables_path_assigned(hashed_tables_path_from_python);

    Chain test_chain(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned, sigma_assigned, burnin_assigned, verbosity);

    Metropolis_Hastings_run(test_chain);

    return EXIT_SUCCESS;
}

int Hasing(int dspin_from_python)
{
    Hash_21j_symbols(dspin_from_python);

    return EXIT_SUCCESS;
}

int Test_tables()
{
    std::string fastwig_tables_folder = "./data_folder/fastwig_tables/";

    test_fastwig(fastwig_tables_folder);

    return EXIT_SUCCESS;
}
