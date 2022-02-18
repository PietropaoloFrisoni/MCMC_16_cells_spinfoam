#include "mcmc.h"
#include "hash_21j_symbols.h"

int MH_run(char *draws_store_path_from_python, char *hashed_tables_path_from_python, int dspin_assigned, int length_assigned,
           double sigma_assigned, int burnin_assigned, int verbosity, int thread_id)
{
    std::string store_path_assigned(draws_store_path_from_python);
    std::string hashed_tables_path_assigned(hashed_tables_path_from_python);

    Chain test_chain(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned,
                     sigma_assigned, burnin_assigned, verbosity, thread_id);

    Metropolis_Hastings_run(test_chain);

    return EXIT_SUCCESS;
}

int Hashing_21j_symbols(char *hash_tables_store_path_from_python, int dspin_from_python)
{
    std::string hash_tables_store_path_assigned(hash_tables_store_path_from_python);

    std::string fastwig_tables_folder = "./data_folder/fastwig_tables";

    init(fastwig_tables_folder, 0);

    Hash_21j_symbols(hash_tables_store_path_assigned, dspin_from_python);

    release();

    return EXIT_SUCCESS;
}

int Test_tables()
{
    std::string fastwig_tables_folder = "./data_folder/fastwig_tables/";

    test_fastwig(fastwig_tables_folder);

    return EXIT_SUCCESS;
}
