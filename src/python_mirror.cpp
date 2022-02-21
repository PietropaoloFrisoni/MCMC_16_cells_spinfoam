#include "mcmc.h"
#include "hash_21j_symbols.h"
#include "omp.h"

/*
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
*/

int MH_parallel_run(std::string store_path_assigned, std::string hashed_tables_path_assigned, int dspin_assigned, int length_assigned,
                    double sigma_assigned, int burnin_assigned, int verbosity, int number_of_threads)
{

    omp_set_num_threads(number_of_threads);

    char tmp[1024];
    sprintf(tmp, "j_%.8g", ((double)(dspin_assigned) / 2.0));

    std::string path = hashed_tables_path_assigned + "/hashed_21j_symbols_" + std::string(tmp);

    auto hash = std::make_shared<Chain::Hash>();
    phmap::BinaryInputArchive ar_in(path.c_str());
    hash->phmap_load(ar_in);

#pragma omp parallel for
    for (int thread_id = 0; thread_id < number_of_threads; thread_id++)
    {

        Chain test_chain(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned,
                         sigma_assigned, burnin_assigned, verbosity, thread_id, hash);

        Metropolis_Hastings_run(test_chain);
    }

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
