#include "mcmc.h"
#include "hash_21j_symbols.h"

int MH_parallel_run(char *store_path_assigned, char *hashed_tables_path_assigned, int dspin_assigned, int length_assigned,
                    double sigma_assigned, int burnin_assigned, int verbosity, int number_of_threads)
{

    omp_set_num_threads(number_of_threads);

    // loading the hash table of 21j symbols into memory
    strcat(hashed_tables_path_assigned, "/hashed_21j_symbols_");
    char tmp[1024];
    sprintf(tmp, "j_%.8g", ((double)(dspin_assigned) / 2.0));
    strcat(hashed_tables_path_assigned, tmp);

    auto hash = std::make_shared<Chain::Hash>();
    phmap::BinaryInputArchive ar_in(hashed_tables_path_assigned);
    hash->phmap_load(ar_in);

    // TODO: add warning (or better error)
    // if verbosity > 0 and there are multiple threads

#pragma omp parallel for
    for (auto thread_id = 1; thread_id <= number_of_threads; thread_id++)
    {
        Chain Markov_chain(store_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned,
                           sigma_assigned, burnin_assigned, verbosity, thread_id, hash);

        Metropolis_Hastings_run(Markov_chain);
    }

    return EXIT_SUCCESS;
}

int Hashing_21j_symbols(char *hash_tables_store_path_assigned_ext, int dspin_assigned)
{
    std::string hash_tables_store_path_assigned(hash_tables_store_path_assigned_ext);

    std::string fastwig_tables_folder = "./data_folder/fastwig_tables";

    init(fastwig_tables_folder, 0);

    Hash_21j_symbols(hash_tables_store_path_assigned, dspin_assigned);

    release();

    return EXIT_SUCCESS;
}

int Test_tables()
{
    std::string fastwig_tables_folder = "./data_folder/fastwig_tables/";

    test_fastwig(fastwig_tables_folder);

    return EXIT_SUCCESS;
}
