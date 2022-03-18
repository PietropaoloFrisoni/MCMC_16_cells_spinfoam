#include "python_mirror.h"

int MH_parallel_run(char *data_folder_path_assigned_char, char *hashed_tables_path_assigned_char, int dspin_assigned, int length_assigned,
                    double sigma_assigned, int burnin_assigned, int verbosity, int number_of_threads)
{

    // TODO: clean this mess between strings and char pointers

    std::string data_folder_path_assigned(data_folder_path_assigned_char);
    std::string hashed_tables_path_assigned(hashed_tables_path_assigned_char);

    char tmp[1024];
    sprintf(tmp, "j_%.8g", ((double)(dspin_assigned) / 2.0));
    std::string path = hashed_tables_path_assigned + "/hashed_21j_symbols_" + std::string(tmp);

    if (!std::filesystem::exists(path))
    {
        error("No hash table of 21j symbols with same provided spin found in folder");
    }
    else
    {
        auto hash = std::make_shared<Chain::Hash>();
        phmap::BinaryInputArchive ar_in(path.c_str());
        hash->phmap_load(ar_in);

#pragma omp parallel for
        for (auto thread_id = 0; thread_id < number_of_threads; thread_id++)
        {

            Chain test_chain(data_folder_path_assigned, hashed_tables_path_assigned, dspin_assigned, length_assigned,
                             sigma_assigned, burnin_assigned, verbosity, thread_id + 1, hash);

            Metropolis_Hastings_run(test_chain);
        }
    }

    return EXIT_SUCCESS;
}

int Hashing_21j_symbols(char *hash_tables_store_path_assigned_ext, char *fastwig_tables_path_assigned_ext, int dspin_assigned)
{
    std::string hash_tables_store_path_assigned(hash_tables_store_path_assigned_ext);

    std::string fastwig_tables_folder(fastwig_tables_path_assigned_ext);

    init(fastwig_tables_folder);

    Hash_21j_symbols(hash_tables_store_path_assigned, dspin_assigned);

    release();

    return EXIT_SUCCESS;
}
