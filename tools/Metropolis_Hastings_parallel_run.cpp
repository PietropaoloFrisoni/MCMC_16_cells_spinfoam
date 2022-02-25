#include "python_mirror.h"

int main(int argc, char **argv)

{
    if (argc < 9)
    {
        fprintf(stderr, "\n\nError: too few arguments.\nRun the script as: %s [NUMBER_OF_THREADS] [DSPIN] [LENGTH] [SIGMA] [BURNIN] [VERBOSITY] [DRAWS_STORE_PATH] [HASH_TABLES_STORE_PATH]\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // flag parsing
    // TODO: add parameter checks
    int const number_of_threads = std::stoi(argv[1]);
    int const dspin_assigned = std::stoi(argv[2]);
    int const length_assigned = std::stoi(argv[3]);
    double const sigma_assigned = std::stof(argv[4]);
    int const burnin_assigned = std::stoi(argv[5]);
    int const verbosity_assigned = std::stoi(argv[6]);
    std::string draws_store_path_assigned = argv[7];
    std::string hash_tables_path_assigned = argv[8];

    MH_parallel_run(&draws_store_path_assigned[0], &hash_tables_path_assigned[0], dspin_assigned, length_assigned,
                    sigma_assigned, burnin_assigned, verbosity_assigned, number_of_threads);

    return EXIT_SUCCESS;
}
