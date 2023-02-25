#pragma once

#include "mcmc.h"
#include "hash_21j_symbols.h"
#include "omp.h"

// starts an independent Metropolis-Hastings Markov chain for each thread 
int MH_parallel_run(char *data_folder_path_assigned_char, char *hashed_tables_path_assigned_charp, int dspin_assigned, int length_assigned,
                    double sigma_assigned, int burnin_assigned, int verbosity, int number_of_threads);

// hashes 21j symbols 
int Hashing_21j_symbols(char *hash_tables_store_path_assigned_ext, char *fastwig_tables_path_assigned_ext, int dspin_assigned);

