#pragma once

#include "mcmc.h"
#include "hash_21j_symbols.h"
#include "omp.h"

// starts an independent Metropolis-Hastings Markov chain for each thread 
int MH_parallel_run(char *, char *, int, int, double, int, int, int);

// hashes 21j symbols 
int Hashing_21j_symbols(char *, char *, int);

