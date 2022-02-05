#pragma once

#include <string>
#include <iostream>

#include "fastwigxj.h"
#include "wigxjpf.h"
#include "mcmc.h"



void funza();

int square_2(int i);

int test_fastwig(std::string fastwing_tables_folder);

double Wigner_21j_symbol(int* key_21j, Chain &chain);
