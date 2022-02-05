#include <string>
#include <iostream>

#include "mcmc.h"



int main(int argc, char *argv[])
{

std::cout << "\n\n\nTESTING NON-THREADED SAMPLER...\n\n\n" << std::endl;

std::string fastwig_tables_folder = argv[1];
std::string store_path_assigned = argv[1];
int dspin_assigned = 2;
int length_assigned = 100;  
double sigma_assigned  = 0.4; 
double burnin_assigned = 10;

Chain test_chain (dspin_assigned, length_assigned, sigma_assigned, burnin_assigned, store_path_assigned);
  
}