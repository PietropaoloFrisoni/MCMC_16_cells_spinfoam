#pragma once

#include <iostream>
#include <string>
#include <vector>

class Chain {      

  public:       

    int dspin;
    int length; 
    std::vector<int> chain;
    std::vector<double> amplitudes;
    int indices[16];
    double sigma;
    int burnin;
    std::string store_path;

    Chain(int dspin_assigned, int length_assigned, double sigma_assigned, int burnin_assigned, std::string store_path_assigned) 
        : dspin(dspin_assigned), length(length_assigned), sigma(sigma_assigned), burnin(burnin_assigned), store_path(store_path_assigned) {};
        
   ~Chain(){      
        chain.clear();
        amplitudes.clear();
        // delete [] indices; SEGFAULT IF NOT COMM
        std::cout << "chain with dspin " << dspin << " destroyed"  << std::endl;
           };

};

void dmc_run(Chain &chain);