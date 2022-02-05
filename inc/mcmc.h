#pragma once

#include <iostream>
#include <string>

class Chain {      

  public:       

    int dspin;
    int length;
    int * chain;
    int * indices;
    int burnin;
    double * amplitudes;
    double sigma;
    std::string store_path;

    Chain(int dspin_assigned, int length_assigned, double sigma_assigned, int burnin_assigned, std::string store_path_assigned) {
          dspin = dspin_assigned;
          indices = new int [16];
          length = length_assigned;  
          chain = new int [length];
          sigma = sigma_assigned;
          burnin = burnin_assigned;
          store_path = store_path_assigned;
          std::cout << "chain with dspin " << dspin << " built"  << std::endl;
         };        

   ~Chain(){      
          delete [] chain;
          delete [] indices;
          std::cout << "chain with dspin " << dspin << " destroyed"  << std::endl;
         };

};

void dmc_run(Chain chain);