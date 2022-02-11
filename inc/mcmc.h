#pragma once

#include <iostream>
#include <string>
#include <vector>

class Chain
{

public:
  // Old and safe style (dynamic arrays)
  double *amplitudes;
  int **draws;

  int dspin;
  int length;
  int indices[16];
  double sigma;
  int burnin;
  std::string store_path;

  Chain(int dspin_assigned, int length_assigned, double sigma_assigned, int burnin_assigned, std::string store_path_assigned)
      : dspin(dspin_assigned), length(length_assigned), sigma(sigma_assigned), burnin(burnin_assigned), store_path(store_path_assigned)
  {

    draws = new int *[length];

    // This allocates several little chunks of memory instead of one single block
    // We pre-allocated more than the required memory, so we don't need to re-allocate
    for (int i = 0; i < length; i++)
    {
      // 16 indices + integer multiplicity of the draw
      draws[i] = new int[17];
    }
    
    amplitudes = new double[length];
  };

  ~Chain()
  {

    for (int i = 0; i < length; i++)
    {
      delete[] draws[i];
    }

    delete[] draws;

    delete[] amplitudes;

    std::cout << "chain with dspin " << dspin << " destroyed" << std::endl;
  };
};

void dmc_run(Chain &chain);