#pragma once

#include <iostream>
#include <string>
#include <vector>

class Chain
{

public:
  // set dimensionality
  static constexpr int BIN_SIZE = 16;

  // coefficient for truncated proposal
  double **Ct;

  int ti_max;
  int i_max;
  int dim_intertw_space;

  // containers of draws == indices, with space for molteplicity
  int draw[BIN_SIZE + 1];
  int prop_draw[BIN_SIZE + 1];
  int gaussian_draw[BIN_SIZE];

  // old safe style (dynamic arrays)
  double *collected_amplitudes;
  int **collected_draws;
  int dspin;
  int length;
  double sigma;
  int burnin;
  std::string store_path;
  int verbosity;

  // accepted moves
  int accepted_moves = 0;

  Chain(int dspin_assigned, int length_assigned, double sigma_assigned, int burnin_assigned, std::string store_path_assigned, int verbosity_assigned)
      : dspin(dspin_assigned), length(length_assigned), sigma(sigma_assigned), burnin(burnin_assigned), store_path(store_path_assigned), verbosity(verbosity_assigned)
  {

    ti_max = 2*dspin;
    i_max = 0.5*ti_max;
    dim_intertw_space = (ti_max - 0) / 2 + 1;

    collected_draws = new int *[length];

    // This allocates several little chunks of memory instead of one single block
    // We pre-allocated more than the required memory, so we don't need to re-allocate
    for (int i = 0; i < length; i++)
    {
      // 16 indices + integer multiplicity of the draw
      collected_draws[i] = new int[BIN_SIZE + 1];
    }

    collected_amplitudes = new double[length];

    // set initial molteplicity to 1
    draw[16] = 1;
  };

  // Prints truncated coefficients
  static inline void trunc_coeff_print(double **Ct, int dspin)
  {
    for (int i = 0; i < BIN_SIZE; i++)
    {
      for (int tk = 0; tk <= 2 * dspin; tk += 2)
      {
        std::cout << "Ct[" << i << "][(" << tk << " - " << 0 << ")/ 2] = " << Ct[i][(tk - 0) / 2] << std::endl;
      }
    }
  }

  // Prints a draw
  static inline void draw_print(int *draw)
  {
    for (size_t i = 0; i < BIN_SIZE; i++)
    {
      std::cout << draw[i] << ' ';
    }
    std::cout << "\t" << draw[16] << std::endl;
  }

  // Prints the amplitude
  static inline void ampl_print(double *ampl)
  {
    std::cout << *ampl << std::endl;
  }

  // Prints all draws and multeplicity
  static inline void print_draws(int **matrix, int rows)
  {
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < BIN_SIZE; j++)
      {
        std::cout << matrix[i][j] << " ";
      }

      std::cout << "\t" << matrix[i][16] << std::endl;
    }
  }

  ~Chain()
  {

    for (int i = 0; i < length; i++)
    {
      delete[] collected_draws[i];
    }

    delete[] collected_draws;

    delete[] collected_amplitudes;

    std::cout << "chain with dspin " << dspin << " destroyed" << std::endl;
  };
};

void dmc_run(Chain &chain);