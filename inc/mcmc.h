#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <filesystem>
#include <getopt.h>
#include <omp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "phmap.h"
#include "phmap_dump.h"
#include "common.h"
#include "hash_21j_symbols.h"
#include "progressbar.h"

// the PDF of a gaussian rounded to the integers
inline double pdf_gaussian_discrete(const uint8_t&, const double&);

// the CDF of a gaussian rounded to the integers
inline double cdf_gaussian_discrete(const uint8_t, const uint8_t, const double&);

class Chain
{

private:
  // boundary spins
  uint8_t ti_1, ti_2, ti_3, ti_4,
      ti_5, ti_6, ti_7, ti_8,
      ti_9, ti_10, ti_11, ti_12,
      ti_13, ti_14, ti_15, ti_16;

  // bulk virtual spins
  uint8_t tb_1, tb_2, tb_3, tb_4, tb_5; // blue spins
  uint8_t tpn_1, tpn_2, tps_1, tps_2;   // purple spins

  // bounds
  uint8_t tb1_min, tb1_max;
  uint8_t tb2_min, tb2_max;
  uint8_t tb3_min, tb3_max;
  uint8_t tb4_min, tb4_max;
  uint8_t tb5_min, tb5_max;
  uint8_t tpn1_min, tpn1_max;
  uint8_t tpn2_min, tpn2_max;
  uint8_t tps1_min, tps1_max;
  uint8_t tps2_min, tps2_max;

  // path in which draws will be stored
  std::string store_path;

  // path in which hashed tables of 21j symbols are stored
  std::string hashed_tables_path;

  // amplitude with compensated summation
  double ampl, c, y, t;

  double df, gdf, ph;

  // build half NORTH
  double aN, cN, yN, tN;
  double aNW, aNE;

  // build half SOUTH
  double aS, cS, yS, tS;
  double aSW, aSE;

  // used as path shortcut
  char tmp[1024];

public:
  using Hash = phmap::flat_hash_map<MyKey, double>;
  std::shared_ptr<Hash> h;

  // set dimensionality
  static constexpr int BIN_SIZE = 16;

  // test if the RW actually moved
  bool RW_monitor;

  // coefficient for truncated proposal
  double **Ct;

  int const ti_max;
  int const i_max;
  int const dim_intertw_space;

  // containers for draws, with space at the end for multeplicity
  int draw[BIN_SIZE];
  int prop_draw[BIN_SIZE];
  int gaussian_draw[BIN_SIZE];

  double *collected_amplitudes;
  int **collected_draws;
  int const dspin;
  int const length;
  double const sigma;
  int const burnin;
  int const verbosity;
  int const chain_id;

  // accepted moves
  int acceptance_ratio;
  int accepted_draws;
  int molteplicity;
  float run_time;

  // constructor for the Markov chain
  Chain(std::string store_path_assigned, std::string hashed_tables_path_assigned, const int dspin_assigned,
        const int length_assigned, const double sigma_assigned, const int burnin_assigned,
        const int verbosity_assigned, const int thread_id_assigned, const std::shared_ptr<Hash> hash);

  // Prints truncated coefficients
  static inline void trunc_coeff_print(double **Ct, int dspin);

  // Prints a generic given draw
  void draw_print(int *draw);

  // Prints a generic given amplitude
  void ampl_print(double *ampl);

  // Prints all draws and multeplicity
  void print_collected_draws();

  void print_statistics();

  void store_draws();

  double spinfoam_16_cell_amplitude();

  ~Chain();
  
};

void Metropolis_Hastings_run(Chain &chain);
