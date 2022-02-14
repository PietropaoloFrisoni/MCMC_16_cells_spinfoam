#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "phmap.h"
#include "phmap_dump.h"
#include "common.h"


// this must be included even if we don't hash because we need "struct MyKey"
// TODO change this
#include "hash_21j_symbols.h"

class Chain
{

private:
  // boundary spins
  uint8_t ti_1, ti_2, ti_3, ti_4,
      ti_5, ti_6, ti_7, ti_8,
      ti_9, ti_10, ti_11, ti_12,
      ti_13, ti_14, ti_15, ti_16;

  // bulk virtual spins
  uint8_t tb1, tb2, tb3, tb4, tb5; // blue spins
  uint8_t tpn1, tpn2, tps1, tps2;  // purple spins

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

  // keys used to recover 21js from hash tables
  // TODO this will be multi-thread
  uint8_t key_21j[9];

  // amplitude with compensated summation
  double ampl, c, y, t;

  double df, gdf, ph;

  // build half NORTH
  double aN, cN, yN, tN;
  double aNW, aNE;

  // build half SOUTH
  double aS, cS, yS, tS;
  double aSW, aSE;

public:
  using Hash = phmap::flat_hash_map<MyKey, double>;
  Hash h{};

  // set dimensionality
  static constexpr int BIN_SIZE = 16;

  // test if the RW actually moved
  bool RW_monitor;

  // coefficient for truncated proposal
  double **Ct;

  int ti_max;
  int i_max;
  int dim_intertw_space;

  // containers for draws, with space at the end for multeplicity
  int draw[BIN_SIZE + 1];
  int prop_draw[BIN_SIZE + 1];
  int gaussian_draw[BIN_SIZE];

  double *collected_amplitudes;
  uint8_t **collected_draws;
  int const dspin;
  int const length;
  double const sigma;
  int const burnin;
  int const verbosity;

  // accepted moves
  int acceptance_ratio;

  int molteplicity;

  Chain(std::string store_path_assigned, std::string hashed_tables_path_assigned, const int dspin_assigned,
        const int length_assigned, const double sigma_assigned, const int burnin_assigned, const int verbosity_assigned)
      : store_path(store_path_assigned), hashed_tables_path(hashed_tables_path_assigned), dspin(dspin_assigned),
        length(length_assigned), sigma(sigma_assigned), burnin(burnin_assigned), verbosity(verbosity_assigned)
  {
    hashed_tables_path_assigned = hashed_tables_path_assigned + "/Hashed_21j_symbols_dspin_" + std::to_string(dspin);

    phmap::BinaryInputArchive ar_in(&hashed_tables_path_assigned[0]);
    h.phmap_load(ar_in);

    ti_max = 2 * dspin;
    i_max = 0.5 * ti_max;
    dim_intertw_space = (ti_max - 0) / 2 + 1;

    collected_draws = new uint8_t *[length];

    // This allocates several little chunks of memory instead of one single block
    // We pre-allocated more than the required memory, so we don't need to re-allocate
    for (int i = 0; i < length; i++)
    {
      // 16 indices + integer multiplicity of the draw
      collected_draws[i] = new uint8_t[BIN_SIZE + 1];
    }

    collected_amplitudes = new double[length];

    tb1_min = tb5_min = 0;
    tb1_max = tb5_max = 2 * dspin;

    tpn1_min = 0;
    tpn1_max = 2 * dspin;

    tps1_min = 0;
    tps1_max = 2 * dspin;
  };

  // Prints truncated coefficients
  static inline void trunc_coeff_print(double **Ct, int dspin)
  {
    for (int i = 0; i < BIN_SIZE; i++)
    {
      for (int tk = 0; tk <= 2 * dspin; tk += 2)
      {
        std::cout << "Ct[" << i << "][(" << tk << " - " << 0 << ")/ 2] = " << Ct[i][(tk - 0) / 2] << "\t";
      }

      std::cout << std::endl;
    }
  }

  // Prints a draw
  static inline void draw_print(int *draw)
  {
    for (int i = 0; i < BIN_SIZE; i++)
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
  static inline void print_draws(uint8_t **matrix, int rows)
  {
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < BIN_SIZE; j++)
      {
        std::cout << unsigned(matrix[i][j]) << " ";
      }

      std::cout << "\t" << unsigned(matrix[i][16]) << std::endl;
    }
  }

  double pce_amplitude_c16()
  {

    // TODO: remove all castings to uint8_t

    // boundary data
    // spins are to be read counterclockwise
    // starting from top
    ti_1 = prop_draw[0];
    ti_2 = prop_draw[1];
    ti_3 = prop_draw[2];
    ti_4 = prop_draw[3];

    ti_5 = prop_draw[4];
    ti_6 = prop_draw[5];
    ti_7 = prop_draw[6];
    ti_8 = prop_draw[7];

    ti_9 = prop_draw[8];
    ti_10 = prop_draw[9];
    ti_11 = prop_draw[10];
    ti_12 = prop_draw[11];

    ti_13 = prop_draw[12];
    ti_14 = prop_draw[13];
    ti_15 = prop_draw[14];
    ti_16 = prop_draw[15];

    // I have to assemble the two halves, the top one and the bottom one
    // I sum the quartes NW and NE over purple spins tpn1, tpn2
    // to make the top half
    // I sum the quartes SW and SE over purple spins tps1, tps2
    // to make the bottom half
    // then I sum the two halves over the 5 blue spins tb1..5

    ampl = c = 0;

    for (tb1 = tb1_min; tb1 <= tb1_max; tb1 += 2)
    {

      tb2_min = abs(tb1 - dspin);
      tb2_max = tb1 + dspin;

      for (tb5 = tb5_min; tb5 <= tb5_max; tb5 += 2)
      {

        tb4_min = abs(tb5 - dspin);
        tb4_max = tb5 + dspin;

        for (tb2 = tb2_min; tb2 <= tb2_max; tb2 += 2)
        {
          for (tb4 = tb4_min; tb4 <= tb4_max; tb4 += 2)
          {

            tb3_min = max(abs(tb2 - dspin), abs(tb4 - dspin));
            tb3_max = min(tb2 + dspin, tb4 + dspin);

            for (tb3 = tb3_min; tb3 <= tb3_max; tb3 += 2)
            {

              // build half NORTH
              aN = cN = 0;

              for (tpn1 = tpn1_min; tpn1 <= tpn1_max; tpn1 += 2)
              {

                tpn2_min = max(abs(tpn1 - dspin), abs(tb3 - dspin));
                tpn2_max = min(tpn1 + dspin, tb3 + dspin);

                for (tpn2 = tpn2_min; tpn2 <= tpn2_max; tpn2 += 2)
                {

                  // NW 21j

                  key_21j[0] = (uint8_t)ti_1;
                  key_21j[1] = (uint8_t)ti_2;
                  key_21j[2] = (uint8_t)ti_3;
                  key_21j[3] = (uint8_t)ti_4;
                  key_21j[4] = (uint8_t)tb1;
                  key_21j[5] = (uint8_t)tb2;
                  key_21j[6] = (uint8_t)tb3;
                  key_21j[7] = (uint8_t)tpn1;
                  key_21j[8] = (uint8_t)tpn2;

                  aNW = h[MyKey{key_21j[0], key_21j[1], key_21j[2], key_21j[3], key_21j[4],
                                key_21j[5], key_21j[6], key_21j[7], key_21j[8]}];

                  // NE 21j
                  // reflect from left

                  key_21j[0] = (uint8_t)ti_16;
                  key_21j[1] = (uint8_t)ti_15;
                  key_21j[2] = (uint8_t)ti_14;
                  key_21j[3] = (uint8_t)ti_13;
                  key_21j[4] = (uint8_t)tb5;
                  key_21j[5] = (uint8_t)tb4;
                  key_21j[6] = (uint8_t)tb3;
                  key_21j[7] = (uint8_t)tpn1;
                  key_21j[8] = (uint8_t)tpn2;

                  aNE = h[MyKey{key_21j[0], key_21j[1], key_21j[2], key_21j[3], key_21j[4],
                                key_21j[5], key_21j[6], key_21j[7], key_21j[8]}];

                  df = DIM(tpn1) * DIM(tpn2);

                  comp_sum(df * aNW * aNE, aN, cN, yN, tN);

                } // tpn2
              }   // tpn1

              // build half SOUTH
              aS = cS = 0;

              for (tps1 = tps1_min; tps1 <= tps1_max; tps1 += 2)
              {

                tps2_min = max(abs(tps1 - dspin), abs(tb3 - dspin));
                tps2_max = min(tps1 + dspin, tb3 + dspin);

                for (tps2 = tps2_min; tps2 <= tps2_max; tps2 += 2)
                {

                  // SW 21j

                  key_21j[0] = (uint8_t)ti_8;
                  key_21j[1] = (uint8_t)ti_7;
                  key_21j[2] = (uint8_t)ti_6;
                  key_21j[3] = (uint8_t)ti_5;
                  key_21j[4] = (uint8_t)tb1;
                  key_21j[5] = (uint8_t)tb2;
                  key_21j[6] = (uint8_t)tb3;
                  key_21j[7] = (uint8_t)tps1;
                  key_21j[8] = (uint8_t)tps2;

                  aSW = h[MyKey{key_21j[0], key_21j[1], key_21j[2], key_21j[3], key_21j[4],
                                key_21j[5], key_21j[6], key_21j[7], key_21j[8]}];

                  // SE 21j

                  key_21j[0] = (uint8_t)ti_9;
                  key_21j[1] = (uint8_t)ti_10;
                  key_21j[2] = (uint8_t)ti_11;
                  key_21j[3] = (uint8_t)ti_12;
                  key_21j[4] = (uint8_t)tb5;
                  key_21j[5] = (uint8_t)tb4;
                  key_21j[6] = (uint8_t)tb3;
                  key_21j[7] = (uint8_t)tps1;
                  key_21j[8] = (uint8_t)tps2;

                  aSE = h[MyKey{key_21j[0], key_21j[1], key_21j[2], key_21j[3], key_21j[4],
                                key_21j[5], key_21j[6], key_21j[7], key_21j[8]}];

                  df = DIM(tps1) * DIM(tps2);

                  // phase from reflecting back to stored 21j
                  // ph = real_negpow(
                  //         (dspin + dspin + tb1) + (tb1 + dspin + tb2) + (tb2 + dspin + tb3) +   // SW 21j
                  //         (dspin + dspin + tb1) + (tb1 + dspin + tb2) + (tb2 + dspin + tb3) +   // SE 21j ...
                  //         (dspin + dspin + tps1) + (tps1 + dspin + tps2) + (tps2 + dspin + tb3) //
                  // );

                  // simplified
                  ph = real_negpow(2 * tps1 + 2 * tps2 + 3 * tb3);

                  comp_sum(ph * df * aSW * aSE, aS, cS, yS, tS);

                } // tps2
              }   // tps1

              df = DIM(tb1) * DIM(tb2) * DIM(tb3) * DIM(tb4) * DIM(tb5);

              // two halves computed, assemble
              comp_sum(ph * df * aN * aS, ampl, c, y, t);

            } // tb3
          }   // tb4
        }     // tb2
      }       // tb5
    }         // tb1

    // global dimensional factors
    gdf = sqrt((double)DIM(ti_1) *
               (double)DIM(ti_2) *
               (double)DIM(ti_3) *
               (double)DIM(ti_4) *
               (double)DIM(ti_5) *
               (double)DIM(ti_6) *
               (double)DIM(ti_7) *
               (double)DIM(ti_8) *
               (double)DIM(ti_9) *
               (double)DIM(ti_10) *
               (double)DIM(ti_11) *
               (double)DIM(ti_12) *
               (double)DIM(ti_13) *
               (double)DIM(ti_14) *
               (double)DIM(ti_15) *
               (double)DIM(ti_16));

    return gdf * ampl;
  }

  ~Chain()
  {

    for (int i = 0; i < length; i++)
    {
      delete[] collected_draws[i];
    }

    delete[] collected_draws;

    delete[] collected_amplitudes;

    // free the coefficients for truncated proposals
    for (int i = 0; i < BIN_SIZE; i++)
    {
      delete[] Ct[i];
    }
    delete[] Ct;

    std::cout << "chain with dspin " << dspin << " destroyed" << std::endl;
  };
};

// the PDF of a gaussian rounded to the integers
static inline double pdf_gaussian_discrete(int n, double s)
{

  return gsl_cdf_gaussian_P((double)n + 0.5, s) -
         gsl_cdf_gaussian_P((double)n - 0.5, s);
}

// the CDF of a gaussian rounded to the integers
static inline double cdf_gaussian_discrete(int n1, int n2, double s)
{

  double r = 0.0;
  int x = n1;
  while (x <= n2)
  {
    r += pdf_gaussian_discrete(x, s);
    x++;
  }

  return r;
}

void Metropolis_Hastings_run(Chain &chain);
