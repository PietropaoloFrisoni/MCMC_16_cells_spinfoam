#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <filesystem>
#include <getopt.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "phmap.h"
#include "phmap_dump.h"
#include "common.h"
#include "hash_21j_symbols.h"
#include "progressbar.h"

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

  Chain(std::string store_path_assigned, std::string hashed_tables_path_assigned, const int dspin_assigned,
        const int length_assigned, const double sigma_assigned, const int burnin_assigned,
        const int verbosity_assigned, const int thread_id_assigned)
      : store_path(store_path_assigned), hashed_tables_path(hashed_tables_path_assigned), dspin(dspin_assigned),
        length(length_assigned), sigma(sigma_assigned), burnin(burnin_assigned),
        verbosity(verbosity_assigned), chain_id(thread_id_assigned)
  {

    sprintf(tmp, "j_%.8g", ((double)(dspin) / 2.0));

    hashed_tables_path_assigned = hashed_tables_path_assigned + "/hashed_21j_symbols_" + std::string(tmp);

    phmap::BinaryInputArchive ar_in(&hashed_tables_path_assigned[0]);
    h.phmap_load(ar_in);

    ti_max = 2 * dspin;
    i_max = 0.5 * ti_max;
    dim_intertw_space = (ti_max - 0) / 2 + 1;

    collected_draws = new int *[length];

    // This allocates several little chunks of memory instead of one single block
    for (int i = 0; i < length; i++)
    {
      // 16 indices + integer multiplicity of the draw
      collected_draws[i] = new int[BIN_SIZE + 1];
    }

    collected_amplitudes = new double[length];

    tb1_min = tb5_min = 0;
    tb1_max = tb5_max = 2 * dspin;

    tpn1_min = 0;
    tpn1_max = 2 * dspin;

    tps1_min = 0;
    tps1_max = 2 * dspin;

    // Coefficients for truncated gaussian proposal
    Ct = (double **)malloc(BIN_SIZE * sizeof(double *));
    for (int i = 0; i < BIN_SIZE; i++)
    {
      double *Cti = (double *)malloc(dim_intertw_space * sizeof(double));
      Ct[i] = Cti;

      double Cxk;
      int k;
      for (int tk = 0; tk <= ti_max; tk += 2)
      {
        k = tk * 0.5;
        Cxk = cdf_gaussian_discrete(0 - k, i_max - k, sigma);
        Cti[(tk - 0) / 2] = Cxk;
      }
    }
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

  // Prints a generic given draw
  static inline void draw_print(int *draw)
  {
    for (int i = 0; i < BIN_SIZE; i++)
    {
      std::cout << draw[i] << ' ';
    }
    std::cout << std::endl;
  }

  // Prints a generic given amplitude
  static inline void ampl_print(double *ampl)
  {
    std::cout << *ampl << std::endl;
  }

  // Prints all draws and multeplicity
  inline void print_collected_draws()
  {
    for (int i = 0; i < accepted_draws; i++)
    {
      for (int j = 0; j < BIN_SIZE; j++)
      {
        std::cout << collected_draws[i][j] << " ";
      }

      std::cout << "\t" << collected_draws[i][BIN_SIZE] << std::endl;
    }
  }

  void print_statistics()
  {

    std::cout << (acceptance_ratio * 100 / length) << "% of draws have been accepted" << std::endl;
  }

  void store_draws()
  {

    char tmp[1024];

    sprintf(tmp, "/j_%.8g", ((double)(dspin) / 2.0));

    store_path += std::string(tmp);

    std::filesystem::create_directories(store_path);

    sprintf(tmp, "/N_%d__sigma_%.8g__burnin_%.d_chain_%.d.csv", length, sigma, burnin, chain_id);

    store_path += std::string(tmp);

    std::ofstream out(store_path);

    for (int i = 0; i < BIN_SIZE; i++)
    {
      out << "intertwiner " << std::to_string(i + 1) << ',';
    }

    out << "draw multeplicity" << ',' << "draw amplitude" << ','
        << "total accept. draws" << ',' << "total accept. rate" << ',' << "total run time" << '\n';

    for (int j = 0; j < BIN_SIZE; j++)
    {
      out << collected_draws[0][j] / 2 << ',';
    }

    out << collected_draws[0][BIN_SIZE] << ',' << collected_amplitudes[0] << ','
        << accepted_draws << ',' << acceptance_ratio << "%," << run_time << " s\n";

    for (int i = 1; i < accepted_draws; i++)
    {
      for (int j = 0; j < BIN_SIZE; j++)
      {
        out << collected_draws[i][j] / 2 << ',';
      }

      out << collected_draws[i][BIN_SIZE] << ',' << collected_amplitudes[i] << '\n';
    }
  }

  double pce_amplitude_c16()
  {
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

              MyKey kNW{ti_1, ti_2, ti_3, ti_4, tb1,
                        tb2, tb3, 0, 0};

              MyKey kNE{ti_16, ti_15, ti_14, ti_13, tb5,
                        tb4, tb3, 0, 0};

              for (tpn1 = tpn1_min; tpn1 <= tpn1_max; tpn1 += 2)
              {

                tpn2_min = max(abs(tpn1 - dspin), abs(tb3 - dspin));
                tpn2_max = min(tpn1 + dspin, tb3 + dspin);

                kNW[7] = kNE[7] = tpn1;

                for (tpn2 = tpn2_min; tpn2 <= tpn2_max; tpn2 += 2)
                {

                  kNW[8] = kNE[8] = tpn2;

                  // NW 21j

                  //     key_21j[0] = (uint8_t)ti_1;
                  //     key_21j[1] = (uint8_t)ti_2;
                  //     key_21j[2] = (uint8_t)ti_3;
                  //     key_21j[3] = (uint8_t)ti_4;
                  //     key_21j[4] = (uint8_t)tb1;
                  //     key_21j[5] = (uint8_t)tb2;
                  //     key_21j[6] = (uint8_t)tb3;
                  //     key_21j[7] = (uint8_t)tpn1;
                  //     key_21j[8] = (uint8_t)tpn2;

                  aNW = h[kNW];

                  // NE 21j
                  // reflect from left

                  //     key_21j[0] = (uint8_t)ti_16;
                  //     key_21j[1] = (uint8_t)ti_15;
                  //     key_21j[2] = (uint8_t)ti_14;
                  //     key_21j[3] = (uint8_t)ti_13;
                  //     key_21j[4] = (uint8_t)tb5;
                  //     key_21j[5] = (uint8_t)tb4;
                  //     key_21j[6] = (uint8_t)tb3;
                  //     key_21j[7] = (uint8_t)tpn1;
                  //     key_21j[8] = (uint8_t)tpn2;

                  aNE = h[kNE];

                  df = DIM(tpn1) * DIM(tpn2);

                  comp_sum(df * aNW * aNE, aN, cN, yN, tN);

                } // tpn2
              }   // tpn1

              // build half SOUTH
              aS = cS = 0;

              MyKey kSW{ti_8, ti_7, ti_6, ti_5, tb1,
                        tb2, tb3, 0, 0};

              MyKey kSE{ti_9, ti_10, ti_11, ti_12, tb5,
                        tb4, tb3, 0, 0};

              for (tps1 = tps1_min; tps1 <= tps1_max; tps1 += 2)
              {

                tps2_min = max(abs(tps1 - dspin), abs(tb3 - dspin));
                tps2_max = min(tps1 + dspin, tb3 + dspin);

                kSW[7] = kSE[7] = tps1;

                for (tps2 = tps2_min; tps2 <= tps2_max; tps2 += 2)
                {

                  kSW[8] = kSE[8] = tps2;

                  // SW 21j

                  //    key_21j[0] = (uint8_t)ti_8;
                  //    key_21j[1] = (uint8_t)ti_7;
                  //    key_21j[2] = (uint8_t)ti_6;
                  //    key_21j[3] = (uint8_t)ti_5;
                  //    key_21j[4] = (uint8_t)tb1;
                  //    key_21j[5] = (uint8_t)tb2;
                  //    key_21j[6] = (uint8_t)tb3;
                  //    key_21j[7] = (uint8_t)tps1;
                  //    key_21j[8] = (uint8_t)tps2;

                  aSW = h[kSW];

                  // SE 21j

                  //    key_21j[0] = (uint8_t)ti_9;
                  //    key_21j[1] = (uint8_t)ti_10;
                  //    key_21j[2] = (uint8_t)ti_11;
                  //    key_21j[3] = (uint8_t)ti_12;
                  //    key_21j[4] = (uint8_t)tb5;
                  //    key_21j[5] = (uint8_t)tb4;
                  //    key_21j[6] = (uint8_t)tb3;
                  //    key_21j[7] = (uint8_t)tps1;
                  //    key_21j[8] = (uint8_t)tps2;

                  aSE = h[kSE];

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

    h.clear();
  };
};

void Metropolis_Hastings_run(Chain &chain);
