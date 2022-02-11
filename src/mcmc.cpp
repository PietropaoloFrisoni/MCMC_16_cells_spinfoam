#include <iostream>
#include <typeinfo>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "fastwigxj.h"
#include "wigxjpf.h"

#include "mcmc.h"
#include "jsymbols.h"
#include "utilities.h"
#include "common.h"

// the PDF of a gaussian rounded to the integers
static double pdf_gaussian_discrete(int n, double s)
{

    return gsl_cdf_gaussian_P((double)n + 0.5, s) -
           gsl_cdf_gaussian_P((double)n - 0.5, s);
}

// the CDF of a gaussian rounded to the integers
static double cdf_gaussian_discrete(int n1, int n2, double s)
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

void dmc_run(Chain &chain)
{

    std::cout << "Starting sampling chain with parameters:\n\nspin = " << chain.dspin << "\nlength = " << chain.length
              << "\nsigma = " << chain.sigma << "\nstore path = " << chain.store_path << std::endl;

    print_matrix(chain.draws, chain.length, 17);

    // intialize global seed for random()
    srandom(time(NULL));

    // set dimensionality
    const int BIN_SIZE = 16;

    // accepted moves
    int accepted_moves = 0;

    // initializes the PRNG
    gsl_rng *ran;
    ran = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(ran, (uint64_t)(random()));

    // precompute the coefficients for truncated proposals
    double **Ct = (double **)malloc(BIN_SIZE * sizeof(double *));
    size_t dimi;
    int ti_min, ti_max;
    int i_min, i_max;
    float Di;

    for (int i = 0; i < BIN_SIZE; i++)
    {

        Di = chain.sigma;
        ti_min = 0;
        ti_max = 2 * chain.dspin;
        i_min = ti_min * 0.5;
        i_max = ti_max * 0.5;

        dimi = (ti_max - ti_min) / 2 + 1;
        double *Cti = (double *)malloc(dimi * sizeof(double));
        Ct[i] = Cti;

        double Cxk;
        int k;
        for (int tk = ti_min; tk <= ti_max; tk += 2)
        {

            k = tk * 0.5;
            Cxk = cdf_gaussian_discrete(i_min - k, i_max - k, Di);
            Cti[(tk - ti_min) / 2] = Cxk;
        }
    }

    printf("Di = %f\n", Di);
    std::cout << "Printing Cx coefficients" << std::endl;
    for (int q = 0; q < BIN_SIZE; q++)
    {
        for (int tk = ti_min; tk <= ti_max; tk += 2)
        {
          std::cout << "Ct[" << q << "][(" << tk << " - " << ti_min << ")/ 2] = " << Ct[q][(tk - ti_min) / 2] << std::endl;
        }
    }
}