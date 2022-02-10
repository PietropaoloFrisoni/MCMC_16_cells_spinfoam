#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "fastwigxj.h"
#include "wigxjpf.h"
#include "mcmc.h"
#include "jsymbols.h"
#include "utilities.h"

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

    int C[9] = {};

    double sym = Wigner_21j_symbol(C, chain);

    printf("21J  %#25.15f\n", sym);
}