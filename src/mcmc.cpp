#include "mcmc.h"

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

    std::cout << "Starting sampling chain with parameters:\n\ndspin = " << chain.dspin << "\nlength = " << chain.length
              << "\nsigma = " << chain.sigma << "\nstore path = " << chain.store_path << std::endl;

    // intialize global seed for random()
    // TODO move this inside class constructor
    srandom(time(NULL));

    // initializes the PRNG
    // TODO move this inside class constructor
    gsl_rng *ran;
    ran = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(ran, (uint64_t)(random()));

    double rd;
    double r0, interval;

    for (int i = 0; i < chain.BIN_SIZE; i++)
    {
        interval = (double)(2 * chain.dspin - 0.0) / 2;
        rd = gsl_rng_uniform(ran);
        r0 = rd * interval;
        chain.draw[i] = 2 * ((int)round(r0));
    }

    // precompute the coefficients for truncated proposals
    // TODO move this inside class constructor
    chain.Ct = (double **)malloc(chain.BIN_SIZE * sizeof(double *));
    for (int i = 0; i < chain.BIN_SIZE; i++)
    {
        double *Cti = (double *)malloc(chain.dim_intertw_space * sizeof(double));
        chain.Ct[i] = Cti;

        double Cxk;
        int k;
        for (int tk = 0; tk <= chain.ti_max; tk += 2)
        {
            k = tk * 0.5;
            Cxk = cdf_gaussian_discrete(0 - k, chain.i_max - k, chain.sigma);
            Cti[(tk - 0) / 2] = Cxk;
        }
    }

    double sym = chain.pce_amplitude_c16();

    if (chain.verbosity > 1)
    {
        std::cout << "Printing Cx coefficients" << std::endl;
        chain.trunc_coeff_print(chain.Ct, chain.dspin);

        std::cout << "Initial draw is:" << std::endl;
        chain.draw_print(chain.draw);

        std::cout << "Initial amplitude is:" << std::endl;
        chain.ampl_print(&sym);
    }
}