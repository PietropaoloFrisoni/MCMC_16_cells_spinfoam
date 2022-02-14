#include "mcmc.h"

void Metropolis_Hastings_run(Chain &chain)
{

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

    // to this since the amplitude is computed by default with prop_draw
    // TODO: change
    for (int i = 0; i < chain.BIN_SIZE; i++)
    {

        chain.prop_draw[i] = chain.draw[i];
    }

    double amp = chain.pce_amplitude_c16();

    chain.draw[16] = 1;

    if (chain.verbosity > 1)
    {
        std::cout << "Printing Cx coefficients" << std::endl;
        chain.trunc_coeff_print(chain.Ct, chain.dspin);

        std::cout << "Initial draw is:" << std::endl;
        chain.draw_print(chain.draw);

        std::cout << "Initial amplitude is:" << std::endl;
        chain.ampl_print(&amp);
    }

    chain.molteplicity = 1;
    chain.acceptance_ratio = 0;

    double prop_amp = 0;

    double p = 0;

    double z = 0;

    // start moving

    for (int step = 0; step < chain.length; step++)
    {

        chain.RW_monitor = true;

        double draw_double_sample;

        double Cx = 1.0;
        double Cx_prop = 1.0;

        for (int i = 0; i < chain.BIN_SIZE; i++)
        {
            do
            {
                // sample from a GAUSSIAN with mu = 0 and sigma = D
                draw_double_sample = gsl_ran_gaussian_ziggurat(ran, chain.sigma);

                chain.gaussian_draw[i] = 2 * (int)round(draw_double_sample);
                chain.prop_draw[i] = chain.draw[i] + chain.gaussian_draw[i];

            } while (!(0 <= chain.prop_draw[i] && chain.prop_draw[i] <= chain.ti_max));

            if (chain.gaussian_draw[i] != 0)
            {
                chain.RW_monitor = false;
            }

            Cx *= chain.Ct[1][chain.draw[i]];
            Cx_prop *= chain.Ct[1][chain.prop_draw[i]];
        }

        if (chain.verbosity > 1)
        {
            std::cout << "----------------------------------------" << std::endl;

            std::cout << "current draw is:" << std::endl;
            chain.draw_print(chain.draw);

            std::cout << "gaussian draw is:" << std::endl;
            chain.draw_print(chain.gaussian_draw);

            std::cout << "proposed draw is:" << std::endl;
            chain.draw_print(chain.prop_draw);

            std::cout << "----------------------------------------" << std::endl;
        }

        if (chain.RW_monitor == false)
        {
            prop_amp = chain.pce_amplitude_c16();

            p = fmin(1.0, (pow(prop_amp, 2) / pow(amp, 2)) * (Cx / Cx_prop));

            if (isnan(p))
            {
                error("got NaN while computing densities ratio")
            }

            // move or stay
            z = gsl_rng_uniform(ran);

            if (z < p) // accept
            {
                if (chain.verbosity > 1)
                {
                    std::cout << "prop draw accepted as prop amp is " << prop_amp << " and current amp is " << amp << std::endl;
                }

                // CONTINUA DA QUI
            }
            else // reject
            {
            }
        }

        else
        {
            if (chain.verbosity > 1)
            {
                std::cout << "\nThe prop_draw is equal to the current draw, so the molteplicity of the current draw is raised to "
                          << chain.molteplicity + 1 << std::endl;
            }

            chain.acceptance_ratio += 1;
            chain.molteplicity += 1;
            chain.draw[16] += 1;

            continue;
        }
    }
}
