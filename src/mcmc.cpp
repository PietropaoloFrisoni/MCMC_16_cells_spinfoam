#include "mcmc.h"

void Metropolis_Hastings_run(Chain& chain)
{

    // intialize global seed for random()
    srandom(time(NULL));

    // initializes the PRNG
    gsl_rng *ran;
    ran = gsl_rng_alloc(gsl_rng_taus2);

    // without "+ chain.chain_id" all the chains are exactly the same
    gsl_rng_set(ran, (uint64_t)(random() + chain.chain_id)); 

    double rd;
    double r0, interval;

    for (int i = 0; i < chain.BIN_SIZE; i++)
    {
        interval = (double)(2 * chain.dspin - 0.0) / 2;
        rd = gsl_rng_uniform(ran);
        r0 = rd * interval;
        chain.draw[i] = 2 * ((int)round(r0));
    }

    // to this since the amplitude is computed by default with prop_draw
    // TODO: change
    for (int i = 0; i < chain.BIN_SIZE; i++)
    {

        chain.prop_draw[i] = chain.draw[i];
    }

    double amp = chain.spinfoam_16_cell_amplitude();

    if (chain.verbosity >= 2)
    {
        std::cout << "Printing Cx coefficients" << std::endl;
        chain.trunc_coeff_print(chain.Ct, chain.dspin);

        std::cout << "Initial draw is:" << std::endl;
        chain.draw_print(chain.draw);

        std::cout << "Initial amplitude is:" << std::endl;
        chain.ampl_print(&amp);
    }

    chain.accepted_draws = 0;
    chain.molteplicity = 1;
    chain.acceptance_ratio = 0;

    double prop_amp = 0;
    double p = 0;
    double z = 0;

    // for measuring timings
    struct timespec start, stop;
    double dtime;

    clock_gettime(CLOCK_MONOTONIC, &start);

    progressbar bar(chain.length);
    bar.set_done_char("█");

    // start moving

    for (int step = 0; step <= chain.length; step++)
    {

        if (chain.verbosity >= 1)
        {
            bar.update();
        }

        chain.RW_monitor = true;

        double draw_double_sample;

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
        }

        if (chain.verbosity >= 2)
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
            prop_amp = chain.spinfoam_16_cell_amplitude();

            size_t x_ind, xn_ind;
            double Cx, Cx_new;
            Cx = Cx_new = 1.0;

            for (int i = 0; i < chain.BIN_SIZE; i++)
            {

                x_ind = (chain.draw[i]) / 2;
                xn_ind = (chain.prop_draw[i]) / 2;

                Cx *= chain.Ct[i][x_ind];
                Cx_new *= chain.Ct[i][xn_ind];
            }

            p = fmin(1.0, (pow(prop_amp, 2) / pow(amp, 2)) * (Cx / Cx_new));

            if (isnan(p))
            {
                error("got NaN while computing densities ratio")
            }

            // move or stay
            z = gsl_rng_uniform(ran);

            if (z < p) // accept
            {
                if (chain.verbosity >= 2)
                {
                    std::cout << "prop draw accepted as prop amp is " << prop_amp << " and current amp is " << amp << std::endl;
                    std::cout << "z is " << z << " and p is " << p << std::endl;
                }

                if (step > chain.burnin)
                {

                    for (int i = 0; i < chain.BIN_SIZE; i++)
                    {
                        chain.collected_draws[chain.accepted_draws][i] = chain.draw[i];
                    }

                    chain.collected_draws[chain.accepted_draws][chain.BIN_SIZE] = chain.molteplicity;

                    chain.collected_amplitudes[chain.accepted_draws] = amp;

                    chain.accepted_draws += 1;

                    if (chain.verbosity >= 2)
                    {
                        std::cout << "\nThe old draw:" << std::endl;
                        chain.draw_print(chain.draw);
                        std::cout << "has been stored with molteplicity " << chain.molteplicity << std::endl;
                        std::cout << "The corresponding amplitude " << amp << " has been stored as well" << std::endl;
                        std::cout << "There are " << chain.accepted_draws << " draws stored so far" << std::endl;
                    }
                }

                chain.molteplicity = 1;

                for (int i = 0; i < chain.BIN_SIZE; i++)
                {
                    chain.draw[i] = chain.prop_draw[i];
                }

                amp = prop_amp;
                chain.acceptance_ratio += 1;

                if (chain.verbosity >= 2)
                {
                    std::cout << "Now the new draw is:" << std::endl;
                    chain.draw_print(chain.draw);
                    std::cout << "the new amp is:" << amp << std::endl;
                }
            }
            else // reject
            {

                chain.molteplicity += 1;

                if (chain.verbosity >= 2)
                {
                    std::cout << "Prop_draw:" << std::endl;
                    chain.draw_print(chain.prop_draw);
                    std::cout << "was rejected, since p " << p << " and z " << z << std::endl;
                    std::cout << "The current draw:" << std::endl;
                    chain.draw_print(chain.draw);
                    std::cout << "remains the same and its molteplicity is " << chain.molteplicity << "\namp remains " << amp << std::endl;
                }
            }
        }

        else
        {
            if (chain.verbosity >= 2)
            {
                std::cout << "\nThe prop_draw is equal to the current draw, so the molteplicity of the current draw is raised to "
                          << chain.molteplicity + 1 << std::endl;
            }

            chain.acceptance_ratio += 1;
            chain.molteplicity += 1;

            continue;
        }

        // final storage
        if (step == chain.length)
        {

            for (int i = 0; i < chain.BIN_SIZE; i++)
            {
                chain.collected_draws[chain.accepted_draws][i] = chain.draw[i];
            }

            chain.collected_draws[chain.accepted_draws][chain.BIN_SIZE] = chain.molteplicity;

            chain.collected_amplitudes[chain.accepted_draws] = amp;

            chain.accepted_draws += 1;

            if (chain.verbosity >= 2)
            {
                std::cout << "\n\n\nThe final draw:" << std::endl;
                chain.draw_print(chain.draw);
                std::cout << "has been stored with molteplicity " << chain.molteplicity << std::endl;
                std::cout << "The corresponding amplitude " << amp << " has been stored as well" << std::endl;
                std::cout << "There are " << chain.accepted_draws << " draws stored in total" << std::endl;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &stop);
    dtime = (double)(stop.tv_sec - start.tv_sec);
    dtime += (double)(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
    chain.run_time = dtime;

    if (chain.verbosity >= 1)
    {
        printf("\n\nChain completed.\nTime elapsed: %f seconds.\n", dtime);
        chain.print_statistics();
    }

    chain.store_draws();
}
