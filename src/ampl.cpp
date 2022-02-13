#include <iostream>

#include "utilities.h"
#include "ampl.h"
#include "common.h"
#include "jsymbols.h"


double pce_amplitude_c16(Chain &chain)
{

    // TODO TUTTO ANCORA DA OTTIMIZZARE

    int key_21j[9];

    // boundary data
    // spins are to be read counterclockwise
    // starting from top
    int ti_1 = chain.draw[0];
    int ti_2 = chain.draw[1];
    int ti_3 = chain.draw[2];
    int ti_4 = chain.draw[3];

    int ti_5 = chain.draw[4];
    int ti_6 = chain.draw[5];
    int ti_7 = chain.draw[6];
    int ti_8 = chain.draw[7];

    int ti_9 = chain.draw[8];
    int ti_10 = chain.draw[9];
    int ti_11 = chain.draw[10];
    int ti_12 = chain.draw[11];

    int ti_13 = chain.draw[12];
    int ti_14 = chain.draw[13];
    int ti_15 = chain.draw[14];
    int ti_16 = chain.draw[15];

    // I have to assemble the two halves, the top one and the bottom one
    // I sum the quartes NW and NE over purple spins tpn1, tpn2
    // to make the top half
    // I sum the quartes SW and SE over purple spins tps1, tps2
    // to make the bottom half
    // then I sum the two halves over the 5 blue spins tb1..5

    // bulk virtual spins
    int tb1, tb2, tb3, tb4, tb5; // blue spins
    int tpn1, tpn2, tps1, tps2;  // purple spins

    // bounds
    int tb1_min, tb1_max;
    int tb2_min, tb2_max;
    int tb3_min, tb3_max;
    int tb4_min, tb4_max;
    int tb5_min, tb5_max;
    int tpn1_min, tpn1_max;
    int tpn2_min, tpn2_max;
    int tps1_min, tps1_max;
    int tps2_min, tps2_max;

    // amplitude with compensated summation
    double ampl, c, y, t;
    ampl = c = 0;

    double df, gdf, ph;

    tb1_min = tb5_min = 0;
    tb1_max = tb5_max = 2 * chain.dspin;

    // khint_t s;
    // HashTable21j_key_t key;

    for (tb1 = tb1_min; tb1 <= tb1_max; tb1 += 2)
    {

        tb2_min = abs(tb1 - chain.dspin);
        tb2_max = tb1 + chain.dspin;

        for (tb5 = tb5_min; tb5 <= tb5_max; tb5 += 2)
        {

            tb4_min = abs(tb5 - chain.dspin);
            tb4_max = tb5 + chain.dspin;

            for (tb2 = tb2_min; tb2 <= tb2_max; tb2 += 2)
            {
                for (tb4 = tb4_min; tb4 <= tb4_max; tb4 += 2)
                {

                    tb3_min = max(abs(tb2 - chain.dspin), abs(tb4 - chain.dspin));
                    tb3_max = min(tb2 + chain.dspin, tb4 + chain.dspin);

                    for (tb3 = tb3_min; tb3 <= tb3_max; tb3 += 2)
                    {

                        // build half NORTH
                        double aN, cN, yN, tN;
                        aN = cN = 0;

                        double aNW, aNE;

                        tpn1_min = 0;
                        tpn1_max = 2 * chain.dspin;

                        for (tpn1 = tpn1_min; tpn1 <= tpn1_max; tpn1 += 2)
                        {

                            tpn2_min = max(abs(tpn1 - chain.dspin), abs(tb3 - chain.dspin));
                            tpn2_max = min(tpn1 + chain.dspin, tb3 + chain.dspin);

                            for (tpn2 = tpn2_min; tpn2 <= tpn2_max; tpn2 += 2)
                            {

                                // NW 21j

                                key_21j[0] = ti_1;
                                key_21j[1] = ti_2;
                                key_21j[2] = ti_3;
                                key_21j[3] = ti_4;
                                key_21j[4] = tb1;
                                key_21j[5] = tb2;
                                key_21j[6] = tb3;
                                key_21j[7] = tpn1;
                                key_21j[8] = tpn2;

                                aNW = Wigner_21j_symbol(key_21j, chain); // CHECK MEMORY ALLOCATION

                                // NE 21j
                                // reflect from left

                                key_21j[0] = ti_16;
                                key_21j[1] = ti_15;
                                key_21j[2] = ti_14;
                                key_21j[3] = ti_13;
                                key_21j[4] = tb5;
                                key_21j[5] = tb4;
                                key_21j[6] = tb3;
                                key_21j[7] = tpn1;
                                key_21j[8] = tpn2;

                                aNE = Wigner_21j_symbol(key_21j, chain);

                                df = DIM(tpn1) * DIM(tpn2);

                                comp_sum(df * aNW * aNE, aN, cN, yN, tN);

                            } // tpn2
                        }     // tpn1

                        // build half SOUTH
                        double aS, cS, yS, tS;
                        aS = cS = 0;

                        double aSW, aSE;

                        tps1_min = 0;
                        tps1_max = 2 * chain.dspin;

                        for (tps1 = tps1_min; tps1 <= tps1_max; tps1 += 2)
                        {

                            tps2_min = max(abs(tps1 - chain.dspin), abs(tb3 - chain.dspin));
                            tps2_max = min(tps1 + chain.dspin, tb3 + chain.dspin);

                            for (tps2 = tps2_min; tps2 <= tps2_max; tps2 += 2)
                            {

                                // SW 21j

                                key_21j[0] = ti_8;
                                key_21j[1] = ti_7;
                                key_21j[2] = ti_6;
                                key_21j[3] = ti_5;
                                key_21j[4] = tb1;
                                key_21j[5] = tb2;
                                key_21j[6] = tb3;
                                key_21j[7] = tps1;
                                key_21j[8] = tps2;

                                aSW = Wigner_21j_symbol(key_21j, chain); // CHECK MEMORY ALLOCATION

                                // SW 21j

                                key_21j[0] = ti_9;
                                key_21j[1] = ti_10;
                                key_21j[2] = ti_11;
                                key_21j[3] = ti_12;
                                key_21j[4] = tb5;
                                key_21j[5] = tb4;
                                key_21j[6] = tb3;
                                key_21j[7] = tps1;
                                key_21j[8] = tps2;

                                aSE = Wigner_21j_symbol(key_21j, chain); // CHECK MEMORY ALLOCATION

                                df = DIM(tps1) * DIM(tps2);

                                // phase from reflecting back to stored 21j
                                // ph = real_negpow(
                                //         (chain.dspin + chain.dspin + tb1) + (tb1 + chain.dspin + tb2) + (tb2 + chain.dspin + tb3) +   // SW 21j
                                //         (chain.dspin + chain.dspin + tb1) + (tb1 + chain.dspin + tb2) + (tb2 + chain.dspin + tb3) +   // SE 21j ...
                                //         (chain.dspin + chain.dspin + tps1) + (tps1 + chain.dspin + tps2) + (tps2 + chain.dspin + tb3) //
                                // );

                                // simplified
                                ph = real_negpow(2 * tps1 + 2 * tps2 + 3 * tb3);

                                comp_sum(ph * df * aSW * aSE, aS, cS, yS, tS);

                            } // tps2
                        }     // tps1

                        df = DIM(tb1) * DIM(tb2) * DIM(tb3) * DIM(tb4) * DIM(tb5);

                        // two halves computed, assemble
                        comp_sum(ph * df * aN * aS, ampl, c, y, t);

                    } // tb3
                }     // tb4
            }         // tb2
        }             // tb5
    }                 // tb1

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