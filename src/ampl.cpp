#include <iostream>

#include "utilities.h"
#include "ampl.h"

double pce_amplitude_c16(Chain chain){

    // boundary data
    // spins are to be read counterclockwise
    // starting from top
    int ti_1 = chain.indices[0];
    int ti_2 = chain.indices[1];
    int ti_3 = chain.indices[2];
    int ti_4 = chain.indices[3];

    int ti_5 = chain.indices[4];
    int ti_6 = chain.indices[5];
    int ti_7 = chain.indices[6];
    int ti_8 = chain.indices[7];

    int ti_9 = chain.indices[8];
    int ti_10 = chain.indices[9];
    int ti_11 = chain.indices[10];
    int ti_12 = chain.indices[11];

    int ti_13 = chain.indices[12];
    int ti_14 = chain.indices[13];
    int ti_15 = chain.indices[14];
    int ti_16 = chain.indices[15];

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
    
}