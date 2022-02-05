#pragma once

///////////////////////////////////////////////////////////
// Common includes, macros and defines for the library.
///////////////////////////////////////////////////////////

#include <stddef.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <complex.h>

#include "config.h"

// Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////
// Root folder. Set at library initialization.
///////////////////////////////////////////////////////////////

extern char* DATA_ROOT;

// Dimension of the spin-j representation.
#define DIM(two_j) ((long)(two_j) + 1)

// Converts a dspin to corresponding spin.
#define SPIN(two_j) ((spin)(two_j) * 0.5)

// Divide an INTEGER dspin by 2 (exact integer divison).
#define DIV2(two_j) ((int)((two_j) >> 1))

///////////////////////////////////////////////////////////////
// Spin labeling and recoupling
///////////////////////////////////////////////////////////////

#define ASSIGN_SPINS(two_js) \
    two_j12 = two_js[0]; \
    two_j13 = two_js[1]; \
    two_j14 = two_js[2]; \
    two_j15 = two_js[3]; \
    two_j23 = two_js[4]; \
    two_j24 = two_js[5]; \
    two_j25 = two_js[6]; \
    two_j34 = two_js[7]; \
    two_j35 = two_js[8]; \
    two_j45 = two_js[9];

#define CHECK_NULL_INTW(intw, err) \
    if (two_##intw##_min_allowed > two_##intw##_max_allowed) { \
        warning("intertwiner " #intw " has empty range"); \
        err = true; \
    } \

#define CHECK_ALLOWED_INTW(intw, err) \
if (two_##intw < two_##intw##_min_allowed || two_##intw > two_##intw##_max_allowed) { \
        warning("intertwiner " #intw " must be in [%d %d]", DIV2(two_##intw##_min_allowed), DIV2(two_##intw##_max_allowed)); \
        err = true; \
    } \

#define CHECK_ALLOWED_INTW_RANGE(intw, err) \
    if (two_##intw##_min < two_##intw##_min_allowed || two_##intw##_max > two_##intw##_max_allowed) { \
        warning("intertwiner " #intw " range must be in [%d %d]", DIV2(two_##intw##_min_allowed), DIV2(two_##intw##_max_allowed)); \
        err = true; \
    } \
