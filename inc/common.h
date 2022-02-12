#pragma once

///////////////////////////////////////////////////////////
// Common includes, macros and defines for the library.
///////////////////////////////////////////////////////////

#include <stddef.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <complex.h>


// Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Dimension of the spin-j representation.
#define DIM(two_j) ((long)(two_j) + 1)

// Divide an INTEGER dspin by 2 (exact integer divison).
#define DIV2(two_j) ((int)((two_j) >> 1))

// Squaring macro.
#define SQ(d) ((d)*(d))

