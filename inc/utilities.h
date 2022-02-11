#pragma once

#include <iostream>

#include "error.h"

#define S3J(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3) \
	fw3jja6(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3)
#define S6J(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6) \
	fw6jja(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6)

////////////////////////////////////////////////////////////////////////
// Simple numerical utilities.
////////////////////////////////////////////////////////////////////////

// Checks for triangular inequalities.
#define check_triangle_return0(tj1, tj2, tj3)           \
	if (!(abs(tj1 - tj2) <= tj3 && tj3 <= (tj1 + tj2))) \
	{                                                   \
		return 0;                                       \
	}

#define check_triangle_continue(tj1, tj2, tj3)          \
	if (!(abs(tj1 - tj2) <= tj3 && tj3 <= (tj1 + tj2))) \
	{                                                   \
		continue;                                       \
	}

#define triangle(tj1, tj2, tj3) \
	(abs(tj1 - tj2) <= tj3 && tj3 <= (tj1 + tj2))

// Checks for integrity.
#define check_integer_2sum_return0(tl1, tl2) \
	{                                        \
		if ((tl1 + tl2) % 2 != 0)            \
		{                                    \
			return 0;                        \
		}                                    \
	}

#define check_integer_3sum_return0(tl1, tl2, tl3) \
	{                                             \
		if ((tl1 + tl2 + tl3) % 2 != 0)           \
		{                                         \
			return 0;                             \
		}                                         \
	}

#define check_integer_2sum_continue(tl1, tl2) \
	{                                         \
		if ((tl1 + tl2) % 2 != 0)             \
		{                                     \
			continue;                         \
		}                                     \
	}

#define check_integer_3sum_continue(tl1, tl2, tl3) \
	{                                              \
		if ((tl1 + tl2 + tl3) % 2 != 0)            \
		{                                          \
			continue;                              \
		}                                          \
	}

// Checks that the dspin value corresponds to an integer spin.
#define ensure_integer_spin(tj)            \
	{                                      \
		if ((tj) % 2 != 0)                 \
		{                                  \
			error("integer check failed"); \
		}                                  \
	}

// Compensated summation.
// TODO: some testing reveals that using a LONG DOUBLE (80-bit) var
//       to accumulate the sum has similar precision but it's
//       ~3 times faster. switch to LONG DOUBLE then
#define comp_sum(inp, sum, c, y, t) \
	{                               \
		y = inp - c;                \
		t = sum + y;                \
		c = (t - sum) - y;          \
		sum = t;                    \
	}

// Returns true if spin is integer.
#define is_integer(tj) ((tj % 2) == 0)

// Returns true if spin is semi-integer.
#define is_semi_integer(tj) ((tj % 2) == 1)

// Computes (-1)^j for integer j.
// For semi-integer j it returns the real part.
static inline int real_negpow(int tj)
{

	ensure_integer_spin(tj);

	int j = tj / 2;

	if (j % 2 == 0)
	{
		return 1;
	}
	return -1;
}

// Returns the maximum of two integers.
static inline int max(int n1, int n2)
{

	if (n1 > n2)
	{
		return n1;
	}
	return n2;
}

// Returns the minimum of two integers.
static inline int min(int n1, int n2)
{

	if (n1 < n2)
	{
		return n1;
	}
	return n2;
}

// display a matrix for given rows and columns

static inline void print_matrix(int **matrix, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			std::cout << matrix[i][j] << " ";
		}
		
		std::cout << std::endl;
	}
}
