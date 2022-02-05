#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "mcmc.h"
#include "fastwigxj.h"
#include "wigxjpf.h"
#include "jsymbols.h"
#include "omp.h"


void dmc_run(Chain &chain) {

     std::cout << "Starting sampling chain with parameters:\n\nspin = " << chain.dspin << "\nlength = " << chain.length 
     << "\nsigma = " << chain.sigma << "\nstore path = " << chain.store_path << std::endl;

     double val3j, val6j;

  val3j = fw3jja(2*  5 , 2*  7 , 2*  5 ,
		 2*(-3), 2*  5);

  printf ("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

  val3j = fw3jja(2* 10 , 2* 15 , 2* 10 ,
		 2*(-3), 2* 12);

  printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val3j = fw3jja6(2*  5 , 2*  7 , 2*  5 ,
		  2*(-3), 2*  5 , 2*(-2));

  printf ("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

  val3j = fw3jja6(2* 10 , 2* 15 , 2* 10 ,
		  2*(-3), 2* 12 , 2*(-9));

  printf ("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val6j = fw6jja(2*  3 , 2*  4 , 2*  2 ,
		 2*  2 , 2*  2 , 2*  3 );

  printf ("6J{ 3   4   2;  2   2   3}:      %#25.15f\n", val6j);

  val6j = fw6jja(2* 10 , 2* 15 , 2* 10 ,
		 2*  7 , 2*  7 , 2*  9 );

  printf ("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);

}