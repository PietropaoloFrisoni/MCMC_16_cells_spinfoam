# Computation of 16 cell spinfoam amplitude

Code for the computation of the 16 cell spinfoam amplitude, boundary observables and correlation functions.

## Dependencies

The library depends on:

1. GNU GSL
2. quadmath (GCC extension)
3. OpenMP
4. [wigxjpf](http://fy.chalmers.se/subatom/wigxjpf/) and [fastwigxj](http://fy.chalmers.se/subatom/fastwigxj/)
5. [parallel hashmap](https://github.com/greg7mdp/parallel-hashmap) 
6. [progressbar](https://github.com/gipert/progressbar) (optional) 
7. [atpbar](https://github.com/alphatwirl/atpbar)
8. Python (≥ 3.6)

## Compilation

The C++ source code can be compiled with `make`. The makefile automatically compiles _wigxjpf_ and _fastwigxj_ as well. 
This creates the _spinfoam_ shared library, the object files and the binary files. Type `make DEBUG=1` in order to build the debug version.
The shared library searches for specific folders and files depending on the type of computation (see `Usage` below for details). 
The _spinfoam_ shared library has been tested with GCC version 8.1 or greater.


## Usage

The C++ code can be executed with a shell script (see `scripts`) or via Python (see the notebook *python_interface*). 
It runs a Markov Chain for each provided thread and stores the states in a folder chosen by the user. 
The code style is hybrid between modern C++ and old C. Unfortunately, I do not have time to update and improve it, as this project as been realized as final assignment of a Scientific Computing course with a restrictive deadline. 

The states stored during the random walk Metropolis-Hastings algorithm are used to compute spinfoam observables and correlations functions in the notebook *Operators* (see `notebooks`). The statistical and data analysis is performed using _pandas_. 
