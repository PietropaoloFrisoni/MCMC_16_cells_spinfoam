# Computation of 16 cell spinfoam amplitude

Code for the computation of the 16 cell spinfoam amplitude, boundary observables and correlation functions.

## Dependencies

The library depends on:

1. GNU GSL
2. quadmath (GCC extension)
3. OpenMP
4. _wigxjpf_ and _fastwigxj_ [Johansson et al., 2015]
5. [parallel hashmap](https://github.com/greg7mdp/parallel-hashmap) 
6. [progressbar](https://github.com/gipert/progressbar) (optional) 
7. [atpbar](https://github.com/alphatwirl/atpbar)
8. Python

## Compilation

The C++ source code can be compiled with `make` after that both _wigxjpf_ and _fastwigxj_ have been compiled as well. 

This creates the _spinfoam_ shared library, the object files and the binary files. Type `make DEBUG=1` in order to build the debug version.

The shared library searches for specific folders and files depending on the type of computation (see `Usage` below for details). 

The _spinfoam_ shared library has been tested with GCC version 8.1 or greater.


## Usage

The C++ code can be executed with a shell script (see `scripts`) or via Python (see the notebook *Python_interface*).
