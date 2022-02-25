# Computation of 16 cell spinfoam amplitude

(...)

## Dependencies

The library depends on:

1. GNU GSL
2. quadmath (GCC extension)
3. _wigxjpf_ and _fastwigxj_ [Johansson et al., 2015]
4. OpenMP
5. _parallelhashmap_  [Git repository](https://github.com/greg7mdp/parallel-hashmap) 
6. _progressbar_ (optional) [Git repository](https://github.com/gipert/progressbar) 
7. Python (optional)
8. Python modules: ctypes, os, threading (optional)

## Compilation

The programs can be compiled typing `make`. This compiles the shared library and the tools programs. There are additional flags that can be provided. For example:

- type `make DEBUG=1` to build the debug version

The shared library searches for folders `wigxjpf`, `fastwigxj` and `parallel_hashmap` under `ext/`. It also searches for folders `fastwig_tables` and `hashed_21j` under `data_folder/`.
The library has been tested with GCC version 8.1 or greater.


## Usage (C++ interface)

(...)

## Usage (Python interface)

See the Jupyter notebook **Python_interface**
