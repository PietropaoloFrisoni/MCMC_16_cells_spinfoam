from ctypes import *
import os
import threading
import warnings
import pandas as pd
import numpy as np
import logging
import re
import glob

# ignore everything except the message
def custom_formatwarning(msg, *args, **kwargs):
    return str(msg) + '\n'

warnings.formatwarning = custom_formatwarning

# compile for optimized performance
def spinfoam_compile():
    os.system("make")

# compile in debug mode
def spinfoam_compile_debug():
    os.system("make DEBUG=1")

# clean binary, object and shared library folders
def spinfoam_clean():
    os.system("make clean")

# remove a folder with all its content
def remove_folder(folder):
    os.system(f"rm -rf {folder}")

# starts an independent Markov chain for every assigned thread
def Metropolis_Hastings_parallel_run(data_folder, hash_tables_path, spin,
                                     length, sigma, burnin, verbosity,
                                     number_of_threads):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    print(f'Starting {number_of_threads} independent Markov chains...')
    spinfoam_functions._Z15MH_parallel_runPcS_iidiii(data_folder.encode(),
                                                     hash_tables_path.encode(),
                                                     int(2 * spin), length,
                                                     c_double(sigma), burnin,
                                                     verbosity,
                                                     number_of_threads)
    print(f'Completed! All draws have been stored')

# hashing of 21j Wigner symbols
def Hashing_21j_symbols(hash_tables_path, fastwig_tables_path, spin):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    print(f'Hashing all 21j symbols with j <= {spin}...')
    spinfoam_functions._Z19Hashing_21j_symbolsPcS_i(
        hash_tables_path.encode(), fastwig_tables_path.encode(), int(2 * spin))
    print(f'Completed! The hash table has been stored')
