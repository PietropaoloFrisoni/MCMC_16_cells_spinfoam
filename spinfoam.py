from ctypes import * 
import os
import threading


def spinfoam_compile():
    os.system("make")


    
def spinfoam_compile_debug():
    os.system("make DEBUG=1")


    
def spinfoam_clean():
    os.system("make clean")
    

    
def Metropolis_Hastings_parallel_run(draws_path, hash_tables_path, spin, length, sigma, burnin, verbosity, number_of_threads):
    spinfoam_functions = CDLL("./lib/libspinfoam.so") 
    print(f'Starting {number_of_threads} independent Markov chains...')
    spinfoam_functions._Z15MH_parallel_runPcS_iidiii(draws_path.encode(),
                                          hash_tables_path.encode(),
                                          int(2 * spin), length,
                                          c_double(sigma), burnin, verbosity, number_of_threads)     
    print(f'Completed! All draws have been stored')

    
    
def Hashing_21j_symbols(hash_tables_path, fastwig_tables_path, spin):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    print(f'Hashing all 21j symbols with j <= {spin}...')
    spinfoam_functions._Z19Hashing_21j_symbolsPcS_i(hash_tables_path.encode(), fastwig_tables_path.encode(), int(2 * spin))           
    print(f'Completed! The hash table has been stored')    
