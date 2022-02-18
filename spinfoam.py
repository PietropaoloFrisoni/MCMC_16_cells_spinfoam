from ctypes import *
import os


def spinfoam_compile():
    os.system("make")


def spinfoam_compile_debug():
    os.system("make DEBUG=1")


def spinfoam_clean():
    os.system("make clean")


def Metropolis_Hastings_run(draws_path, spin, length, sigma, burnin, verbosity, thread_id):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    spinfoam_functions._Z6MH_runPcS_iidiii(draws_path.encode(),
                                          b"./data_folder/hashed_21j",
                                          int(2 * spin), length,
                                          c_double(sigma), burnin, verbosity, thread_id)


def Hashing_21j_symbols(hash_tables_path, spin):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    spinfoam_functions._Z19Hashing_21j_symbolsPci(hash_tables_path.encode(), int(2 * spin))                                         
