from ctypes import *
import os


def spinfoam_compile():
    os.system("make")


def spinfoam_compile_debug():
    os.system("make DEBUG=1")


def spinfoam_clean():
    os.system("make clean")


def Metropolis_Hastings_run(string, spin, length, sigma, burnin, verbosity):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    spinfoam_functions._Z6MH_runPcS_iidii(string.encode(),
                                          b"./data_folder/hashed_21j",
                                          int(2 * spin), length,
                                          c_double(sigma), burnin, verbosity)
