from ctypes import *


def load_spinfoam_library(path_to_spinfoam_library):
    spinfoam_functions = CDLL(path_to_spinfoam_library)


def Metropolis_Hastings_run(spin, length, sigma, burnin, verbosity):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    spinfoam_functions._Z6MH_runPcS_iidii(b"./data_folder/collected_draws", b"./data_folder/hashed_21j", spin, length, c_double(sigma), burnin, verbosity)


