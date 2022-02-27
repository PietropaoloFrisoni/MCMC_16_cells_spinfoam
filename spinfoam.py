from ctypes import * 
import os
import threading
import warnings
import pandas as pd
import numpy as np

def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
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
    
# takes an intertwiner and returns the corresponding angle eigenvalue        
def from_intertwiner_to_angle(matrix_element, spin):
    return ((matrix_element*(matrix_element + 1) - 2*spin*(spin + 1))/(2*spin*(spin + 1)))        
        
# takes a draw and returns the following array: 
# [average of all the boundary angles operator, total_accept_rate, total_run_time]
def from_draw_to_angles_average(draws_folder, spin, length, sigma, burnin, chain_id):
    draw_path = f"{draws_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
    if (os.path.isfile(draw_path)):
      df = pd.read_csv(draw_path, low_memory=False)
      multeplicity = df[['draw multeplicity']].to_numpy().astype(int)  
      total_accept_draws = int(df['total accept. draws'][0])
      total_accept_rate = float(df['total accept. rate'][0].strip('%'))
      total_run_time = float(df['total run time'][0].strip(' s'))
      df = df.drop(columns=['draw multeplicity', 'draw amplitude', 'total accept. draws', 'total accept. rate', 'total run time']) 
      angles_matrix = np.matrix(df.values.transpose()).astype(float)
      angles_matrix = np.vectorize(from_intertwiner_to_angle)(angles_matrix, spin)
      angles_matrix = np.matmul(angles_matrix, multeplicity) 
      angles_matrix = angles_matrix.sum(axis=1)/(length-burnin)  
      angles_matrix = angles_matrix.sum(axis=0)/16
      return [angles_matrix[0,0], total_accept_rate, total_run_time]
    else:
      warnings.warn("Warning: the draw %s was not found" % (draw_path))    

# starts an independent Markov chain for every thread assigned 
def Metropolis_Hastings_parallel_run(draws_path, hash_tables_path, spin, length, sigma, burnin, verbosity, number_of_threads):
    spinfoam_functions = CDLL("./lib/libspinfoam.so") 
    print(f'Starting {number_of_threads} independent Markov chains...')
    spinfoam_functions._Z15MH_parallel_runPcS_iidiii(draws_path.encode(),
                                          hash_tables_path.encode(),
                                          int(2 * spin), length,
                                          c_double(sigma), burnin, verbosity, number_of_threads)     
    print(f'Completed! All draws have been stored')

# hashing of 21j Wigner symbols     
def Hashing_21j_symbols(hash_tables_path, fastwig_tables_path, spin):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    print(f'Hashing all 21j symbols with j <= {spin}...')
    spinfoam_functions._Z19Hashing_21j_symbolsPcS_i(hash_tables_path.encode(), fastwig_tables_path.encode(), int(2 * spin))           
    print(f'Completed! The hash table has been stored')    
