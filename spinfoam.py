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


def spinfoam_compile():
    os.system("make")


    
def spinfoam_compile_debug():
    os.system("make DEBUG=1")


    
def spinfoam_clean():
    os.system("make clean")
    
    
def old_find_draws(draws_folder, spin, length, sigma, burnin, number_of_chains_to_process):
    for chain_id in range(1,number_of_chains_to_process+1):
      draw = f"{draws_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
      if (os.path.isfile(draw)):
        print(f"esiste per chain_id {chain_id}")
      else:
        warnings.warn("Warning: chain id %d not found" % (chain_id))  
        
def path_to_draw(draws_folder, spin, length, sigma, burnin, chain_id):
      draw_path = f"{draws_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
      if (os.path.isfile(draw_path)):
        return draw_path
      else:
        return False        

        
def from_intertwiner_to_angle(matrix_element,spin):
    return ((matrix_element*(matrix_element + 1) - 2*spin*(spin + 1))/(2*spin*(spin + 1)))        
        
# takes an existing draw path with some parameters and returns the average of all the boundary angles operator 
def from_draw_to_angles_average(draw, spin, length, burnin):
    df = pd.read_csv(draw, low_memory=False)
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
    return angles_matrix[0,0]

    
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
