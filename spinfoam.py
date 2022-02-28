from ctypes import * 
import os
import threading
import warnings
import pandas as pd
import numpy as np
import logging
import re
import glob

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
    
# remove a folder with all its content    
def remove_folder(folder):
    os.system(f"rm -rf {folder}")
    
# takes an intertwiner and returns the corresponding angle eigenvalue        
def from_intertwiner_to_angle(matrix_element, spin):
    return ((matrix_element*(matrix_element + 1) - 2*spin*(spin + 1))/(2*spin*(spin + 1)))  


# takes a draw computes the averaged angles: 
def from_draw_to_angles_average(draws_folder, spin, length, sigma, burnin, chain_id, angles_folder):
    draw_path = f"{draws_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
    path_collection = []
    if (os.path.isfile(draw_path)):
      path_collection.append(draw_path)  
    else:
      warnings.warn("Warning: the draw %s was not found" % (draw_path))   
    
    if (len(path_collection) != 0):
      os.makedirs(f"{angles_folder}/j_{spin}/tmp", exist_ok=True)    
      for draw_path in path_collection:
         df = pd.read_csv(draw_path, low_memory=False)
         multeplicity = df[['draw multeplicity']].to_numpy().astype(int)  
         #total_accept_draws = int(df['total accept. draws'][0])
         total_accept_rate = float(df['total accept. rate'][0].strip('%'))
         total_run_time = float(df['total run time'][0].strip(' s'))
         df = df.drop(columns=['draw multeplicity', 'draw amplitude', 'total accept. draws', 'total accept. rate', 'total run time']) 
         angles_matrix = np.matrix(df.values.transpose()).astype(float)
         angles_matrix = np.vectorize(from_intertwiner_to_angle)(angles_matrix, spin)
         angles_matrix = np.matmul(angles_matrix, multeplicity) 
         angles_matrix = angles_matrix.sum(axis=1)/(length-burnin)  
         angles_matrix = angles_matrix.sum(axis=0)/16
         df = pd.DataFrame({'angle average': [angles_matrix[0,0]],
                            'accept. rate': [total_accept_rate],
                            'run time': [total_run_time]}, 
                            index =[f'chain {chain_id}'])
         tmp_angle_path = f"{angles_folder}/j_{spin}/tmp/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
         df.to_csv(tmp_angle_path, index=True) 
    else:
      warnings.warn("I cannot combine chains since there are no chains available")   
    
     
def funza(draw_path, tmp_angle_path):
    parameters = re.findall(r"[-+]?(?:\d*\.\d+|\d+)", draw_path) 
    spin = float(parameters[0])
    length = int(parameters[1])
    sigma = float(parameters[2])
    burnin = int(parameters[3])
    chain_id = int(parameters[4])
    df = pd.read_csv(draw_path, low_memory=False)
    multeplicity = df[['draw multeplicity']].to_numpy().astype(int)  
    #total_accept_draws = int(df['total accept. draws'][0])
    total_accept_rate = float(df['total accept. rate'][0].strip('%'))
    total_run_time = float(df['total run time'][0].strip(' s'))
    df = df.drop(columns=['draw multeplicity', 'draw amplitude', 'total accept. draws', 'total accept. rate', 'total run time']) 
    angles_matrix = np.matrix(df.values.transpose()).astype(float)
    angles_matrix = np.vectorize(from_intertwiner_to_angle)(angles_matrix, spin)
    angles_matrix = np.matmul(angles_matrix, multeplicity) 
    angles_matrix = angles_matrix.sum(axis=1)/(length-burnin)  
    angles_matrix = angles_matrix.sum(axis=0)/16
    df = pd.DataFrame({'angle average': [angles_matrix[0,0]],
                       'accept. rate': [total_accept_rate],
                       'run time': [total_run_time]}, 
                       index = [f'chain {chain_id}'])
    tmp_angle_path = f"{tmp_angle_path}/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
    df.to_csv(tmp_angle_path, index=True) 
    
    
    
def from_draw_to_angles_average_multithreading(draws_folder, spin, length, sigma, burnin, angles_folder, number_of_threads):
    path_collection = []
    for chain_id in range(1, number_of_threads+1):
       draw_path = f"{draws_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_chain_{chain_id}.csv"
       if (os.path.isfile(draw_path)):
          path_collection.append(draw_path)  
       else:
          warnings.warn("Warning: the draw %s was not found" % (draw_path))   
    
    if (len(path_collection) != 0):   
       tmp_angle_path = f"{angles_folder}/j_{spin}/tmp"
       os.makedirs(tmp_angle_path, exist_ok=False)    
       threads = []
       for draw_path in path_collection:
          t = threading.Thread(target=funza, args=(draw_path, tmp_angle_path,))
          threads.append(t)
          t.start()

       # wait for the threads to complete
       for t in threads:
          t.join()   
                
       combined_anges = pd.DataFrame()
       combined_chains = 0
        
       for file_name in glob.glob(f"{tmp_angle_path}/" + '*.csv'):
          combined_chains += 1 
          x = pd.read_csv(file_name, index_col=0, low_memory=False)
          combined_anges = pd.concat([combined_anges, x], ignore_index=False)
    
       os.makedirs(f"{angles_folder}/j_{spin}/", exist_ok=True)
       angle_path = f"{angles_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_combined_chains_{combined_chains}.csv"  
       combined_anges.to_csv(angle_path, index=True)  
       os.system(f"rm -rf {tmp_angle_path}") 
        
    else:
       warnings.warn("I cannot combine chains since there are no chains available")    
    
def retrieve_combined_angles(angles_folder, spin, length, sigma, burnin, combined_chains):   
    path = f"{angles_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}_combined_chains_{combined_chains}.csv"  
    df = pd.read_csv(path, index_col=0)
    return df
    
         
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
