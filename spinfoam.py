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


# takes an intertwiner and returns the corresponding angle eigenvalue
def from_intertwiner_to_angle(matrix_element, spin):
    return ((matrix_element * (matrix_element + 1) - 2 * spin * (spin + 1)) /
            (2 * spin * (spin + 1)))


def from_draws_to_angles(folder_prefix, spin, length, burnin, angle_path, chain_id): 

    # load in memory the stored draws
    draw_path = f"{folder_prefix}/draws/draws_chain_{chain_id}.csv"
    df = pd.read_csv(draw_path, low_memory=False) 
  
    # retrieving relevant parameters
    multeplicity = df[['draw multeplicity']].to_numpy().astype(int)
    total_accept_draws = int(df['total accept. draws'][0])
    total_accept_rate = float(df['total accept. rate'][0].strip('%'))
    total_run_time = float(df['total run time'][0].strip(' s'))
    
    # dropping columns
    df = df.drop(columns=[
        'draw multeplicity', 'draw amplitude', 'total accept. draws',
        'total accept. rate', 'total run time'
    ])
    
    # from csv to matrix in order to use the numpy optimized routines
    angles_matrix = np.matrix(df.values.transpose()).astype(float)
    
    # from intertwiners to angles 
    angles_matrix = np.vectorize(from_intertwiner_to_angle)(angles_matrix, spin)    
    angles_sq_matrix = np.power(angles_matrix, 2)
    
    # from intertwiners to angles sq 
    angles_matrix = np.matmul(angles_matrix, multeplicity)
    angles_sq_matrix = np.matmul(angles_sq_matrix, multeplicity)
    
    # average of angles
    angles_matrix = angles_matrix.sum(axis=1) / (length - burnin)
    angles_sq_matrix = angles_sq_matrix.sum(axis=1) / (length - burnin)
    
    nodes_angles_avg_matrix = angles_matrix.sum(axis=0) / 16
    nodes_angles_sq_avg_matrix = angles_sq_matrix.sum(axis=0) / 16
    
    df = pd.DataFrame(
        {
            'node 1': [angles_matrix[0,0], angles_sq_matrix[0,0]],
            'node 2': [angles_matrix[1,0], angles_sq_matrix[1,0]],
            'node 3': [angles_matrix[2,0], angles_sq_matrix[2,0]],
            'node 4': [angles_matrix[3,0], angles_sq_matrix[3,0]],
            'node 5': [angles_matrix[4,0], angles_sq_matrix[4,0]],
            'node 6': [angles_matrix[5,0], angles_sq_matrix[5,0]],
            'node 7': [angles_matrix[6,0], angles_sq_matrix[6,0]],
            'node 8': [angles_matrix[7,0], angles_sq_matrix[7,0]],
            'node 9': [angles_matrix[8,0], angles_sq_matrix[8,0]],
            'node 10': [angles_matrix[9,0], angles_sq_matrix[9,0]],
            'node 11': [angles_matrix[10,0], angles_sq_matrix[10,0]],
            'node 12': [angles_matrix[11,0], angles_sq_matrix[11,0]],
            'node 13': [angles_matrix[12,0], angles_sq_matrix[12,0]],
            'node 14': [angles_matrix[13,0], angles_sq_matrix[13,0]],
            'node 15': [angles_matrix[14,0], angles_sq_matrix[14,0]],
            'node 16': [angles_matrix[15,0], angles_sq_matrix[15,0]],
            'nodes avg': [nodes_angles_avg_matrix[0,0], nodes_angles_sq_avg_matrix[0,0]]            
        },
        index=[f'angle average', f'angle sq average'])
        
    angle_path_chain = f"{angle_path}/angles_chain_{chain_id}.csv"
    df.to_csv(angle_path_chain, index=True)
    
    
def from_draws_to_angles_correlations(folder_prefix, spin, length, burnin, angle_correlations_path, chain_id=1): 

    # load in memory the stored draws
    draw_path = f"{folder_prefix}/draws/draws_chain_{chain_id}.csv"
    df = pd.read_csv(draw_path, low_memory=False) 
  
    # retrieving relevant parameters
    multeplicity = df[['draw multeplicity']].to_numpy().astype(int)
    total_accept_draws = int(df['total accept. draws'][0])
    total_accept_rate = float(df['total accept. rate'][0].strip('%'))
    total_run_time = float(df['total run time'][0].strip(' s'))
    
    # dropping columns
    df = df.drop(columns=[
        'draw multeplicity', 'draw amplitude', 'total accept. draws',
        'total accept. rate', 'total run time'
    ])
    
    # from csv to matrix in order to use the numpy optimized routines
    angles_matrix = np.matrix(df.values).astype(float)
    
    # from intertwiners to angles 
    angles_matrix = np.vectorize(from_intertwiner_to_angle)(angles_matrix, spin) 
 
    indices_collection = []
    values_collection = []    
        
    for node_1 in range(16):    
        for node_2 in range(node_1, 16):
            corr = 0.0
            for n in range(total_accept_draws):
                corr += angles_matrix[n, node_1]*angles_matrix[n, node_2]*multeplicity[n]
            corr /= (length - burnin)    
            indices_collection.append(f"<O({node_1+1}),O({node_2+1})>")
            values_collection.append(corr[0])
    
    df = pd.DataFrame(
        {                
             f'j={spin}': values_collection[:],  
        },
        index = indices_collection)
            
    angle_correlations_path_chain = f"{angle_correlations_path}/angles_correlations_chain_{chain_id}.csv"
    df = df.T
    df.to_csv(angle_correlations_path_chain, index=True)       
      
    
def angles_correlations_compute(data_folder, spin, length,
                                sigma, burnin, number_of_threads):   
                   
    folder_prefix = f"{data_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}"
    chain_id_collection = []
    
    for chain_id in range(1, number_of_threads + 1):
        draw_path = f"{folder_prefix}/draws/draws_chain_{chain_id}.csv"
        if (os.path.isfile(draw_path)):
            chain_id_collection.append(chain_id)
        else:
            warnings.warn("Warning: the draw %s was not found" % (draw_path))
            
    chains_to_assemble = len(chain_id_collection)        

    if (chains_to_assemble != 0):
    
        angle_correlations_path = f"{folder_prefix}/operators/angles_correlations"
        os.makedirs(angle_correlations_path, exist_ok=True)
        
        print(f'Converting {chains_to_assemble} chains from draws to angles correlations...')
        
        threads = []
        for chain_id in chain_id_collection:
            t = threading.Thread(target=from_draws_to_angles_correlations,
                                 args=(folder_prefix, 
                                       spin, length, burnin,
                                       angle_correlations_path, 
                                       chain_id,
                                      ))
            threads.append(t)
            t.start()

        # wait for the threads to complete
        for t in threads:
            t.join()    
            
        print(f'Completed! All draws have been processed')     
        
        
        

def angles_compute(data_folder, spin, length,
                   sigma, burnin, number_of_threads):   
                   
    folder_prefix = f"{data_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}"
    chain_id_collection = []
    
    for chain_id in range(1, number_of_threads + 1):
        draw_path = f"{folder_prefix}/draws/draws_chain_{chain_id}.csv"
        if (os.path.isfile(draw_path)):
            chain_id_collection.append(chain_id)
        else:
            warnings.warn("Warning: the draw %s was not found" % (draw_path))
            
    chains_to_assemble = len(chain_id_collection)        

    if (chains_to_assemble != 0):
    
        angle_path = f"{folder_prefix}/operators/angles"
        os.makedirs(angle_path, exist_ok=True)
        
        print(f'Converting {chains_to_assemble} chains from draws to averaged angles...')
        
        threads = []
        for chain_id in chain_id_collection:
            t = threading.Thread(target=from_draws_to_angles,
                                 args=(folder_prefix, 
                                       spin, length, burnin,
                                       angle_path, 
                                       chain_id,
                                      ))
            threads.append(t)
            t.start()

        # wait for the threads to complete
        for t in threads:
            t.join()
            
        print(f'Completed! All draws have been processed')    
        
        print(f'Assembling {chains_to_assemble} chains...')
        
        matrix = np.zeros((chains_to_assemble, 3))
        
        for chain_id in range(chains_to_assemble):
        
            angle_path_chain = f"{angle_path}/angles_chain_{chain_id+1}.csv"
            df = pd.read_csv(angle_path_chain, index_col=0, low_memory=False)
            matrix[chain_id, 0] = df["nodes avg"]['angle average']
            matrix[chain_id, 1] = df["nodes avg"]['angle sq average']
            matrix[chain_id, 2] = np.sqrt(matrix[chain_id, 1] - matrix[chain_id, 0]*matrix[chain_id, 0])
          
        values = np.mean(matrix, axis=0)
            
        angles_fluc = np.std(matrix, axis=0)[0]     
        angles_computed = values[0]
        angles_sq_computed = values[1]
        angles_quantum_spread_computed = values[2]
        
        df = pd.DataFrame(
        {
            'angles_average': [angles_computed],
            'angles_std': [angles_fluc],
            'angles_quantum_spread': [angles_quantum_spread_computed],
            'chains_length' : [length],
            'chains_assembled' : [chains_to_assemble],
            'sigma' : [sigma],
            'burnin' : [burnin],
            
        },
        index=[f'j={spin}'])
        
        assembled_angle_path = f"{data_folder}/final_tables"
        os.makedirs(assembled_angle_path, exist_ok=True)
        
        assembled_angle_path = f"{assembled_angle_path}/angles_j={spin}_chains_combined_{chains_to_assemble}.csv"
        df.to_csv(assembled_angle_path, index=True)
        
        print(f'Done')

    else:
        warnings.warn("I can't compute angles since there are no chains available")


def retrieve_combined_angles(data_folder, spin, chains_combined):
    path = f"{data_folder}/final_tables/angles_j={spin}_chains_combined_{chains_combined}.csv"
    df = pd.read_csv(path, index_col=0)
    return df


def retrieve_correlations(data_folder, spin, length,
                                sigma, burnin, chain_id=1):
    
    folder_prefix = f"{data_folder}/j_{spin}/N_{length}__sigma_{sigma}__burnin_{burnin}"
    angle_correlations_path = f"{folder_prefix}/operators/angles_correlations"
    angle_correlations_path_chain = f"{angle_correlations_path}/angles_correlations_chain_{chain_id}.csv"
    
    df = pd.read_csv(angle_correlations_path_chain, index_col=0)
    return df



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
