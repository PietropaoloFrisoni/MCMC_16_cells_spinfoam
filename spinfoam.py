from ctypes import * 
import os
import threading


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
    
    
def Metropolis_Hastings_multithread_run(draws_path, spin, length, sigma, burnin, verbosity, number_of_threads):
    threads = []
    for thread_id in range(1, number_of_threads+1):
     t = threading.Thread(target=Metropolis_Hastings_run, args=(draws_path, spin, length, sigma, burnin, verbosity, thread_id,))
     threads.append(t)
     print(f'Starting run from thread {thread_id}...')   
     t.start()
    
    for t in threads:   
     t.join()     
     
    print(f'All threads completed the run') 
    threads.clear()
    
    
def Hashing_21j_symbols(hash_tables_path, spin):
    spinfoam_functions = CDLL("./lib/libspinfoam.so")
    spinfoam_functions._Z19Hashing_21j_symbolsPci(hash_tables_path.encode(), int(2 * spin))                                         
