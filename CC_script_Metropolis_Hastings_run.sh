#!/bin/bash
#SBATCH -A def-vidotto
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=3-0:00:00
#SBATCH --job-name=MH_parallel
#SBATCH --output=MH_parallel.log
#SBATCH --error=MH_parallel.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=pfrisoni@uwo.ca


make clean

module load gsl

make


echo "Running on: $SLURM_NODELIST"
echo

export LD_LIBRARY_PATH="./lib":$LD_LIBRARY_PATH

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

FASTWIG_TABLES_PATH=./data_folder/fastwig_tables/

HASH_TABLES_STORE_PATH=./data_folder/hashed_21j/

DRAWS_STORE_PATH=./data_folder/collected_draws

DSPIN=4

LENGTH=1000000

SIGMA=0.35

BURNIN=10

VERBOSITY=0

now=$(date)
echo
echo "Starting Metropolis-Hastings algorithm... [ NUM_OF_THREADS = ${OMP_NUM_THREADS} DSPIN = ${DSPIN} LENGTH = ${LENGTH} ]... (now: $now)"

bin/Metropolis_Hastings_parallel_run  $OMP_NUM_THREADS  $DSPIN  $LENGTH  $SIGMA  $BURNIN  $VERBOSITY  $DRAWS_STORE_PATH  $HASH_TABLES_STORE_PATH
            
now=$(date)
echo
echo "... completed! (now: $now)"

