#!/bin/bash
#SBATCH -A def-vidotto
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3-0:00:00
#SBATCH --job-name=hash_21j
#SBATCH --output=hash_21j.log
#SBATCH --error=hash_21j.err
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=pfrisoni@uwo.ca


echo "Running on: $SLURM_NODELIST"
echo

# start commands

export LD_LIBRARY_PATH="./lib":$LD_LIBRARY_PATH

FASTWIG_TABLES_PATH=./data_folder/fastwig_tables/

HASH_TABLES_STORE_PATH=./data_folder/hashed_21j/

DSPIN=6

VERBOSITY=0

bin/Hashing_21j    $DSPIN    $HASH_TABLES_STORE_PATH    $FASTWIG_TABLES_PATH    $VERBOSITY
  
echo
echo "completed."
