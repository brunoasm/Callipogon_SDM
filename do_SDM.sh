#!/bin/bash
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2-00:00:00 # Runtime in minutes
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH -o SDM_Callipogon.out # stdout and stderr go to this file

module load R/3.2.2-fasrc03 geos/3.4.2-fasrc01 jdk/1.8.0_45-fasrc01 gdal/1.11.1-fasrc01 proj/4.9.2-fasrc01 

export R_LIBS_USER=$HOME/R:$R_LIBS_USER #for some reason, when running from slurm script it can't find my library

Rscript --no-restore SDM_Callipogon.R 
