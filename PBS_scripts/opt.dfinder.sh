#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=150gb,walltime=24:00:00
#PBS -m e
#PBS -N dfinder_HTO1.2_opt_expected_10per
#PBS -M danaco@jax.org
#(pbs) #-q #high_mem

cd $PBS_O_WORKDIR

module load Anaconda/4.2.0
source activate ML

module load R/3.5.1
module load gcc/4.9.2

R --no-save -f Optimized.DF.analysis.R
