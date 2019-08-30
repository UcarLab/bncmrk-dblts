#!/bin/bash

#PBS -l nodes=1:ppn=16,walltime=72:00:00
#PBS -m e
#PBS -N hto.smm
#PBS -M onur.danaci@jax.org 

cd $PBS_O_WORKDIR
cd /projects/ucar-lab/danaco/Ground

module unload R
#module load Anaconda/4.2.0
#source activate scSplit
module unload Anaconda/4.2.0
module load R/3.5.1
module load gcc/4.9.2

R --no-save -f HTO.GMM.analysis.R
