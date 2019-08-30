#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=100gb,walltime=72:00:00
#PBS -m e
#PBS -N scrublet_out
#PBS -M onur.danaci@jax.org 

cd $PBS_O_WORKDIR
cd /projects/ucar-lab/danaco/bncmrk-dblts/Python

module unload R
#module load Anaconda/4.2.0
#source activate scSplit
module load Anaconda/4.2.0
source activate scSplit

#module load R/3.5.1
module load gcc/4.9.2

python scr_pipe_script.py
