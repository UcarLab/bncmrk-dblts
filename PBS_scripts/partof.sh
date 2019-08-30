#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=150gb,walltime=72:00:00
#PBS -m e
#PBS -N seurat_with_hto_doubs
#PBS -M onur.danaci@jax.org 

cd $PBS_O_WORKDIR
cd /projects/ucar-lab/danaco/Ground

module unload R
#module load Anaconda/4.2.0
#source activate scSplit
module load Anaconda/4.2.0
source activate ML

module load R/3.5.1
module load gcc/4.9.2

R --no-save -f Partof.R
