#!/bin/bash

#PBS -l nodes=1:ppn=16,walltime=72:00:00
#PBS -m e
#PBS -M danaco@jax.org
#PBS -t 1-10

cd $PBS_O_WORKDIR

cd /projects/ucar-lab/danaco/PBMC_scRNAseqv

module load gcc/4.9.2
module load zlib/1.2.8
module load htslib/1.8
module load R

R --no-save -f Demux.par.R
