#!/bin/bash

#PBS -l nodes=1:ppn=16,mem=150gb,walltime=72:00:00
#PBS -m e
#PBS -N Soupour1-10
#PBS -M danaco@jax.org
#PBS -q dev_centos7


cd $PBS_O_WORKDIR

singularity exec /projects/ucar-lab/danaco/Software/Soupour/souporcell.sif souporcell_pipeline.py\
 -i /projects/ucar-lab/danaco/lawlon/PBMC_scRNAseq/BAM_Files/1-RNA_pooled_and_merged_Cellranger_Out_50k/outs/possorted_genome_bam.bam\
 -b /projects/ucar-lab/danaco/PBMC_scRNAseqv/GroundBarcodes/1-HTO.all.cell.barcodes.for.demuxlet.txt\
 -f /projects/ucar-lab/danaco/Software/Soupour/GRCh38_full_analysis_set_plus_decoy_hla.fa\
 -t 8 -o /projects/ucar-lab/danaco/Ground/Soupor/1-RNA -k 4
