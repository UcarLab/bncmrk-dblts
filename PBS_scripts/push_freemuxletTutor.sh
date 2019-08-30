#!/bin/bash
#PBS -l nodes=1:ppn=16 walltime=24:00:00
#PBS -m e
#PBS -N freemuxTutor
#PBS -M onur.danaci@jax.org

cd $PBS_O_WORKDIR

module load gcc/4.9.2
module load zlib/1.2.8
module load htslib/1.8

#/home/danaco/Software/popscle/bin/popscle dsc-pileup --sam /home/danaco/data/jurkat_293t_downsampled_n500_full_bam.bam --vcf /home/danaco/data/jurkat_293t_exons_only.vcf.withAF.vcf.gz --out /home/danaco/result/dsc_pileup.jurkat_293t.pooled
/home/danaco/Software/popscle/bin/popscle freemuxlet --plp /home/danaco/result/dsc_pileup.jurkat_293t.pooled --nsample 2 --out /home/danaco/result/jurkat_293t_freemuxlet.pooled
