#!bin/bash

#PBS -l nodes=1:ppn=16 walltime=24:00:00
#PBS -m e
#PBS -N demuxTRY
#PBS -M onur.danaci@jax.org

cd $PBS_O_WORKDIR

module load gcc/4.9.2
module load zlib/1.2.8
module load htslib/1.8


#/home/danaco/Software/popscle/bin/popscle dsc-pileup --sam /home/danaco/data/jurkat_293t_downsampled_n500_full_bam.bam --vcf /home/danaco/data/jurkat_293t_exons_only.vcf.withAF.vcf --out /home/danaco/data/pileup

#/home/danaco/Software/popscle/bin/popscle demuxlet --plp /home/danaco/data/pileup --vcf /home/danaco/data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.snps_shared.coding_exon_v25.reorder.autosomes.maf1.vcf --field PL --out /home/danaco/data/demuxTRY
