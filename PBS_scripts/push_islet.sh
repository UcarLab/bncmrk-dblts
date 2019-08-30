#!/bin/bash

#PBS -l nodes=1:ppn=16,walltime=48:00:00
#PBS -m e
#PBS -N freemux_islet
#PBS -M onur.danaci@jax.org

cd $PBS_O_WORKDIR

module load gcc/4.9.2
module load zlib/1.2.8
module load htslib/1.8

# generate pileups

/home/danaco/Software/popscle/bin/popscle dsc-pileup --sam /projects/ucar-lab/danaco/Multiplexed/possorted_genome_bam.bam --vcf /home/danaco/data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.snps_shared.coding_exon_v25.reorder.autosomes.maf1.vcf \
--group-list  /projects/ucar-lab/danaco/Multiplexed/MS19006_MS19007.barcodes.txt \
 --out /home/danaco/result/islet/Freemuxlet.Output

#/projects/stitzel-lab/lawlon/Software/popscle/bin/popscle dsc-pileup --sam possorted_genome_bam.bam --vcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.snps_shared.coding_exon_v25.reorder.autosomes.maf1.vcf \
#--group-list MS19006_MS19007.barcodes.txt --out Freemuxlet.Output

# run freemux, to demultiplex identity of pooled data (consisting of 2 samples)
printf "Pileup Complete\n"

/home/danaco/Software/popscle/bin/popscle freemuxlet --plp /home/danaco/result/islet/Freemuxlet.Output \
 --group-list /projects/ucar-lab/danaco/Multiplexed/MS19006_MS19007.barcodes.txt \
 --out  /home/danaco/result/islet/Freemuxlet.Output --nsample 2 --aux-files

