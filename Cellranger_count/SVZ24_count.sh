#! /bin/bash
#SBATCH -A LSENS2018-3-3
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH -J Run_SVZ24_count
#SBATCH -o SVZ24count%j.out
#SBATCH -e SVZ24count%j.err
#SBATCH -p dell
module --force purge
module load SoftwareTree/Haswell
module load cellranger/6.0.0
cellranger count --id=SVZ24_count \
--fastqs=/projects/fs1/mtorbati/Torbati/CellRanger_V7/mkfastq/firstRun/mkfastq_1stRun/outs/fastq_path/Chromium_20180912,/projects/fs1/mtorbati/Torbati/CellRanger_V7/mkfastq/secondRun/mkfastq_2ndRun/outs/fastq_path/Chromium_20180912 \
--sample=SVZ24 \
--transcriptome=/projects/fs1/mtorbati/Torbati/cellranger_custom_RefGen/custom_ref_gene/Mus.musculus_genome_EGFP
exit 0
