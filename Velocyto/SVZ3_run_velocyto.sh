#!/bin/sh
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -t 20:00:00
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -J velocyto
#SBATCH -o SVZ3_velocyto.%j.out
#SBATCH -e SVZ3_velocyto.%j.err
module purge
module load GCC/10.2.0
module load velocyto/0.17.17
module load SAMtools/1.12
repeats="/projects/fs1/mtorbati/Torbati/CellRanger_V6/runVelocyto/mm10_rmsk.gtf"
transcriptome="/projects/fs1/mtorbati/ref_gene_mouse/EGFP.gtf"
cellranger_output="/projects/fs1/mtorbati/Torbati/CellRanger_V6/cellranger_count/SVZ3/SVZ3_count"
velocyto run10x -m $repeats \
                $cellranger_output \
                $transcriptome
exit 0
