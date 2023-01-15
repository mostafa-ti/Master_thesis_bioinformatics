# Investigation of age-related changes in neuroblast populations using RNA velocity


# Project overview
README file for a master project in Bioinformatics
* BINP52 - Master project - 60 ECTS
* Master programme in Bioinformatics
* Department of Biology, Faculty of Science, Lund University, Sweden

This work aimed to understand the aging signature on the neurogenesis process in the brain. For this purpose, three regions of mice brain, the Dentate Gyrus (DG), Subventricular zone (SVZ), and Olfactory Bulb (OB), in three age sample groups young (3 months), adults (14 months), and aged (24 months) were compared. This project provided an RNA velocity map of neurogenic niches of mice, together with substantial new knowledge about the biology of newly born neuroblasts and how they are affected by aging. Overall, our data advanced our understanding of the aging signature on the neurogenesis process in the brain.

# People

* Student: Mostafa Torbati
* Supervisors: Henrik Ahlenius (henrik.ahlenius@med.lu.se) and Jonas Fritze (jonas.fritze@med.lu.se)
* Institute: Department of Clinical Sciences Lund, Division of Neurology, Lund Stem Cell Center, Lund University, Sweden 
* Corresponding author: Mostafa Torbati

# Hardware specifications
The following hardware was used:
* Apple MacBook Pro (Mid 2014, macOS Big Sur Version 11.7) used as a local computer.
* [LUNARC's Aurora service HPC Desktop](https://lunarc-documentation.readthedocs.io/en/latest/using_hpc_desktop/ "LUNARC Documentation") (High-performance computing), is used for preprocessing and visualizing the raw data.

# Software versions
* Cell Ranger v6.0
* Scanpy v1.7.2
* Velocyto v0.17.17
* ScVelo v0.2.4

# Installation
On this project, we're working with large data and software such as Cell Ranger, which need High-performance computing (HPC) power. All the bioinformatics workflows for this project have been done on LUNARC Aurora service HPC Desktop ( the center for scientific and technical computing at Lund University).
All the software is pre-installed on the LUNARC clusters so that users can access and load software packages based on the project's requirements.
To efficiently management of the compute resources, we need to follow the [SLURM](https://lunarc-documentation.readthedocs.io/en/latest/batch_system/ "SLURM - the batch system at LUNARC")
 job scheduler, which in simple words, is a bash script that loads all the necessary packages and sends it to the system backend.
Here is an example of a `bash script` for submitting a job on Aurora:
```bash
#! /bin/bash
#SBATCH -A LSENS2018-3-3 # the ID of our Aurora project
#SBATCH -n 20 # how many processor cores to use
#SBATCH -N 1 # how many processors to use
#SBATCH -t 24:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J [JOB_NAME]# name of the job
#SBATCH -o [OUTPUT_report]%j.out # stdout log file
#SBATCH -e [ERROR_report]%j.err # stderr log file
#SBATCH -p dell # which partition to use
# load the necessary modules
module purge #remove all currently loaded modules 
module load [packages/software] #load packages
exit 0
```
The script is then submitted to the scheduler using the command:
```bash
sbatch my_script.sh
```

# Flowchart

# Cell Ranger workflow
## Build a Custom Reference
For this project, we used transgenic mice in which Enhanced Green Fluorescent Protein (EGFP) is expressed under the control of the DCX promoter. As EGFP is not included in the pre-built reference genome, we first created a custom reference genome and added EGFP to an existing mouse reference.

The cellranger mkref command is used to build a custom reference for use with the Cell Ranger pipeline.

### Prepare input files
* The following files are required to build the reference:
* FASTA file containing the genome sequences
* GTF file containing the gene annotation

**get the FASTA file**
```bash
wget ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus//dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
```
**get the GTF file**
```bash
wget ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus//Mus_musculus.GRCm38.102.chr.gtf.gz
```
#### filter the GTF file
GTF files can contain entries for non-polyA transcripts that overlap with protein-coding gene models.To remove these entries from the GTF, we need to use `cellranger mkgtf` nad add the filter argument `--attribute=gene_biotype:protein_coding` to the mkgtf command:
```bash
#! /bin/bash
#SBATCH -A LSENS2018-3-3 # the ID of our Aurora project
#SBATCH -n 20 # how many processor cores to use
#SBATCH -N 1 # how many processors to use
#SBATCH -t 24:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J filter_gtf_Mus # name of the job
#SBATCH -o filter_gtf_Mus%j.out # stdout log file
#SBATCH -e filter_gtf_Mus%j.err # stderr log file
#SBATCH -p dell # which partition to use
module purge
module load cellranger/6.0.0 \
cellranger mkgtf [PATH_TO_GTFfile]\
Mus_musculus.filtered_gtf.gtf \ #This will output the file Mus_musculus.filtered_gtf.gtf, which will be used in the cellranger mkref step.
--attribute=gene_biotype:protein_coding
exit 0
```

### Setup the command for cellranger mkref
```bash
#! /bin/bash
#SBATCH -A LSENS2018-3-3 # the ID of our Aurora project
#SBATCH -n 20 # how many processor cores to use
#SBATCH -N 1 # how many processors to use
#SBATCH -t 24:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J mkref_Mus_musculus # name of the job
#SBATCH -o mkref_Mus_musculus%j.out # stdout log file
#SBATCH -e mkref_Mus_musculus%j.err # stderr log file
#SBATCH -p dell # which partition to use
module purge
module load cellranger/6.0.0
cellranger mkref --genome=Mus.musculus_genome \ #name of the custom reference folder
--fasta=[replace with the path of FASTA file containing the genome sequences] \
--genes=[replace with the path  of _Mus_musculus.filtered_gtf.gtf_ created in the previous step]
exit 0
```

### Add a marker gene (EGFP) to the FASTA and GTF
 **get the EGFP gene sequence**
 [EGFP complete sequence](https://www.ncbi.nlm.nih.gov/nuccore/U55762.1?report=fasta "EGFP complete sequence")
>Edit the header of EGFP fasta file tp make it more informatice
```bash
>EGFP
TAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAA
CTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATG
TTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCA
CTTGGCAGTACATCAAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCC
GCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCA
TCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGG
ATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCA
AAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAA
GCAGAGCTGGTTTAGTGAACCGTCAGATCCGCTAGCGCTACCGGACTCAGATCTCGAGCTCAAGCTTCGA
ATTCTGCAGTCGACGGTACCGCGGGCCCGGGATCCACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGC
TGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTC
```
We need to count the number of bases in the sequence:
```bash
cat EGFP.fa | grep -v "^>" | tr -d "\n" | wc -c
```
The results of this command shows there are 4733 bases. This is important to know for the creating custom GTF for EGFP.

Now to make a custom GTF for EGFP with the following command:
>Note: we need to insert the tabs that separate the 9 columns of information required for GTF.
```bash
echo -e 'EGFP\tunknown\texon\t1\t4733\t.\t+\t.\tgene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";' > EGFP.gtf
```

This is what the EGFP.gtf file looks like with the `cat EGFP.gtf` command:
```bash
EGFP	unknown	exon	1	4733	.	+	.	gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
```

Next, add the EGFP.fa to the end of the M. Musculus genome FASTA. But first, make a copy so that the original is unchanged.
```bash
cp Mus_musculus.filtered_gtf.gtf Mus_musculus.filtered_gtf_EGFP.gtf
```

Now append `EGFP.fa` to `Mus_musculus.filtered_gtf_EGFP.gtf` as following:
```bash
cat EGFP.fa >> Mus_musculus.filtered_gtf_EGFP.gtf
```

###Use cellranger mkref command to create custom reference directory with EGFP added
Now use the genome_Mus_musculus_EGFP.fa and Mus_musculus.filtered_gtf_EGFP.gtf files as inputs to the cellranger mkref pipeline:
```bash
#! /bin/bash
#SBATCH -A LSENS2018-3-3 # the ID of our Aurora project
#SBATCH -n 20 # how many processor cores to use
#SBATCH -N 1 # how many processors to use
#SBATCH -t 24:00:00 # kill the job after ths hh::mm::ss time
#SBATCH -J addGFP_Mus_musculus # name of the job
#SBATCH -o addGFP_Mus_musculus%j.out # stdout log file
#SBATCH -e addGFP_Mus_musculus%j.err # stderr log file
#SBATCH -p dell # which partition to use
module purge
module load cellranger/6.0.0
cellranger mkref --genome=Mus.musculus_genome_EGFP \
--fasta=genome_Mus_musculus_EGFP.fa \
-genes=Mus_musculus.filtered_gtf_EGFP.gtf
exit 0
```
This outputs a custom reference directory called `Mus.musculus_genome_EGFP/`.


