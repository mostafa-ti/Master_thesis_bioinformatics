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
 job scheduler, which in simple words, is a bash script that loads all the necessary packages and sends it to the backend system.
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

