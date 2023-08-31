#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=gather_genomes_6.17
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=schlenke
### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --nodes=1
#SBATCH --ntasks=30
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=100gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=40:00:00

#USAGE: sbatch gather_genomes.slurm NCBI_taxID_# output_dir_path #threads
bash ~/Path_to_gits/OrthoPhyl/utils/gather_filter_asms.restart.sh $1 $2 $3
