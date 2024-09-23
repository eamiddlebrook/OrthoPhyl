#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=gather_genomes_8.20.24
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
#SBATCH --partition=fast
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --nodes=1
#SBATCH --ntasks=28
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=100gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=40:00:00

#USAGE: sbatch gather_genomes.slurm NCBI_taxID_# local_output_dir #threads

bash ~/gits/OrthoPhyl/utils/gather_filter_asms.sh $1 $2 $3
