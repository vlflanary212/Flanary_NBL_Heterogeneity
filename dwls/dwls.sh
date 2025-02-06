#!/bin/bash

##### SLURM #####
#SBATCH --job-name=dwls
#SBATCH --partition=medium,long,largemem
#SBATCH --time=49:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL

##### PACKAGES ######
module load Singularity

##### VARIABLES #####
sif="/home/flanary/sif_files/deconvolution_2.4.sif"
script="/home/flanary/Projects/Dry_Lab/NB_ITH/PCGP_RNA/dwls/dwls.R"

##### COMMANDS #####
singularity exec $sif Rscript $script

##### END #####
echo "done"
