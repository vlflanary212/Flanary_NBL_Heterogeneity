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
module load Anaconda3
conda activate r-env

##### VARIABLES #####
script="/data/project/sen-lab/Flanary_NBL_Heterogeneity/DWLS/dwls.R"

##### COMMANDS #####
Rscript $script

##### END #####
echo "done"
