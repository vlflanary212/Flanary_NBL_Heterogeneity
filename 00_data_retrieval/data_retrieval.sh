#!/bin/bash

##### SLURM #####
#SBATCH --job-name=nbatlas
#SBATCH --partition=express
#SBATCH --time=1:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=flanary@uab.edu
#SBATCH --mail-type=ALL


##### VARIABLES #####
data_dir="/home/flanary/Dry_Lab/NB_ITH/NBAtlas/data"
subset_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/0a569381-3a0c-4eec-863a-e20544b686ed/file_downloaded"
alldata_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/e970f25e-de9f-4956-8d18-7dfa2e8bb6fe/file_downloaded"

##### COMMANDS #####
# Create directories if they don't exist
mkdir -p "$data_dir/subset" "$data_dir/alldata"

# Download the 50k cell subset
echo "Downloading 50k cell subset..."
cd "$data_dir/subset"
wget -O 00_init_obj.rds "$subset_url" && echo "50k subset downloaded successfully."

# Download the entire dataset
echo "Downloading entire dataset..."
cd "$data_dir/alldata"
wget -O 00_init_obj.rds "$alldata_url" && echo "Entire dataset downloaded successfully."

##### END #####
echo "NBAtlas download completed."