#!/bin/bash
#SBATCH --partition=<partition-name>
#SBATCH --job-name=HUVEC
#SBATCH --ntasks=2
#SBATCH --mem=64G
#SBATCH --mail-user=<your-email>
#SBATCH --time=03-12:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

source /path/to/baseenv/bin/activate re

# select and edit the config file accordingly
snakemake --configfile=config/endothelial-cell-of-umbilical-vein.yaml all --cores 2 --unlock
snakemake --configfile=config/endothelial-cell-of-umbilical-vein.yaml all --cores 2 --rerun-incomplete