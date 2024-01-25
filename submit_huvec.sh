#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=RE_huvec
#SBATCH --ntasks=4
#SBATCH --mem=80G
#SBATCH --mail-user=simone.riva@imm.ox.ac.uk
#SBATCH --time=06-23:59:59
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

source /project/Wellcome_Discovery/sriva/mambaforge/bin/activate re

# select and edit the config file accordingly
snakemake --configfile=config/huvec/huvec.yaml all --cores 6 --unlock
snakemake --configfile=config/huvec/huvec.yaml all --cores 6 --rerun-incomplete