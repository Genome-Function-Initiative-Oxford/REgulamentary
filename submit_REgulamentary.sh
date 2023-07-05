#!/bin/sh

#SBATCH --job-name=RE_3
#SBATCH --time=1-02:00:00
#SBATCH --ntasks=5
#SBATCH --mem=60G
#SBATCH --partition=jumbo
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

#source /ceph/project/Wellcome_Discovery/sriva/mambaforge/bin/activate upstream
# conda activate upstream

snakemake --configfile=config/paper_huvec.yaml all --cores 3 --unlock
snakemake --configfile=config/paper_huvec.yaml all --cores 3 --rerun-incomplete

snakemake --configfile=config/paper_keratinocyte.yaml all --cores 3 --unlock
snakemake --configfile=config/paper_keratinocyte.yaml all --cores 3 --rerun-incomplete