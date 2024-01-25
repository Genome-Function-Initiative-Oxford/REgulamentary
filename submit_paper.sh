#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=RE_1
#SBATCH --ntasks=4
#SBATCH --mem=80G
#SBATCH --mail-user=simone.riva@imm.ox.ac.uk
#SBATCH --time=06-23:59:59
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

source /project/Wellcome_Discovery/sriva/mambaforge/bin/activate re

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/astrocyte.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/astrocyte.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/B-cell.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/B-cell.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/cardiac-muscle-cell.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/cardiac-muscle-cell.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/CD14-positive-monocyte.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/CD14-positive-monocyte.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/CD4-positive.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/CD4-positive.yaml all --cores 6 --rerun-incomplete

# select and edit the config file accordingly
snakemake --configfile=config/paper/endothelial-cell-of-umbilical-vein.yaml all --cores 6 --unlock
snakemake --configfile=config/paper/endothelial-cell-of-umbilical-vein.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/fibroblast-of-dermis.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/fibroblast-of-dermis.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/fibroblast-of-lung.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/fibroblast-of-lung.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/keratinocyte.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/keratinocyte.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/mammary-epithelial-cell.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/mammary-epithelial-cell.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/natural-killer-cell.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/natural-killer-cell.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/osteoblast.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/osteoblast.yaml all --cores 6 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/paper/skeletal-muscle-myoblast.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/skeletal-muscle-myoblast.yaml all --cores 6 --rerun-incomplete


# # select and edit the config file accordingly
# snakemake --configfile=config/paper/osteoblast.yaml all --cores 6 --unlock
# snakemake --configfile=config/paper/osteoblast.yaml all --cores 6 --rerun-incomplete

##### Ery

# # select and edit the config file accordingly
# snakemake --configfile=config/erythroid/Don001.yaml all --cores 3 --unlock
# snakemake --configfile=config/erythroid/Don001.yaml all --cores 3 --rerun-incomplete

# # select and edit the config file accordingly
# snakemake --configfile=config/erythroid/Don002.yaml all --cores 3 --unlock
# snakemake --configfile=config/erythroid/Don002.yaml all --cores 3 --rerun-incomplete