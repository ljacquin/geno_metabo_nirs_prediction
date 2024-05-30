#!/bin/sh
#=====================================================#
#  script for launching genomic prediction for traits #
#=====================================================#
### Requirements
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64
#SBATCH --cpus-per-task=20

### Email
#SBATCH --mail-user=laval.jacquin@inrae.fr
#SBATCH --mail-type=ALL

R -q --vanilla < genomic_prediction_and_analysis_BLUPaudpc_eau_EUB04_18.R