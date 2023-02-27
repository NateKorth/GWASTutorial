#!/bin/sh
#SBATCH --time=167:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name=RunR
#SBATCH --error=RunR.err
#SBATCH --output=RunR.out
#SBATCH --nodes=1

ml R/4.0

R CMD BATCH CalculateBLUEs.R

mv Phenotypes.csv ../input
