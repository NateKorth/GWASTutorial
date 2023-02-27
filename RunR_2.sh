#!/bin/sh
#SBATCH --time=80:00:00
#SBATCH --mem-per-cpu=60000
#SBATCH --job-name=GWASPrep
#SBATCH --error=GWASPrep.err
#SBATCH --output=GWASPrep.out
#SBATCH --nodes=1

ml R/4.0
R CMD BATCH PrepGenoForMVP.R
