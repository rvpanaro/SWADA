#!/bin/bash
# Submit the published simulation from the simulation/ directory:
#   sbatch fastR/r_job.sh
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/main_%j.out
#SBATCH --error=logs/main_%j.err

module load gcc/14.2.0
module load r/4.4.1
module load cmake
module load gsl
module load intel-oneapi-tbb

cd "$(dirname "$0")/.."
export SIM_MODE="${SIM_MODE:-paper}"
Rscript fastR/main.R
