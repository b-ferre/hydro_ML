#!/bin/bash
#SBATCH  --account=def-aliameli
#SBATCH  --time=30:00
#SBATCH  --cpus-per-task=10
#SBATCH  --mem=300G
#SBATCH  --job-name='clustering_test'

module load StdEnv/2023 r/4.4.0
Rscript ./clustering.R