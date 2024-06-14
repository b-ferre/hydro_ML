#!/bin/bash
#SBATCH  --account=def-aliameli
#SBATCH  --time=180:00
#SBATCH  --mem-per-cpu=500G
#SBATCH  --cpus-per-task=1
#SBATCH  --job-name='clustering'

module load StdEnv/2023 r/4.4.0
Rscript ./clustering.R