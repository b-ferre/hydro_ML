#!/bin/bash
#SBATCH  --array=1-1528:2
#SBATCH  --account=def-aliameli
#SBATCH  --time=120:00
#SBATCH  --mem-per-cpu=10G
#SBATCH  --cpus-per-task=2
#SBATCH  --job-name='climate_testing'

module load StdEnv/2023 r/4.4.0
Rscript ./cc_inference.R $SLURM_ARRAY_TASK_ID