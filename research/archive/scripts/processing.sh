#!/bin/bash
#SBATCH  --array=1-2192:2
#SBATCH  --account=def-aliameli
#SBATCH  --time=3:00:00
#SBATCH  --mem-per-cpu=30G
#SBATCH  --cpus-per-task=4
#SBATCH  --job-name='hail_mary_gamma'

echo "running on 4x30G..."
module load gcc/9.3.0 StdEnv/2020 udunits/2.2.28 gdal/3.5.1 proj/9.0.1 geos/3.10.2 r/4.2.2
Rscript ../R/processing.R $SLURM_ARRAY_TASK_ID