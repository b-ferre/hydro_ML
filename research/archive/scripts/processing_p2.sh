#!/bin/bash
#SBATCH  --array=548-1096:20
#SBATCH  --account=def-aliameli
#SBATCH  --time=30:00:00
#SBATCH  --mem-per-cpu=25G
#SBATCH  --cpus-per-task=25
#SBATCH  --job-name='processing-2'

echo "running on 25x25G..."
module load gcc/9.3.0 StdEnv/2020 udunits/2.2.28 gdal/3.5.1 proj/9.0.1 geos/3.10.2 r/4.2.2
Rscript ../R/processing.R $SLURM_ARRAY_TASK_ID