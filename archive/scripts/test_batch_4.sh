#!/bin/bash
#SBATCH  --array=1-10:2
#SBATCH  --account=def-aliameli
#SBATCH  --time=1:00:00
#SBATCH  --mem-per-cpu=50G
#SBATCH  --cpus-per-task=5
#SBATCH  --job-name='test_batch_2'

echo "running on 5x50G..."
module load gcc/9.3.0 StdEnv/2020 udunits/2.2.28 gdal/3.5.1 proj/9.0.1 geos/3.10.2 r/4.2.2
Rscript ../R/processing.R $SLURM_ARRAY_TASK_ID 
