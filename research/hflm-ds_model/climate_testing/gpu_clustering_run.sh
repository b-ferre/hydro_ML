#!/bin/bash
#SBATCH  --account=def-aliameli
#SBATCH --job-name=gpu_clustering
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem= 200G

# Load modules
module load python/3.8
module load cuda/11.2
module load cudnn/8.1
module load gcc/9.3.0

# Run the Python script
python gpu_clustering.py