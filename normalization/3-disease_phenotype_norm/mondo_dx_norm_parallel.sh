#!/bin/bash
#SBATCH --time=5-00
#SBATCH --array=0-3  # This creates 4 jobs (job IDs: 0, 1, 2, 3)
#SBATCH --output=dx_norm_slurm_log/job_output_%A_%a.out
#SBATCH --error=dx_norm_slurm_log/job_error_%A_%a.err
#SBATCH -c 6
#SBATCH --mem=16G
#SBATCH -p gpuq
#SBATCH --gres=gpu:1

#conda activate vllm 

mkdir -p dx_norm_slurm_log

# Run Python script with job index
python llm_dx_norm_mondo.py $SLURM_ARRAY_TASK_ID
