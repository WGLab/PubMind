#!/bin/bash
#SBATCH --time=5-00
#SBATCH --array=1-8
#SBATCH --output=/mnt/isilon/wang_lab/pengwang/projects/LLM/PubMind/run_log/slurm-%A_%a.out
#SBATCH --error=/mnt/isilon/wang_lab/pengwang/projects/LLM/PubMind/run_log/slurm-%A_%a.err
#SBATCH -c 6
#SBATCH --mem=12G
#SBATCH -p gpu-xe9680q
#SBATCH --gres=gpu:h100:2


#source activate vllm

mkdir -p /mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/run_log

batch_files=(
    "run_all_pmc_01"
    "run_all_pmc_02"
    "run_all_pmc_03"
    "run_all_pmc_04"
    "run_all_pmc_05"
    "run_all_pmc_06"
    "run_all_pmc_07"
    "run_all_pmc_08"
)


# Get the batch file name based on the array task ID
batch_file=${batch_files[$SLURM_ARRAY_TASK_ID - 1]}

echo "Running job $batch_file"
python /mnt/isilon/wang_lab/pengwang/projects/LLM/PubMind/pmc_fulltext-llm_inference.py "$batch_file"
