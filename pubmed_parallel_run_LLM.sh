#!/bin/bash
#SBATCH --time=5-00
#SBATCH --array=1-4
#SBATCH --output=/mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/run_log/slurm-%A_%a.out
#SBATCH --error=/mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/run_log/slurm-%A_%a.err
#SBATCH -c 6
#SBATCH --mem=12G
#SBATCH -p gpuq
#SBATCH --gres=gpu:a100:2


#source activate vllm

mkdir -p /mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/run_log

batch_files=(
    "run_pubmed_01"
    "run_pubmed_02"
    "run_pubmed_03"
    "run_pubmed_04"
)

# Get the batch file name based on the array task ID
batch_file=${batch_files[$SLURM_ARRAY_TASK_ID - 1]}

echo "Running job $batch_file"
python /mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/pubmed_abstract_extract-automation.py "$batch_file"