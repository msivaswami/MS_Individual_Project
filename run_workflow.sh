#!/bin/bash
#SBATCH --account=weberrj01
#SBATCH --nodes 1
#SBATCH --ntasks 5
#SBATCH --time 6:00:0
#SBATCH --qos bbdefault
#SBATCH --mem-per-cpu 4096M
#SBATCH --mail-type ALL

# bluebear modules
set -e

module purge
module load bear-apps/2022b
module load Miniconda3/22.11.1-1

env_path=/rds/projects/2015/viantm-01/users/weberrj/run-mspurity-metfrag/workflows/envs
env_name=mspurity_metfrag_workflow_UNIX

CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh  # bb

# Activate and deactivate conda environment
conda activate "${env_path}/${env_name}"


# Rscript 1_run_xcms_mspurity_workflow.R
${env_path}/${env_name}/bin/python 2_metfrag.py

conda deactivate
