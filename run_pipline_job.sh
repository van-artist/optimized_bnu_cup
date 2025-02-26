#!/bin/bash
#SBATCH --job-name=pipline
#SBATCH --output=./out/pipline/%x_%j.out
#SBATCH --error=./out/pipline/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=G1Part_sce

source /es01/paratera/parasoft/module.sh  
source /es01/software/oneapi/vtune/latest/vtune-vars.sh

conda init
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

bash ./run_pipline.sh