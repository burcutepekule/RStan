#!/usr/bin/env bash
#SBATCH --job-name=model-raw     # create a short name for your job
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=bt6725@princeton.edu

module load anaconda3/2023.3
conda activate RStan01

srun --nodes=1 --ntasks=1 --cpus-per-task=1 --exclusive \
  Rscript --vanilla B0_MODEL_RUN.R config_raw.txt 


