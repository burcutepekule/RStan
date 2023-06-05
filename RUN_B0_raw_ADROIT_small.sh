#!/bin/sh
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=bt6725@princeton.edu

module load anaconda3/2023.3
eval "$(conda shell.bash hook)"
conda activate RStan01

srun --nodes=1 --ntasks=1 --cpus-per-task=4 --mem-per-cpu=1G --time=04:00:00 --job-name=model-raw-small --exclusive \
  Rscript --vanilla B0_MODEL_RUN_ADROIT.R config_raw_small.txt 


