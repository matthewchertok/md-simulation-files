#!/bin/bash
#SBATCH --job-name=production_irreversible_bonding
#SBATCH --output=slurm_outputs/production/irreversible_bonding_%A_%a.out
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=12G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --gres=gpu:1             # number of GPUs
#SBATCH --time=03:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-10799%30      # no more than 30 GPUs at once to avoid using the whole cluster (48). We want to go up to index 10799, which will be the last index of the 10800 jobs we need to run     
#SBATCH --mail-user=mc7345@princeton.edu
module purge
source ~/.bashrc
module load anaconda3/2024.6
conda activate '/home/mc7345/micromamba/envs/matthew'

# irreversible bonding
python main_generate_production_data.py