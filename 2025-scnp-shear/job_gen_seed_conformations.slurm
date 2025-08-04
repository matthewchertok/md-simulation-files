#!/bin/bash
#SBATCH --job-name=generate_seed_conformations
#SBATCH --output=slurm_outputs/init_conformations/%A_%a.out
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=12G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --gres=gpu:1             # number of GPUs
#SBATCH --time=01:02:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-1079%40      # no more than 40 GPUs at once to avoid using the whole cluster (48)     
#SBATCH --mail-user=mc7345@princeton.edu 
module purge
source ~/.bashrc
module load anaconda3/2024.6
conda activate '/home/mc7345/micromamba/envs/matthew'
 
python main_generate_seed_conformations.py