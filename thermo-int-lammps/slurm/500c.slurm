#!/bin/bash
#SBATCH --job-name=500c        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=96     # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=02:00:00          # total run time limit (HH:MM:SS)

module load intel/19.1.1.217
module load intel-mpi/intel/2019.7

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

lammpsexe=/LAMMPS/EXECUTABLE/FILE

cd ../500
for i in {1..12}
do
    cd $i
    srun --nodes=1 --ntasks=8 --exclusive $lammpsexe -in c.in > $i.out &
    cd ..
done
cd ..

wait
