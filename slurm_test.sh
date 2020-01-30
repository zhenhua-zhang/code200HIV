#!/bin/bash
#SBATCH --time=0:1:0
#SBATCH --mem=200M
#SBATCH --cpus=4
#SBATCH --output=slurm_test.log

source /apps/modules/modules.bashrc
module load OpenMPI/1.8.8-GNU-4.9.3-2.25

mpirun hostname > tmp.txt

mpirun -np 4 --host umcg-node011 --oversubscribe ls
