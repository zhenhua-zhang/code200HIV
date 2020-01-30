#!/bin/bash
#SBATCH --time=4:59:0
#SBATCH --mem=15G
#SBATCH --cpus=15
#SBATCH --output=%A_%a-%u-qtlmapping_permutation_mpi_v0.1.0.log
#SBATCH --job-name=qtlmapping_permutation_mpi
#SBATCH --array=1-2

#
## NOTE: 
## 1. The R doesn't support large matrices of which any demension is greater
## than 2^31-1. Provided the number of common SNPs is more than 4 millions, it's
## not applicable to use pbdMPI for some funciton, e.g cumsum, doesn't support 
## vector longer than 2^31 - 1. The main reasons are because of the the linear
## algbra libraries by Fortran are compiled by `INTEGER` 32
## 2. Set HWLOC_HIDE_ERRORS to 1 to suppress an error.
## 3. Flag --oversubscribe to mimic slots for OpenMPI.
#

set -o errexit
set -o errtrace

export HWLOC_HIDE_ERRORS=1

source /apps/modules/modules.bashrc
module load OpenMPI/1.8.8-GNU-4.9.3-2.25
module load R/3.5.1-foss-2015b-bare
module list

p_dir=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV
i_dir=${p_dir}/inputs
o_dir=${p_dir}/outputs
chromId=${SLURM_ARRAY_TASK_ID:=1}

mpirun -np 15 --oversubscribe \
    Rscript ./qtlmapping_mpi_v0.1.0.R \
    --run-flag padding4_mpi_${chromId} \
    --trps-cvrt-dtfm \
    --trps-pntp-dtfm \
    --work-dir ${o_dir}/HIVReservior/HIVReservior_AGC4nHd/padding4_mpi \
    --pntp-file ${i_dir}/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv \
    --target-pntp RNAHIV_CD4LOG,DNAHIV_CD4LOG,RNAvsDNA_CD4LOG \
    --cvrt-file ${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
    --target-cvrt age,gender,HIV_DURATION,CD4_NADIR \
    --gntp-dosage-file ${p_dir}/inputs/dosages/MatrixEQTL/chr${chromId}_dosage.gz \
    --gntp-info-file ${p_dir}/inputs/dosages/MatrixEQTL/chr${chromId}_variantInfo.gz \
    --pm-times 6000

