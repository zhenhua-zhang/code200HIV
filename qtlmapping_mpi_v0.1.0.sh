#!/bin/bash
#SBATCH --time=0:19:0
#SBATCH --mem=10G
#SBATCH --cpus=7
#SBATCh --nodes=4
#SBATCH --output=%j-%u-qtlmapping_mpi_v0.1.0.log
#SBATCH --job-name=qtlmapping_mpi

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

~/tools/bin/mpirun -np 5 Rscript ./qtlmapping_mpi_v0.1.0.R \
    --run-flag padding4_mpi_test \
    --draw-pwcor \
    --trps-cvrt-dtfm \
    --trps-pntp-dtfm \
    --work-dir ${o_dir}/HIVReservior/HIVReservior_AGC4nHd/padding4_mpi_test \
    --pntp-file ${i_dir}/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv \
    --target-pntp RNAHIV_CD4LOG,DNAHIV_CD4LOG,RNAvsDNA_CD4LOG \
    --cvrt-file ${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
    --target-cvrt age,gender,HIV_DURATION,CD4_NADIR \
    --gntp-dosage-file ${p_dir}/inputs/dosages/MatrixEQTL/200HIV_dosage.gz \
    --gntp-info-file ${p_dir}/inputs/dosages/MatrixEQTL/200HIV_variantInfo.gz \
    --pm-times 1000 \
    --mhtn-fig-p-thrd 0.01

