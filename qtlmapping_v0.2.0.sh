#!/bin/bash
#SBATCH --time=0:59:0
#SBATCH --mem=100G
#SBATCH --cpus=1
#SBATCH --output=%j-%u-qtlmapping_v0.2.0.log
#SBATCH --job-name=%j-%u-qtlmapping

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc
module load R/3.5.1-foss-2015b-bare
module list

p_dir=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV
i_dir=${p_dir}/inputs
o_dir=${p_dir}/outputs

Rscript ./qtlmapping_v0.2.0.R \
    --run-flag DNAandRNA_RNAvsDNA_padding4_newGenotype \
    --draw-pwcor \
    --trps-cvrt-dtfm \
    --trps-pntp-dtfm \
    --work-dir ${o_dir}/HIVReservior/HIVReservior_AGC4nHd \
    --pntp-file ${i_dir}/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv \
    --target-pntp RNAHIV_CD4LOG,DNAHIV_CD4LOG,RNAvsDNA_CD4LOG \
    --cvrt-file ${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
    --target-cvrt age,gender,HIV_DURATION,CD4_NADIR \
    --gntp-dosage-file ${p_dir}/inputs/dosages/200HIV_dosage.gz \
    --gntp-info-file ${p_dir}/inputs/dosages/200HIV_variantInfo.gz \
    --mhtn-fig-p-thrd 0.01

