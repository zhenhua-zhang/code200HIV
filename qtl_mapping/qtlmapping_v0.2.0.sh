#!/bin/bash
#SBATCH --time=0:59:0
#SBATCH --mem=55G
#SBATCH --cpus=1
#SBATCH --output=%j-%u-qtlmapping_v0.2.0.log
#SBATCH --job-name=qtlmapping_v0.2.0

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc
module load R/3.5.1-foss-2015b-bare
module list

p_dir=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV
i_dir=${p_dir}/inputs
o_dir=${p_dir}/outputs
run_flag=HIVDNA_P4AGHdCnHrna # HIVDNA_padding4_age_gender_HIVDuration_CD4Nadir_HIVRNALOG4

echo ${run_flag}

    # --target-cvrt age,gender,HIV_DURATION,CD4_NADIR,cARTDuration,CD4Latest
    # --cvrt-file ${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv
    # --apply-log10-on-cvrt age,cARTDuration \
Rscript ./qtlmapping_v0.2.0.R \
    --padding 4 \
    --run-flag ${run_flag} \
    --draw-pwcor \
    --trps-cvrt-dtfm \
    --trps-pntp-dtfm \
    --work-dir ${o_dir}/HIVReservior/HIVReservior_${run_flag} \
    --pntp-file ${i_dir}/datasets/20190524_HIVreservoir_GENT_HIVDNA.tsv \
    --target-pntp DNAHIV_CD4LOG \
    --cvrt-file ${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_HIVRNA_20200302.csv \
    --target-cvrt age,gender,HIV_DURATION,CD4_NADIR,RNAHIV_CD4LOG \
    --gntp-dosage-file ${p_dir}/inputs/dosages/MatrixEQTL/200HIV_dosage.gz \
    --gntp-info-file ${p_dir}/inputs/dosages/MatrixEQTL/200HIV_variantInfo.gz \
    --mhtn-fig-p-thrd 0.01

