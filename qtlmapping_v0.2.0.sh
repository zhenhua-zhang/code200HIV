#!/bin/bash
#SBATCH --time=0:10:0
#SBATCH --mem=10G
#SBATCH --cpus=1
#SBATCH --output=%j-%u-qtlmapping_v0.2.0.log
#SBATCH --job-name=%j-%u-qtlmapping

set -o errexit
set -o errtrace

module load R/3.5.1-foss-2015b-bare
module list

p_dir=/groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV
i_dir=${p_dir}/inputs

Rscript ./qtlmapping_v0.2.0.R \
    --draw-pwcor \
    --trps-cvrt-dtfm \
    --trps-pntp-dtfm \
    --work-dir test_run \
    --pntp-file ${i_dir}/datasets/20190524_HIVreservoir_GENT.csv \
    --target-pntp DNAHIV_CD4 \
    --cvrt-file ${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
    --gntp-dosage-file ${i_dir}/dosage/200HIV_dosages/chr22_dosage.gz \
    --gntp-info-file ${i_dir}/dosage/200HIV_dosages/chr22_variantInfo.gz


