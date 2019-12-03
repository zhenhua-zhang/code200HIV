#!/bin/bash
#SBATCH --time=0:20:0
#SBATCH --cpus=1
#SBATCH --mem=10G
#SBATCH --output=%j-%u-check_cellCounts.log

set -o errexit
set -o errtrace

module load R/3.5.1-foss-2015b-bare
module list

p_dir="~/Documents/projects/200HIV"
o_dir="~/Documents/projects/200HIV/outputs/HIVReservior/200HIV_phenotypeLevelPerGenotype_rs7817587/cellCounts_rs7817589"

Rscript check_genotype_distribution_v0.1.0.R \
    --scale-phenotype \
    --use-genotype-symbols \
    --genotype-dosage-file ${p_dir}/inputs/dosage/200HIV_dosages/chr8_dosage.gz \
    --genotype-info-file ${p_dir}/inputs/dosage/200HIV_dosages/chr8_variantInfo.gz \
    --phenotype-file ${p_dir}/inputs/datasets/totalDataHans200HivWithPercent.csv \
    --covariate-file ${p_dir}/inputs/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
    --target-covariates age,gender,CD4_NADIR,HIV_DURATION \
    --target-snp rs7817589 \
    --output-dir ${o_dir} \
    1>cellCounts_rs7817589_ttest.txt \
    && echo "[INFO] Job exit with $?" || echo "[ERROR] Job exit with $?"
