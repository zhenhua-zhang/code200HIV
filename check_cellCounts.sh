#!/bin/bash
#SBATCH --time=0:20:0
#SBATCH --cpus=1
#SBATCH --mem=10G
#SBATCH --output=%j-%u-check_cellCounts.log

set -o errexit
set -o errtrace

module load R/3.5.1-foss-2015b-bare
module list

p_dir=~/Documents/projects/200HIV
i_dir=${p_dir}/inputs
o_dir=${p_dir}/outputs

# Input
genotype_dosage_file=${i_dir}/dosage/200HIV_dosages/chr8_dosage.gz
genotype_info_file=${i_dir}/dosage/200HIV_dosages/chr8_variantInfo.gz
phenotype_file=${i_dir}/datasets/totalDataHans200HivWithPercent.csv
covariate_file=${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv
target_covariates=age,gender,CD4_NADIR,HIV_DURATION
target_snp=rs7817589

# Outputs
glm_summary_file=${o_dir}/HIVReservior/200HIV_phenotypeLevelPerGenotype_rs7817587/cellCounts_rs7817589/cellCounts_rs7817589_glm.txt
glm_summary_sig_file=${o_dir}/HIVReservior/200HIV_phenotypeLevelPerGenotype_rs7817587/cellCounts_rs7817589/cellCounts_rs7817589_glm_gntp0.05.txt
output_dir=${o_dir}/HIVReservior/200HIV_phenotypeLevelPerGenotype_rs7817587/cellCounts_rs7817589

mkdir -p ${output_dir}/cellCountsAll ${output_dir}/cellCounts0.05

Rscript check_genotype_distribution_v0.2.0.R \
    --genotype-dosage-file ${genotype_dosage_file} \
    --genotype-info-file ${genotype_info_file} \
    --phenotype-file ${phenotype_file} \
    --covariate-file ${covariate_file} \
    --target-covariates ${target_covariates} \
    --target-snp ${target_snp} \
    --output-dir ${output_dir}/cellCountsAll \
    > ${glm_summary_file} \
    && echo "[INFO] Job exit with $?" 1>&2 || echo "[ERROR] Job exit with $?" 1>&2

# Collect the GLM results of which genotype is significant.
for line in $(grep -n gntp ${glm_summary_file} | grep "\*" | cut -d':' -f1); do
    start=$[ $line - 12 ]
    stop=$[ $line + 16 ]
    sed -n "${start},${stop}p" ${glm_summary_file}
done > ${glm_summary_sig_file}

# Collect the box plot of which genotype is significant in GLM regression
for x in $(grep "<<<" ${glm_summary_sig_file}  | cut -d $':' -f 2 | tr -d " "); do
    cp ${output_dir}/cellCountsAll/*${x}_${target_snp}.pdf ${output_dir}/cellCounts0.05
done