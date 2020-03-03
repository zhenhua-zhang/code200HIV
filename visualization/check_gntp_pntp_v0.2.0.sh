#!/bin/bash
#SBATCH --time=0:59:0
#SBATCH --cpus=1
#SBATCH --mem=20G
#SBATCH --output=%j-%u-check-gntp-pntp-v0.2.0.log

set -o errexit
set -o errtrace
source /apps/modules/modules.bashrc

module load R/3.5.1-foss-2015b-bare
module list

p_dir=~/Documents/projects/200HIV
i_dir=${p_dir}/inputs
o_dir=${p_dir}/outputs

# Input
genotype_dosage_file=${i_dir}/dosages/MatrixEQTL/200HIV_dosage.gz
genotype_info_file=${i_dir}/dosages/MatrixEQTL/200HIV_variantInfo.gz
phenotype_file=${i_dir}/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv
covariate_file=${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv
target_covariates=age,gender,CD4_NADIR,HIV_DURATION
target_snp=rs7113204,rs7817589,rs7113255,rs2613996,rs12366210

# Outputs
output_dir=${o_dir}/HIVReservior/HIVReservior_AGC4nHd/significant_QTLs/$target_snp
mkdir -p ${output_dir}

Rscript check_gntp_pntp_v0.2.0.R \
    --gntp-dosage-file $genotype_dosage_file \
    --gntp-info-file $genotype_info_file \
    --pntp-file $phenotype_file \
    --target-pntp RNAvsDNA_CD4LOG,DNAHIV_CD4LOG \
    --cvrt-file $covariate_file \
    --target-cvrt $target_covariates \
    --target-snp $target_snp \
    --output-dir $output_dir \
    && echo "[INFO] Job exit with $?" 1>&2 || echo "[ERROR] Job exit with $?" 1>&2
