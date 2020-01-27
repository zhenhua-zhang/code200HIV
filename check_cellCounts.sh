#!/bin/bash
#SBATCH --time=0:59:0
#SBATCH --cpus=1
#SBATCH --mem=20G
#SBATCH --output=%j-%u-check_cellCounts.log

set -o errexit
set -o errtrace
source /apps/modules/modules.bashrc

module load R/3.5.1-foss-2015b-bare
module list

p_dir=~/Documents/projects/200HIV
i_dir=${p_dir}/inputs
o_dir=${p_dir}/outputs

# Input
genotype_dosage_file=${i_dir}/dosage/200HIV_dosage.gz
genotype_info_file=${i_dir}/dosage/variantInfo.gz
phenotype_file=${i_dir}/datasets/totalDataHans200HivWithPercent.csv
covariate_file=${i_dir}/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv
target_covariates=age,gender,CD4_NADIR,HIV_DURATION

top_snp=rs7113204
if [ $top_snp == "rs7817589" ]; then # rs7817589 p-vaue < 5e-7
    target_snp=rs7817589,rs4739797,rs1834712,rs4523210,rs10504736,rs4357250,rs7818093,rs9298346,rs1427075,rs6473264,rs12542001,rs1427073,rs1863648,rs2081675,rs7015085,rs10958004,rs12548276
elif [ $top_snp == "rs7113204" ]; then #rs7113204 p-value < 4e-7
    target_snp=rs7113204,rs12576389,rs7396778,rs146730949,rs12577368,rs6598020,rs113061269,rs55768561,rs72851107,rs7110722,rs2613996,rs56356502,rs111512069,rs11827672,rs67092853
fi

# Outputs
output_dir=${o_dir}/HIVReservior/cellCountsAnalysis/$top_snp

mkdir -p ${output_dir}/cellCounts{All,0.05}

Rscript check_gntp_pntp_v0.2.0.R \
    --genotype-dosage-file $genotype_dosage_file \
    --genotype-info-file $genotype_info_file \
    --phenotype-file $phenotype_file \
    --covariate-file $covariate_file \
    --target-covariates $target_covariates \
    --target-snp $target_snp \
    --output-dir $output_dir/cellCountsAll \
    && echo "[INFO] Job exit with $?" 1>&2 || echo "[ERROR] Job exit with $?" 1>&2

for snp in $(echo $target_snp | tr "," " "); do
    input_file=$output_dir/cellCounts_$snp.txt
    output_file=${input_file/$snp/$snp.0.05}

    # Collect the GLM results of which genotype is significant.
    for line in $(grep -n gntp\  $input_file | grep "\*" | cut -d':' -f1); do
        start=$[ $line - 12 ]
        stop=$[ $line + 16 ]
        sed -n "${start},${stop}p" $input_file
    done > $output_file

    # Collect the box plot of which genotype is significant in GLM regression
    for x in $(grep "<<< Summary" $output_file  | cut -d $':' -f 2 | tr -d " "); do
        cp $output_dir/cellCountsAll/*_${x}_$snp.pdf $output_dir/cellCounts0.05
    done
done

rm -fr $output_dir/cellCountsAll/*
