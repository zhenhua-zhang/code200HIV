#!/bin/bash
#SBATCH --job-name=QtlMappingPlink
#SBATCH --output=%j-%u-QtlMappingPlink.log
#SBATCH --time=0:59:50
#SBATCH --cpus=1
#SBATCH --mem=10G

source /apps/modules/modules.bashrc
set -o errexit
set -o errtrace
# # Preprocessing VCF files
# module load BCFtools
# module list
# 
# # Filter by Xiaojing Chu
# # bcftools view -i 'INFO/MAF[0] > 0.1 && (INFO/R2[0] > 0.3 || INFO/ER2[0] > 0.3)'
# 
# # Concate VCF file per chromosome.

# # NOTE: Some SNP misses its ID, it's not legal for plink, one should annotate
# # missing SNP id using bcftools annotate. Or set the missing id in Plink using
# # --set-missing-var-ids, when not using --dosage option
cd /groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV/inputs/dosages/Plink
# bcftools concat --threads 5 \
#     /groups/umcg-wijmenga/tmp04/umcg-zzhang/projects/200HIV/annotated/*vcf.gz \
#     | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
#     -O z \
#     -o 200HIV_annotated.vcf.gz
#  
# # Rename the sample names
# bcftools view \
#     -h 200HIV_annotated.vcf.gz \
#     | sed -n -e 's/v_[0-9]\{1,4\}_X/X/gp' -e "/##/p"\
#     > 200HIV_annotated_new_header.txt
# 
# bcftools reheader \
#     -h 200HIV_annotated_new_header.txt \
#     -o 200HIV_annotated_reheader.vcf.gz \
#     200HIV_annotated.vcf.gz
# 
# rm 200HIV_annotated_new_header.txt -fr
# mv 200HIV_annotated_reheader.vcf.gz 200HIV_annotated.vcf.gz -f
# 
# bcftools index \
#     -tf \
#     --threads 5 \
#     200HIV_annotated.vcf.gz
# 
# # DosageConvertor: Convert VCF into Plink dosage
# # NOTE: DosageConvertor is single-thread, so per chrom converting could be
# # better choice.
# module purge
# module load DosageConvertor
# echo "Modules used for DosageConvertor:"
# module list
# 
# DosageConvertor \
#     --vcfDose 200HIV_annotated.vcf.gz \
#     --prefix 200HIV_annotated \
#     --format 1 \
#     --type plink

# # Add sex informatoin:
# # No X chromosome in the VCF file, then one needs update the gender information
# # manually.
# 
# # Create a file including phenotypes, 
# # FID IID HIVDNA_CD4LOG HIVRNA_CD4LOG
# # X1009 X1009 2.9 2.5
# # NOTE: on needs edit the database first to ensure no illegal commas
# grep -wf <(cut -f1 200HIV_annotated.plink.fam; echo "id") ../../datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv \
#     | tr "," "\t" \
#     | cut -d$'\t' -f 1,3,5,7 \
#     | awk '{print $1"\t"$0}' \
#     > 200HIV_phenotypes.txt
# # Edit the file, change header
# 
# # Create a file including covariates
# grep -wf <(cut -f1 200HIV_annotated.plink.fam; echo "id") ../../datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv \
#     | tr "," "\t" \
#     | cut -d$'\t' -f 1,3,4,13,14 \
#     | awk '{print $1"\t"$0}' \
#     > 200HIV_covariates.txt
# # Edit the file, change header


# Plink: QTL mapping and permutation
# Because Plink2.0 implemented multiple threads for the majority of its
# funciton, using Plink2.0 could be a better choice.
module purge
module load plink/1.9-foss-2015b
echo "Modules used for Plink:"
module list

i_dir=/home/umcg-zzhang/Documents/projects/200HIV/inputs/dosages/Plink

plink --dosage ${i_dir}/200HIV_annotated.plink.dosage.gz format=1 \
    --fam ${i_dir}/200HIV_annotated.plink.fam \
    --map ${i_dir}/200HIV_annotated.plink.map \
    --pheno-merge \
    --pheno-name RNAvsDNA_CD4LOG \
    --pheno ${i_dir}/200HIV_phenotypes_RNAvsDNA.txt \
    --allow-no-sex \
    --covar ${i_dir}/200HIV_covariates.txt \
    --covar-name AGE, GENDER, HIV_DURATION, CD4_NADIR
