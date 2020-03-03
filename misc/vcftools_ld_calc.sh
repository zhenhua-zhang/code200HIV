#!/bin/bash
#SBATCH --mem=2G
#SBATCH --cpus=1
#SBATCH --time=0:30:0
#SBATCH --output=%j-%u-vcftools_ld_calc.log

source /apps/modules/modules.bashrc

# The postion of rs7113204 is chr11:500248

module purge
module load VCFtools
module list

vcf_iptf=../../annotated/chr11_annotated.vcf.gz
ld_optf=chr11_rs7113204_win50k
tm_optf=tmp

wincn=500248  # The postion of target SNP.
winsz=50000
winfm=$[ ${wincn} - ${winsz} / 2 ]
winto=$[ ${wincn} + ${winsz} / 2 ]

vcftools --gzvcf ${vcf_iptf} \
    --chr 11 \
    --recode \
    --from-bp ${winfm} \
    --to-bp ${winto} \
    --out ${tm_optf} \
    && vcftools \
    --hap-r2 \
    --vcf ${tm_optf}.recode.vcf \
    --out ${ld_optf}

rm -f ${tm_optf}*
