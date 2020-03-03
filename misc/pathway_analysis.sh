#!/bin/bash
# Do the pathway analysis

# Fetch genes in up- and down- stream
module load BEDTools
bedtools intersect \
    -a <(echo -e "11\t247\t1000248") \
    -b /apps/data/ftp.ensembl.org/pub/release-75/gff/Homo_sapiens.GRCh37.75.gff \
    -wb \
    | sed -n '/gene/ s/.*;Name=\(.*\);.*/\1/p' \
    > genes_around_rs7113204.txt

# Fetch eQTL genes based on SNPs around rs7113204
# Some gene have been remove because of duplication from different transcript ID
zgrep -wf ../significant_QTLs/candidate_chr11_snpids.tsv \
    ../../../../inputs/eQTLGen/cis-eQTL_significant_20181017.txt \
    | cut -f9 \
    | sort \
    | uniq \
    > eQTL_genes_by_SNPs_around_rs7113204.txt

