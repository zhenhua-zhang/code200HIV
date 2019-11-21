#!/usr/bin/env Rscript
# TODO: Commandline arguments to set the splitter.

library(data.table)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Wrong number of arguments, at least 2, but found ", length(args))
} else if (length(args) == 2) {
    target_snp <- "rs7817589"
} else if (length(args) == 3) {
    target_snp <- args[3]
} else {
    warning("More than 3 arguments are given, only pick up first three!\n")
}

genotype_dosage_file <- args[1]
genotype_info_file <- args[2]

genotype_dosage <- fread(genotype_dosage_file, sep = " ", data.table = FALSE)
col_names <- colnames(genotype_dosage)

target_snp_dosage <- genotype_dosage[
    which(genotype_dosage$id == target_snp),
    col_names[which(! col_names %in% c("id"))]
]

target_snp_round <- round(target_snp_dosage)

genotype_info <- fread(genotype_info_file, sep = " ", data.table = FALSE)
target_snp_info <- genotype_info[genotype_info$rsID == target_snp, ]

chrom <- target_snp_info["SequenceName"]
position <- target_snp_info["Position"]
effect_allele <- target_snp_info["EffectAllele"]
alternative_allele <- target_snp_info["AlternativeAllele"]

target_snp_genotype <- sapply(
    X = target_snp_round, FUN = function(e) {
        if (e == 0) {
            return(paste0(alternative_allele, alternative_allele))
        } else if (e == 1) {
            return(paste0(effect_allele, alternative_allele))
        } else if (e == 2) {
            return(paste0(effect_allele, effect_allele))
        } else {
            return(NULL)
        }
    }
)

snp_info <- str_glue(
    "Information SNP {target_snp}: \n  Position: chr{chrom}:{position}\n  Alleles: {effect_allele}>{alternative_allele}"
)

print(snp_info)
print(table(target_snp_genotype))