#!/usr/bin/env Rscript
# TODO: Commandline arguments to set the splitter.

library(ggplot2)
library(stringr)
library(reshape2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# args[1] genotype dosage file
# args[2] genotype information file
# args[3] target SNP
# args[4] phenotype file
# args[5] target phenotype
#

phenotype_file <- NULL
target_phenotype <- NULL
if (length(args) < 2) {
    stop("Wrong number of arguments, at least 2, but found ", length(args))
} else if (length(args) == 2) {
    target_snp <- "rs7817589"
} else if (length(args) == 3) {
    target_snp <- args[3]
} else if (length(args) == 4) {
    target_snp <- args[3]
    phenotype_file <- args[4]
} else if (length(args) == 5){
    target_snp <- args[3]
    phenotype_file <- args[4]
    target_phenotype <- args[5]
} else {
    warning("More than 4 arguments are given, only pick up first 4!\n")
}

genotype_dosage_file <- args[1]
genotype_info_file <- args[2]

genotype_dosage <- fread(genotype_dosage_file, data.table = FALSE)

genotype_col_names <- colnames(genotype_dosage)
target_snp_dosage <- genotype_dosage[
    which(genotype_dosage$id == target_snp),
    genotype_col_names[which(! genotype_col_names %in% c("id"))]
]

target_snp_round <- round(target_snp_dosage)

genotype_info <- fread(genotype_info_file, data.table = FALSE)
target_snp_info <- genotype_info[genotype_info$rsID == target_snp, ]

chrom <- target_snp_info["SequenceName"]
position <- target_snp_info["Position"]
effect_allele <- target_snp_info["EffectAllele"]
alternative_allele <- target_snp_info["AlternativeAllele"]

effect_genotype <- paste0(effect_allele, effect_allele)
heterozygous_genotype <- paste0(effect_allele, alternative_allele)
alternative_genotype <- paste0(alternative_allele, alternative_allele)

target_snp_genotype <- sapply(
    X = target_snp_round, FUN = function(e) {
        if (e == 0) {
            return(alternative_genotype)
        } else if (e == 1) {
            return(heterozygous_genotype)
        } else if (e == 2) {
            return(effect_genotype)
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

if (!is.null(phenotype_file)) {
    phenotype_level <- fread(phenotype_file, data.table = FALSE)
    phenotype_level <- reshape2::dcast(reshape::melt(phenotype_level, variable.name = "phenotype"), ...~id)
    
    rownames(phenotype_level) <- phenotype_level$variable
    phenotype_level <- phenotype_level[, which(!colnames(phenotype_level) %in% c("phenotype"))]
    
    phenotype_col_names <- colnames(phenotype_level)
    common_samples <- intersect(phenotype_col_names, genotype_col_names)
    
    if (is.null(target_phenotype)){
        target_phenotype <- "prcnt_mnct"
    }
    
    phenotype_level_chosen <- phenotype_level[target_phenotype, common_samples]
    target_snp_genotype_chosen <- target_snp_genotype[common_samples]
    
    phenotype_genotype_chosen <- rbind(phenotype = phenotype_level_chosen, genotype = target_snp_genotype_chosen)
    phenotype_genotype_chosen["id"] <- rownames(phenotype_genotype_chosen)
    
    phenotype_genotype_chosen <- reshape2::dcast(reshape2::melt(phenotype_genotype_chosen, id.vars="id"), ...~id)
    phenotype_genotype_chosen$genotype <- factor(
        phenotype_genotype_chosen$genotype,
        levels = sort(unique(phenotype_genotype_chosen$genotype), decreasing = (effect_allele > alternative_allele))
    )


    count_per_genotype <- table(phenotype_genotype_chosen$genotype)
    genotypes <- sort(names(count_per_genotype), decreasing = (effect_allele > alternative_allele))
    x_tick_lables <- paste0(paste(genotypes, count_per_genotype, sep="("), ")")
    
    ylabel <- str_glue("Level of {target_phenotype}")
    xlabel <- str_glue("SNP: {target_snp}; Position: chr{chrom}:{position}; Variant: {effect_allele}>{alternative_allele}")
    ftitle <- str_glue("Phenotype level per genotype")
    
    aov_test <- aov(phenotype~genotype, data=phenotype_genotype_chosen)
    aov_summary <- summary.aov(aov_test)
    print(aov_summary)

    if (aov_summary[[1]]["genotype", "Pr(>F)"] <= 0.9) {
        tukey_test <- TukeyHSD(aov_test, ordered = TRUE)
        print(tukey_test)
    }
    
    g <- ggplot(data = phenotype_genotype_chosen) + theme_bw()
    g <- g + geom_boxplot(aes(x=genotype, y=as.numeric(phenotype), color=genotype))
    g <- g + geom_point(aes(x=genotype, y=as.numeric(phenotype), color=genotype))
    g <- g + labs(title = ftitle, x = xlabel, y = ylabel)
    g <- g + scale_x_discrete(labels = x_tick_lables)

    
    plot_output_file <- str_glue("phenotypeLevelPerGenotype_{target_phenotype}_{target_snp}.pdf")
    ggsave(plot_output_file)
} else {
    warning("No phenotype is given, skipping draw phenotype level per genotype plot...")
}
