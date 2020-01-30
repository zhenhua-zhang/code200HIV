#!/usr/bin/env Rscript

# Create date: 2019-Nov-25
# Last update: 2019-Dec-4
# Version    : 0.2.0
# Author     : Zhenhua Zhang
# Email      : zhenhua.zhang217@gmail.com

# A script to draw the phenotype level per genotype in 200HIV project. It could
# be used in out of the project, for more information please contact the author.

#TODO:
#    1. Cope with multiple target SNP.
#    2. Covariates should be optional

#NOTE:
#    1. The effect alleles are encoded as 2
#    2. If covariates are given, the phenotype level will be substitued by the
#    GLM residuals. The residuals will be scaled if --scale-phenotype is setted.
#    3. The genotype are the dosage value from the genotype array after
#    imputation, not the category genotypes.
#    4. If the flag --use-genotype-symbols is given, the regression will be done using
#    genotype symbols instead of genotype dosage. And consequently, there will
#    no genotype dodage in the work_dtfm. Will this be a problem? According to
#    MatrixEQTL, it tests the association between each SNP and each transcript
#    by modelling the effect of genotype as either additive linear (least
#    squares model) or categorical (ANOVA model). Therefore, using genotype
#    dosage in the GLM model should result consist outcomes as MatrixEQTL does.
#    5. Perhaps a generalized additive model (GAM) helps? 为什么 QTL mapping
#    不用更复杂的模型，比如GLM，GAM等？ additive model, least squares model?

library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)
library(stringr, warn.conflicts = FALSE, quietly = TRUE)
library(optparse, warn.conflicts = FALSE, quietly = TRUE)
library(data.table, warn.conflicts = FALSE, quietly = TRUE)

if (file.exists("./utils.R")) {
    source("./utils.R")
}

# << Options and flags
parser <- OptionParser(
    description = "Draw phenotype level per genotype plot for given SNP"
)

parser <- add_option(
    parser, c("--scale-phenotype"),
    action = "store_true", dest = "scale_phenotype",
    help = "Whether scale phenotype for the phenotype level per genotype box plots."
)

parser <- add_option(
    parser, c("--use-genotype-symbols"),
    action = "store_true", dest = "use_smb",
    help = paste(
        "Whether use genotype symbols in the regression instead of genotype dosage.",
        "If the flag is on, the genotype will be taken as categorical variables.",
        "Otherwise, it's additive linear."
    )
)

parser <- add_option(
    parser, c("-d", "--genotype-dosage-file"),
    action = "store", dest = "genotype_dosage_file", type = "character",
    help = "The genotype dosage file (could be compressed)."
)

parser <- add_option(
    parser, c("--genotype-dosage-idx-col"),
    action = "store", dest = "genotype_dosage_idx_col", type = "character", default = "id",
    help = "The id column in genotype file, usually its the name of column of SNP id. Default: %default"
)

parser <- add_option(
    parser, c("-i", "--genotype-info-file"),
    action = "store", dest = "genotype_info_file", type = "character",
    help = "The genotype information file (could be compressed)"
)

parser <- add_option(
    parser, c("--genotype-info-cols"),
    action = "store", dest = "genotype_info_cols", type = "character",
    default = "rsID,SequenceName,Position,EffectAllele,AlternativeAllele",
    help = "The columns will be used. Default: %default"
)

parser <- add_option(
    parser, c("--target-snp"),
    action = "store", dest = "target_snp", type = "character",
    help = "The SNP which the phenotype level per genotype plot will be plotted."
)

parser <- add_option(
    parser, c("-p", "--phenotype-file"),
    action = "store", dest = "phenotype_file", type = "character",
    help = "The phenotype file"
)

parser <- add_option(
    parser, c("-P", "--target-phenotypes"),
    action = "store", dest = "target_phenotypes", type = "character",
    help = "The phenotype will be used in the plots"
)

parser <- add_option(
    parser, c("--phenotype-idx-col"),
    action = "store", dest = "phenotype_idx_col", type = "character", default = "id",
    help = "The id column in phenotype file. Default: %default"
)

parser <- add_option(
    parser, c("-c", "--covariate-file"),
    action = "store", dest = "covariate_file", type = "character",
    help = "The file including covariates"
)

parser <- add_option(
    parser, c("--target-covariates"),
    action = "store", dest = "target_covariates", type = "character",
    help = "Target covariates"
)

parser <- add_option(
    parser, c("--covariate-idx-col"),
    action = "store", dest = "covariate_idx_col", type = "character", default = "id",
    help = "The id column in covariates file. Default: %default"
)

parser <- add_option(
    parser, c("-O", "--output-dir"),
    action = "store", dest = "output_dir", type = "character", default = "output_dir",
    help = "The output direcotry which will be created if not exists. Default: %default"
)


opts_args <- parse_args2(parser)
opts <- opts_args$options
args <- opts_args$args
# Options and flags >>

genotype_dosage_file <- opts$genotype_dosage_file
if (is.null(genotype_dosage_file)) {
    print_help(parser)
    stop("-d/--genotype-dosage-file is required!")
}

genotype_info_file <- opts$genotype_info_file
if (is.null(genotype_info_file)) {
    print_help(parser)
    stop("-i/--genotype-info-file is required!")
}

genotype_info_cols <- opts$genotype_info_cols
genotype_info_cols_vec <- str_split(genotype_info_cols, pattern = ",")[[1]]
if (length(genotype_info_cols_vec) < 5) {
    print_help(parser)
    stop(
        "The length of --genotype-info-cols should be 5 and splitted by comma.\n",
        "    e.g: rsID,SequenceName,Position,EffectAllele,AlternativeAllele"
    )
}

target_snp <- opts$target_snp
if (is.null(target_snp)) {
    print_help(parser)
    stop("-t/--target-snp is required!")
} else if (grepl(",", target_snp)) {
    target_snp_vec <- str_split(target_snp, ",")[[1]]
} else {
    target_snp_vec <- c(target_snp)
}

output_dir <- opts$output_dir
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
} else {
    warning("The given output direcotry exists, will using it directly!")
}

#
## Genotypes
#
genotype_dosage_idx_col <- opts$genotype_dosage_idx_col
genotype_dosage <- smtread(genotype_dosage_file, idxc = genotype_dosage_idx_col)
genotype_sample_names <- colnames(genotype_dosage)
genotype_info <- smtread(genotype_info_file, idxc = genotype_info_cols_vec[1])

#
## Phenotypes
#
phenotype_file <- opts$phenotype_file
phenotype_idx_col <- opts$phenotype_idx_col
target_phenotypes <- opts$target_phenotypes

common_samples <- NULL
if (! is.null(phenotype_file)) {
    # Phenotype level
    if (is.null(target_phenotypes)) {
        phenotype_level <- smtread(phenotype_file, idxc = phenotype_idx_col)
        target_phenotypes_vec <- colnames(phenotype_level)
    } else {
        target_phenotypes_vec <- str_split(target_phenotypes, pattern = ",")[[1]]
        phenotype_level <- smtread(phenotype_file, idxc = phenotype_idx_col, kpc = target_phenotypes_vec)
    }
    common_samples <- intersect(rownames(phenotype_level), genotype_sample_names)
}

#
## Covariates
#
covariate_file <- opts$covariate_file
covariate_idx_col <- opts$covariate_idx_col
target_covariates <- opts$target_covariates

if (! is.null(covariate_file)) {
    if (is.null(target_covariates)) {
        covariate_level <- smtread(covariate_file, idxc = covariate_idx_col)
        target_covariates_vec <- colnames(covariate_level)
    } else {
        target_covariates_vec <- str_split(target_covariates, pattern = ",")[[1]]
        covariate_level <- smtread(covariate_file, idxc = covariate_idx_col, kpc = target_covariates_vec)
    }

    common_samples <- intersect(common_samples, rownames(covariate_level))

    if (length(target_covariates_vec) > 1) {
        covariate_level_chosen <- covariate_level[common_samples, ]
    } else {
        covariate_level_chosen <- unlist(covariate_level[common_samples, ])
    }
} else {
    covariate_level_chosen <- NULL
}

use_smb <- ifelse(is.null(opts$use_smb), FALSE, opts$use_smb)
scale_phenotype <- ifelse(is.null(opts$scale_phenotype), FALSE, opts$scale_phenotype)

for (target_snp in target_snp_vec) {
    output_file <- str_glue("{output_dir}/check_gntp_pntp_{target_snp}.txt" )

    capture.output(cat(str_glue("Number of samples in common: {length(common_samples)}")), file = output_file, append = TRUE)

    target_snp_dosage <- genotype_dosage[target_snp, ]
    target_snp_round <- round(target_snp_dosage)

    target_snp_info <- genotype_info[target_snp, ]

    chrom <- target_snp_info[genotype_info_cols_vec[2]]
    position <- target_snp_info[genotype_info_cols_vec[3]]
    effect_allele <- target_snp_info[genotype_info_cols_vec[4]]
    alternative_allele <- target_snp_info[genotype_info_cols_vec[5]]

    effect_genotype <- paste0(effect_allele, effect_allele)
    heterozygous_genotype <- paste0(effect_allele, alternative_allele)
    alternative_genotype <- paste0(alternative_allele, alternative_allele)

    genotype_code_vec <- c(effect_genotype, heterozygous_genotype, alternative_genotype)
    genotype_code_num <- length(genotype_code_vec)

    target_snp_genotype <- sapply(
      X = target_snp_round,
      FUN = function(e) {
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

    count_per_genotype <- table(target_snp_genotype)
    genotypes <- sort(names(count_per_genotype), decreasing = (effect_allele > alternative_allele))
    genotype_freq <- paste0(paste(genotypes, count_per_genotype[genotypes], sep = "("), ")")
    snp_info <- str_glue("SNP: {target_snp};", " Chr: {chrom};", " Pos: {position};", " Var: {effect_allele}>{alternative_allele}")
    capture.output(cat("<<<\n", snp_info, "; Genotype freq: ", paste(genotype_freq), "\n>>>\n", sep = ""), file = output_file, append = TRUE)

    for (target_phenotype in target_phenotypes_vec) {
        capture.output(cat("<<< Summary for generalized linear regression for phenotype:", target_phenotype, "\n"), file = output_file, append = TRUE)

        genotype_level_chosen <-  target_snp_dosage[common_samples]
        phenotype_level_chosen <- phenotype_level[common_samples, target_phenotype]
        work_dtfm <- cbind(pntp = unlist(phenotype_level_chosen), gntp = unlist(genotype_level_chosen))

        if (! is.null(covariate_level_chosen)) {
            work_dtfm <- as.data.frame(cbind(work_dtfm, covariate_level_chosen))
        } else {
            warning("[WARN] No covariates is given, therefore skipping adjustment ...")
            work_dtfm <- as.data.frame(work_dtfm)
        }

        work_dtfm["gntp_enc"] <- factor(target_snp_genotype[rownames(work_dtfm)])

        if (use_smb) {
            compare_table <- matrix(NA, nrow = genotype_code_num, ncol = genotype_code_num)
            rowname_vec <- paste0("gntp_enc", genotype_code_vec)
            rownames(compare_table) <- colnames(compare_table) <- rowname_vec

            for (genotype_code in genotype_code_vec) {
                work_dtfm[, "gntp_enc"] <- relevel(work_dtfm[, "gntp_enc"], ref = genotype_code)
                glm_fit <- glm(pntp ~ . - gntp, data = work_dtfm)

                glm_fit_sum <- summary(glm_fit)
                glm_fit_coef <- glm_fit_sum$coefficient
                baseline_rowname <- paste0("gntp_enc", genotype_code)
                rest_rowname_vec <- rowname_vec[! rowname_vec %in% c(baseline_rowname)]
                compare_table[rest_rowname_vec, baseline_rowname] <- glm_fit_coef[rest_rowname_vec, "Pr(>|t|)"]
            }

            rownames(compare_table) <- colnames(compare_table) <- genotype_code_vec
            capture.output(print(compare_table), file = output_file, append = TRUE)

            gntp_beta <- gntp_pval <- "NULL"
        } else {
            glm_fit <- glm(pntp ~ . - gntp_enc, data = work_dtfm)

            glm_fit_sum <- summary(glm_fit)
            gntp_pval_beta <- coef(glm_fit_sum)["gntp", c(1, 4)]
            gntp_beta <- gntp_pval_beta[1]
            gntp_pval <- gntp_pval_beta[2]
        }

        capture.output(print(glm_fit_sum), file = output_file, append = TRUE)
        capture.output(cat(">>> Summary for generalized linear regression for phenotype:", target_phenotype, "\n\n\n"), file = output_file, append = TRUE)

        if (scale_phenotype) {
            work_dtfm["pntp_adj"] <- scale(glm_fit$residuals[rownames(work_dtfm)])
        } else {
            work_dtfm["pntp_adj"] <- glm_fit$residuals[rownames(work_dtfm)]
        }

        work_dtfm["gntp_enc"] <- factor(work_dtfm[, "gntp_enc"], levels = genotypes)

        # Title, x label, x ticklabel, y label
        gntp_pval <- ifelse(is.character(gntp_pval), gntp_pval, signif(gntp_pval, 4))
        gntp_beta <- ifelse(is.character(gntp_beta), gntp_beta, signif(gntp_beta, 4))

        xlabel <- paste0(snp_info, "; p-val: ", gntp_pval, "; beta: ", gntp_beta, 4)

        count_per_genotype <- table(work_dtfm[, "gntp_enc"])
        x_tick_lables <- paste0(paste(genotypes, count_per_genotype[genotypes], sep = "("), ")")

        ylabel <- ifelse(
             is.null(target_covariates),
             str_glue("Level of {target_phenotype}"),
             str_glue("Level of ({target_phenotype} | {target_covariates})")
        )
        ftitle <- str_glue("Boxplot for {target_snp} x {target_phenotype}")

        g <- ggplot(data = work_dtfm) + theme_bw()
        g <- g + geom_boxplot(aes(x = gntp_enc, y = pntp_adj, color = gntp_enc))
        g <- g + geom_point(aes(x = gntp_enc, y = pntp_adj), size = 0.5)
        g <- g + labs(title = ftitle, x = xlabel, y = ylabel, color = "Genotype")
        g <- g + scale_x_discrete(labels = x_tick_lables)
        g <- g + theme(plot.title = element_text(hjust = 0.5))

        opt_path <- str_glue("{output_dir}/boxplot_{target_phenotype}_{target_snp}.pdf")
        ggsave(opt_path, g, width = 10, height = 10)
        capture.output(cat("\n----------------------------------------------------------------------------\n"), file = output_file, append = TRUE)
    }
}
