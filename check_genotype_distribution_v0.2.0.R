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

#NOTE:
#    1. The effect alleles are encoded as 2
#    2. If covariates are given, the phenotype level will be substitued by the
#    GLM residuals. The residuals will be scaled if --scale-phenotype is setted.
#    3. The genotype are the dosage value from the genotype array after
#    imputation, not the category genotypes.
#    4. If the flag --use-genotype-symbols is given, the regression will be done using
#    genotype symbols instead of genotype dosage. And consequently, there will
#    no genotype dodage in the work_dtfm. Will this be a problem?

library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)
library(stringr, warn.conflicts = FALSE, quietly = TRUE)
library(optparse, warn.conflicts = FALSE, quietly = TRUE)
library(data.table, warn.conflicts = FALSE, quietly = TRUE)

# A smarter read function
smtread <- function(file_path, ..., idxc = "id", kpr = NULL, rmr = NULL, kpc = NULL, rmc = NULL,
                    trps = FALSE) {
    # A function to read your file smarter
    # 1. Remove / choose columns and rows after loading the file.
    # 2. Transform the data.frame after loading the file.
    # 3. It does remove and choose first, and then do the transform, if all
    #    of them are requested.
    # 4. You also would like to give the id columns, which will help make life
    #    easier

    dtfm <- fread(
        file_path,
        data.table = FALSE, stringsAsFactors = FALSE, verbose = FALSE, ...
    )

    col_names <- colnames(dtfm)
    if (idxc %in% col_names) {
        rownames(dtfm) <- dtfm[, idxc]
        dtfm <- dtfm[, col_names[!col_names %in% c(idxc)]]
    } else {
        warning("The given `idxc = ", idxc, "` is not in the column names")
    }

    if (!is.null(kpr)) {
        kpt_rows <- kpr
    } else {
        kpt_rows <- rownames(dtfm)
    }

    if (!is.null(rmr)) {
        kpt_rows <- kpt_rows[!kpt_rows %in% rmr]
    }

    if (!is.null(kpc)) {
        kpt_cols <- kpc
    } else {
        kpt_cols <- colnames(dtfm)
    }

    if (!is.null(rmc)) {
        kpt_cols <- kpt_cols[!kpt_cols %in% rmc]
    }

    if (!all(rownames(dtfm) %in% kpt_rows)) {
        if (length(kpt_rows) > 1) {
            dtfm <- dtfm[kpt_rows, ]
        } else {
            singlten <- as.list(dtfm[kpt_rows, ])
            names(singlten) <- colnames(dtfm)
            dtfm <- data.frame(singlten, row.names = kpt_rows)
        }
    }

    if (!all(colnames(dtfm) %in% kpt_cols)) {
        if (length(kpt_cols) > 1) {
            dtfm <- dtfm[, kpt_cols]
        } else {
            singlten <- dtfm[, kpt_cols]
            dtfm <- data.frame(singlten, row.names = rownames(dtfm))
            colnames(dtfm) <- kpt_cols
        }
    }

    if (trps) {
        # Any is character
        anic <- any(sapply(dtfm, function(e) {
            return(typeof(e) == "character")
        }))

        if (anic) {
            stop("There's character in the data.frame, plase remove them then try transform again")
        }
        dtfm <- as.data.frame(t(dtfm))
    }

    return(dtfm)
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
    help = "Whether use genotype symbols in the regression instead of genotype dosage."
)

parser <- add_option(
    parser, c("-d", "--genotype-dosage-file"),
    action = "store", dest = "genotype_dosage_file", type = "character",
    help = "The genotype dosage file (could be compressed)."
)

parser <- add_option(
    parser, c("--genotype-dosage-idx-col"),
    action = "store", dest = "genotype_dosage_idx_col", type = "character", default = "id",
    help = "The id column in genotype file, usually its the name of column of SNP id. Default: id"
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
    help = "The columns will be used. Default: rsID,SequenceName,Position,EffectAllele,AlternativeAllele"
)

parser <- add_option(
    parser, c("-t", "--target-snp"),
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
    help = "The id column in phenotype file"
)

parser <- add_option(
    parser, c("-c", "--covariate-file"),
    action = "store", dest = "covariate_file", type = "character",
    help = "The file including covariates"
)

parser <- add_option(
    parser, c("-C", "--target-covariates"),
    action = "store", dest = "target_covariates", type = "character",
    help = "Target covariates"
)

parser <- add_option(
    parser, c("--covariate-idx-col"),
    action = "store", dest = "covariate_idx_col", type = "character", default = "id",
    help = "The id column in covariates file"
)

parser <- add_option(
    parser, c("-O", "--output-dir"),
    action = "store", dest = "output_dir", type = "character", default = "output_dir",
    help = "The output direcotry which will be created if not exists. Default: output_dir"
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

target_snp <- opts$target_snp
if (is.null(target_snp)) {
    print_help(parser)
    stop("-t/--target-snp is required!")
}

output_dir <- opts$output_dir
if (! dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

genotype_dosage_idx_col <- opts$genotype_dosage_idx_col
genotype_dosage <- smtread(genotype_dosage_file, idxc = genotype_dosage_idx_col)
genotype_sample_names <- colnames(genotype_dosage)

target_snp_dosage <- genotype_dosage[target_snp, ]
target_snp_round <- round(target_snp_dosage)

#<< Add genotype information for target SNP
genotype_info_cols <- opts$genotype_info_cols
genotype_info_cols_vec <- str_split(genotype_info_cols, pattern = ",")[[1]]

if (length(genotype_info_cols_vec) < 5) {
    stop(
        "The length of --genotype-info-cols should be 5 and splitted by comma.\n",
        "    e.g: rsID,SequenceName,Position,EffectAllele,AlternativeAllele"
    )
}

genotype_info <- smtread(genotype_info_file, idxc = genotype_info_cols_vec[1])
target_snp_info <- genotype_info[target_snp, ]

chrom <- target_snp_info[genotype_info_cols_vec[2]]
position <- target_snp_info[genotype_info_cols_vec[3]]
effect_allele <- target_snp_info[genotype_info_cols_vec[4]]
alternative_allele <- target_snp_info[genotype_info_cols_vec[5]]

effect_genotype <- paste0(effect_allele, effect_allele)
heterozygous_genotype <- paste0(effect_allele, alternative_allele)
alternative_genotype <- paste0(alternative_allele, alternative_allele)

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
# Add genotype information for target SNP >>

count_per_genotype <- table(target_snp_genotype)
genotypes <- sort(names(count_per_genotype), decreasing = (effect_allele > alternative_allele))
genotype_freq <- paste0(paste(genotypes, count_per_genotype[genotypes], sep = "("), ")")
snp_info <- str_glue("SNP: {target_snp};", " Chr: {chrom};", " Pos: {position};", " Var: {effect_allele}>{alternative_allele}")
cat("<<<\n", snp_info, "; Genotype freq: ", paste(genotype_freq), "\n>>>\n", sep = "")

#
## Phenotype level per genotype plot, with correction for covariates
#
phenotype_file <- opts$phenotype_file
phenotype_idx_col <- opts$phenotype_idx_col
target_phenotypes <- opts$target_phenotypes

covariate_file <- opts$covariate_file
covariate_idx_col <- opts$covariate_idx_col
target_covariates <- opts$target_covariates

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

    # Covariates, if given
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
    for (target_phenotype in target_phenotypes_vec) {
        cat("<<< Summary for generalized linear regression for phenotype:", target_phenotype, "\n")

        if (use_smb) {
            target_snp_genotype_chosen <- target_snp_genotype[common_samples]
            my_level <- sort(unique(target_snp_genotype_chosen), decreasing = (effect_allele > alternative_allele))
            target_snp_genotype_chosen <- factor(target_snp_genotype_chosen, levels = my_level)
        } else {
            target_snp_genotype_chosen <-  target_snp_dosage[common_samples]
        }

        phenotype_level_chosen <- phenotype_level[common_samples, target_phenotype]
        work_dtfm <- cbind(pntp = unlist(phenotype_level_chosen), gntp = unlist(target_snp_genotype_chosen))

        if (! is.null(covariate_level_chosen)) {
            work_dtfm <- as.data.frame(cbind(work_dtfm, covariate_level_chosen))
        } else {
            warning("[WARN] No covariates is given, therefore skipping adjustment ...")
            work_dtfm <- as.data.frame(work_dtfm)
        }


        glm_fit <- glm(pntp ~ ., data = work_dtfm)

        glm_fit_sum <- summary(glm_fit)
        print(glm_fit_sum)
        cat(">>> Summary for generalized linear regression for phenotype:", target_phenotype, "\n\n\n")

        work_dtfm["gntp_enc"] <- target_snp_genotype[rownames(work_dtfm)]
        my_level <- sort(unique(work_dtfm[, "gntp_enc"]), decreasing = (effect_allele > alternative_allele))
        work_dtfm["gntp_enc"] <- factor(work_dtfm[, "gntp_enc"], levels = my_level)

        if (scale_phenotype) {
            work_dtfm["pntp_adj"] <- scale(glm_fit$residuals[rownames(work_dtfm)])
        } else {
            work_dtfm["pntp_adj"] <- glm_fit$residuals[rownames(work_dtfm)]
        }

        # Title, x label, x ticklabel, y label
        count_per_genotype <- table(work_dtfm["gntp_enc"]) # The genotype distribution of samples chosen for regression
        genotypes <- sort(names(count_per_genotype), decreasing = (effect_allele > alternative_allele))

        gntp_pval_beta <- coef(glm_fit_sum)["gntp", c(1, 4)]
        gntp_beta <- gntp_pval_beta[1]
        gntp_pval <- gntp_pval_beta[2]

        ftitle <- str_glue("Phenotype level per genotype")
        x_tick_lables <- paste0(paste(genotypes, count_per_genotype, sep = "("), ")")
        xlabel <- paste0(snp_info, "; p-val: ", signif(gntp_pval, 4), "; beta: ", signif(gntp_beta, 4))
        ylabel <- ifelse(
            is.null(target_covariates),
            str_glue("Level of {target_phenotype}"),
            str_glue("Level of {target_phenotype} adjusted by {target_covariates}")
        )

        g <- ggplot(data = work_dtfm) + theme_bw()
        g <- g + geom_boxplot(aes(x = gntp_enc, y = pntp_adj, color = gntp_enc))
        g <- g + geom_point(aes(x = gntp_enc, y = pntp_adj), size = 0.5)
        g <- g + labs(title = ftitle, x = xlabel, y = ylabel, color = "Genotype")
        g <- g + scale_x_discrete(labels = x_tick_lables)
        g <- g + theme(plot.title = element_text(hjust = 0.5))

        opt_path <- str_glue("{output_dir}/phenotypeLevelPerGenotype_{target_phenotype}_{target_snp}.pdf")
        ggsave(opt_path, width = 10, height = 10)
        cat("\n----------------------------------------------------------------------------\n")
    }
} else {
    warning("No phenotype is given, skipping draw phenotype level per genotype plot...")
}
