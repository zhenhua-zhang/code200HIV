#!/usr/bin/env Rscript

# Create date: 2019-Nov-25
# Last update: Mon 30 Dec 2019 11:18:41 AM CET
# Version    : 0.2.0
# License    : MIT
# Author     : Zhenhua Zhang
# Email      : zhenhua.zhang217@gmail.com

# NOTE:
#    1. The id column is mandatory for each input file, i.e phenotypes, covariates
#    2. For the phenotype and the covariate file, the script supposes the
#    columns are traits, while each row represents one sample. Therefore, if the
#    --trps-pntp-dtfm or --trps-cvrt-dtfm is given, it means the input file
#    isn't in the formated as the script supposes.
#    3. Dosage = Pr(Het|Data) + 2*Pr(Alt|Data)
#    4. The script was developed under "R/3.5.1", under other version of R it
#    should also function, but not been tested.

# TODO:
#    1. A README.md to descript this shit.


library(qqman)
library(stringr)
library(GGally)
library(ggplot2)
library(optparse)
library(data.table)
library(MatrixEQTL)

if (file.exists("./utils.R")) {
    source("./utils.R")
}

debug <- FALSE

#
## Parsing CLI arguments
#
parser <- OptionParser(description = "A QTL mapping script based on MatrixEQTL")

parser <- add_option(
    parser, c("-w", "--work-dir"),
    action = "store", dest = "work_dir", type = "character", default = "./qtlmapping_opd",
    help = "Output direcotry."
)

parser <- add_option(
    parser, c("--run-flag"),
    action = "store", dest = "run_flag", type = "character", default = "test",
    help = "Running flag which help to discrimnate different runs."
)

# Phenotypes related parameters
parser <- add_option(
    parser, c("-p", "--pntp-file"),
    action = "store", dest = "pntp_file", type = "character",
    help = "Phenotype input file."
)

parser <- add_option(
    parser, c("--padding"),
    action = "store", dest = "padding", type = "integer", default = 4,
    help = "The times of standard deviation as the boundary of outliers."
)

parser <- add_option(
    parser, c("--target-pntp"),
    action = "store", dest = "target_pntp", type = "character",
    help = "The phenotypes will be used, if more than one, using comma as delimiter."
)

parser <- add_option(
    parser, c("--trps-pntp-dtfm"),
    action = "store_true", dest = "trps_pntp_dtfm",
    help = "Whether should transpose the data.frame of phenotypes to make the columns as sample ID."
)

parser <- add_option(
    parser, c("--pntp-idx-col"),
    action = "store", dest = "pntp_idx_col", type = "character", default = "id",
    help = "The id column in phenotype file"
)

# Covariates related parameters
parser <- add_option(
    parser, c("-c", "--cvrt-file"),
    action = "store", dest = "cvrt_file", type = "character",
    help = "Covariates file."
)

parser <- add_option(
    parser, c("--target-cvrt"),
    action = "store", dest = "target_cvrt", type = "character",
    help = "Covariates will be used, if more than one, using comma as delimiter."
)

parser <- add_option(
    parser, c("--trps-cvrt-dtfm"),
    action = "store_true", dest = "trps_cvrt_dtfm",
    help = "Whether should transpose the data.frame of covariates."
)

parser <- add_option(
    parser, c("--cvrt-idx-col"),
    action = "store", dest = "cvrt_idx_col", type = "character", default = "id",
    help = "The id column in covariates file"
)

# Phenotypes and covariates correlation
parser <- add_option(
    parser, c("-t", "--pwcor-trait"),
    action = "store", dest = "pwcor_trait", type = "character",
    help = "Traits will be correlated with in paire-wised way. The options will be ignored if --draw-pwcor isn't given."
)

# Genotypes related parameters
parser <- add_option(
    parser, c("-d", "--gntp-dosage-file"),
    action = "store", dest = "gntp_dosage_file", type = "character",
    help = "The genotype dosage file (could be compressed)."
)

parser <- add_option(
    parser, c("--genotype-dosage-idx-col"),
    action = "store", dest = "gntp_dosage_idx_col", type = "character", default = "id",
    help = "The id column in genotype file, usually its the name of column of SNP id. Default: id"
)

parser <- add_option(
    parser, c("-i", "--gntp-info-file"),
    action = "store", dest = "gntp_info_file", type = "character",
    help = "The genotype information file (could be compressed)."
)

parser <- add_option(
    parser, c("--gntp-info-cols"),
    action = "store", dest = "gntp_info_cols", type = "character",
    default = "rsID,SequenceName,Position,EffectAllele,AlternativeAllele",
    help = "The columns will be used. Default: rsID,SequenceName,Position,EffectAllele,AlternativeAllele"
)

parser <- add_option(
    parser, c("--maf-thrd"),
    action = "store", dest = "maf_thrd", type = "double", default = 0.05,
    help = "Minor allele frequency."
)

# Permutations
parser <- add_option(
    parser, c("--pm-times"),
    action = "store", dest = "pm_times", type = "integer", default = 0,
    help = paste(
        "How many times of permutations should be done. If it's less than 1,",
        "no permutation will be performed but only 'raw' data will be used in",
        "the mapping. Default: 0"
    )
)

parser <- add_option(
    parser, c("--pm-seed"),
    action = "store", dest = "pm_seed", type = "integer", default = 31415,
    help = "The random seed for permutation. Defautl: 31415"
)

# Misc
parser <- add_option(
    parser, c("--mhtn-fig-p-thrd"),
    action = "store", dest = "mhtn_fig_p_thrd", type = "double", default = 0.05,
    help = "The threshold of p-value for Manhattan plot. Default: 0.05"
)

parser <- add_option(
    parser, c("--draw-pwcor"),
    action = "store_true", dest = "draw_pwcor",
    help = "If the flag is given, the script will draw the pair-wise correlation plot for genotypes and covariates."
)

opts_args <- parse_args2(parser)
opts <- opts_args$options
args <- opts_args$args

#
## Visualization of the correlation of phenotypes and covariates
#
pntp_file <- opts$pntp_file
if (is.null(pntp_file)) {
    print_help(parser)
    stop("-p/--pntp-file is mandatory!")
}

gntp_dosage_file <- opts$gntp_dosage_file
if (is.null(gntp_dosage_file)) {
    print_help(parser)
    stop("-d/--gntp-dosage-file is mandatory!")
}

gntp_info_file <- opts$gntp_info_file
if (is.null(gntp_info_file)) {
    print_help(parser)
    stop("-i/--gntp-info-file is mandatory!")
}

gntp_info_cols <- opts$gntp_info_cols
gntp_info_cols_vec <- str_split(gntp_info_cols, pattern = ",")[[1]]
if (length(gntp_info_cols_vec) != 5) {
    print_help(parser)
    stop(
        "The length of --genotype-info-cols should be 5 and splitted by comma.\n",
        "    e.g: rsID,SequenceName,Position,EffectAllele,AlternativeAllele"
    )
}

run_flag <- opts$run_flag
work_dir <- opts$work_dir
if (! dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE)
} else {
    warning("The given work direcotry exists, will using it directly!")
}
setwd(work_dir)

target_pntp <- opts$target_pntp
if (is.null(target_pntp)) {
    target_pntp_vec <- NULL
} else {
    target_pntp_vec <- str_split(target_pntp, ",")[[1]]
}

pntp_idx_col <- opts$pntp_idx_col

trps_pntp_dtfm <- ifelse(is.null(opts$trps_pntp_dtfm), FALSE, opts$trps_pntp_dtfm)
if (trps_pntp_dtfm) {
    pntp_chosen <- smtread(pntp_file, idxc = pntp_idx_col, kpc = target_pntp_vec)
} else {
    # If not transpose the data.frame, meaning the columns are samples id,
    # therefore, pntp_idx_col should be the column refering the names of
    # phenotypes. The similar rule fits covariates input files.
    pntp_chosen <- smtread(pntp_file, idxc = pntp_idx_col, kpr = target_pntp_vec, trps = TRUE)
}

# Read covariates.
cvrt_file <- opts$cvrt_file
if (is.null(cvrt_file)) {
    with_cvrt <- FALSE
    warning("Not covariates is given, there won't be covariates in the model")
} else {
    with_cvrt <- TRUE

    target_cvrt <- opts$target_cvrt
    if (is.null(target_cvrt)) {
        target_cvrt_vec <- NULL
    } else {
        target_cvrt_vec <- str_split(target_cvrt, ",")[[1]]
    }

    cvrt_idx_col <- opts$cvrt_idx_col

    trps_cvrt_dtfm <- ifelse(is.null(opts$trps_cvrt_dtfm), FALSE, opts$trps_cvrt_dtfm)
    if (trps_cvrt_dtfm) {
        cvrt_chosen <- smtread(cvrt_file, idxc = cvrt_idx_col, kpc = target_cvrt_vec)
    } else {
        cvrt_chosen <- smtread(cvrt_file, idxc = cvrt_idx_col, kpr = target_cvrt_vec, trps = TRUE)
    }
}

# Whether do pair-wise correlation analysis for phenotypes and covariates (if
# supplied).
draw_pwcor <- ifelse(is.null(opts$draw_pwcor), FALSE, opts$draw_pwcor)
if (draw_pwcor) {
    if (with_cvrt) {
        pntp_cvrt_chosen <- merge(pntp_chosen, cvrt_chosen, by = "row.names", all = TRUE)
        rownames(pntp_cvrt_chosen) <- pntp_cvrt_chosen[, "Row.names"]
        pntp_cvrt_chosen <- pntp_cvrt_chosen[, -1]
    } else {
        pntp_cvrt_chosen <- pntp_chosen
    }

    pwcor_trait <- opts$pwcor_trait
    if (is.null(pwcor_trait)) {
        pwcor_trait_vec <- colnames(pntp_cvrt_chosen)
    } else {
        pwcor_trait_vec <- str_split(pwcor_trait, ",")[[1]]
    }

    if (dim(pntp_cvrt_chosen)[[2]] < 16) {
        pwcor_plot <- ggpairs(pntp_cvrt_chosen[, pwcor_trait_vec], cardinality_threshold = 16)
        ggsave(str_glue("phenotypes_covariates_pairwise_{run_flag}.pdf"), plot = pwcor_plot, width = 25, height = 25)
    } else {
        warning("More than 16 variables in the data.frame, exceeding cardinality_threshold")
    }
} else {
    warning("Skipping the pair-wise correlation analysis for genotypes and covariates.")
}

# TODO: add a section to draw heatmap of Spearman's rank correlation matrix???

# Remove outliers of phenotypes
# TODO: A CLI options to decide Whether to remove outliers or not.
padding <- opts$padding
pntp_means <- sapply(as.data.frame(pntp_chosen[, target_pntp_vec]), FUN = mean, na.rm = T)
pntp_sd <- sapply(as.data.frame(pntp_chosen[, target_pntp_vec]), FUN = sd, na.rm = T)

upper_bound <- pntp_means + padding * pntp_sd
lower_bound <- pntp_means - padding * pntp_sd

names(upper_bound) <- target_pntp_vec
names(lower_bound) <- target_pntp_vec

pntp_chosen_rownames <- rownames(pntp_chosen)
for (pntp in target_pntp_vec) {
    outlier_row_idx <- which((lower_bound[pntp] > pntp_chosen[, pntp]) | (pntp_chosen[, pntp] > upper_bound[pntp]))
    pntp_chosen[outlier_row_idx, pntp] <- NA
    cat("For ", pntp, " the outlier index out of ", padding, " * SD boundary [", lower_bound[pntp], ", ", upper_bound[pntp], "]: \n ", sep = "")
    cat(pntp_chosen_rownames[outlier_row_idx], "\n")
}

## Phenotypes
#
pntp_4me <- as.data.frame(t(pntp_chosen))
rownames(pntp_4me) <- colnames(pntp_chosen)
colnames(pntp_4me) <- rownames(pntp_chosen)
cat("Dim of phenotype data.frame will be piped into MatrixEQTL:", dim(pntp_4me), "\n")

## Covariates
#
if (with_cvrt) {
    cvrt_4me <- as.data.frame(t(cvrt_chosen))
    rownames(cvrt_4me) <- colnames(cvrt_chosen)
    colnames(cvrt_4me) <- rownames(cvrt_chosen)
    cat("Dim of covariate data.frame will be piped into MatrixEQTL:", dim(cvrt_4me), "\n")
} else {
    cat("[INFO] No covariate is supplied. Skipping covariates")
}

# Preprocessing genotypes
gntp_dosage_idx_col <- opts$gntp_dosage_idx_col
gntp_dosage <- smtread(gntp_dosage_file, idxc = gntp_dosage_idx_col)

maf_thrd <- opts$maf_thrd
gntp_round <- round(gntp_dosage) # For the phenotype level per genotype plot
gntp_af <- apply(gntp_round, 1, sum) / (2 * dim(gntp_dosage)[[2]])
nmvar_idx <- which((gntp_af >= maf_thrd) & (gntp_af <= 1 - maf_thrd))

#
## Geontypes
#
gntp_4me <- gntp_dosage[nmvar_idx, ]
cat("Dim of genotype data.frame will be piped into MatrixEQTL:", dim(gntp_4me), "\n")

pntp_samples <- colnames(pntp_4me)
gntp_samples <- colnames(gntp_4me)

if (with_cvrt) {
    cvrt_samples <- colnames(cvrt_4me)
    shared_samples <- intersect(pntp_samples, intersect(cvrt_samples, gntp_samples))

    # Write the covariates which are used in the QTL mapping into disk.
    cvrt_mtrx_4me <- as.matrix(cvrt_4me[, shared_samples])
    fwrite(cvrt_4me[, shared_samples], file = str_glue("covariates_fedtoMatrixEQTL_{run_flag}.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

    # nolint start
    cvrt <- SlicedData$new()
    cvrt$fileOmitCharacters <- "NA" # denote missing values;
    cvrt$fileSliceSize <- 2000;
    cvrt$CreateFromMatrix(cvrt_mtrx_4me)
    # nolint end
} else {
    shared_samples <- intersect(pntp_samples, gntp_samples)
    cvrt <- NULL
}

n_shared_samples <- length(shared_samples)
cat("Samples (", n_shared_samples, ") will be fed to MatrixEQTL:", shared_samples, "\n")

pntp_mtrx_4me <- as.matrix(pntp_4me[, shared_samples])
fwrite(pntp_4me[, shared_samples], file = str_glue("phenotypes_fedtoMatrixEQTL_{run_flag}.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

gntp_mtrx_4me <- as.matrix(gntp_4me[, shared_samples])

#
## QTL mapping
#
# nolint start
snps <- SlicedData$new()
snps$fileOmitCharacters <- "NA" # denote missing values;
snps$fileSliceSize <- 2000;
snps$CreateFromMatrix(gntp_mtrx_4me)
# nolint end

pntp_mtrx_4me_4pm <- pntp_mtrx_4me # Matrix for the permutations operation
use_model <- MatrixEQTL::modelLINEAR
pv_opt_thrd <- 0.9999  # TODO: add a CLI option for it
err_cov <- numeric()

pm_seed <- opts$pm_seed
set.seed(pm_seed)

pm_times <- opts$pm_times
for (pm in 0:pm_times) {
    # nolint start
    gene <- SlicedData$new()
    gene$fileOmitCharacters <- "NA" # denote missing values;
    gene$fileSliceSize <- 2000
    # nolint end

    if (pm) {
        cat("Permuted phenotype:", pm, "\n")
        qtls_opt_file <- as.character(str_glue("qtls_{run_flag}_permut{pm}.tsv"))
        pntp_mtrx_4me_4pm <- t(apply(pntp_mtrx_4me_4pm, 1, sample, size = n_shared_samples, replace = TRUE))
        gene$CreateFromMatrix(pntp_mtrx_4me_4pm)
    } else {
        cat("Raw phenotype\n")
        qtls_opt_file <- as.character(str_glue("qtls_{run_flag}.tsv"))
        gene$CreateFromMatrix(pntp_mtrx_4me)
    }

    me <- MatrixEQTL::Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = qtls_opt_file,
        pvOutputThreshold = pv_opt_thrd,
        useModel = use_model,
        errorCovariance = err_cov,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    )
}

#
## Report
#

# NOTE: There're repeated SNP IDs in the QTL mapping results due to multiple
# phentoypes (or genes), while the smtread() isn't smart enough to handle it.
qtls <- fread(qtls_opt_file, data.table = FALSE)

gntp_info_idx_col <- gntp_info_cols_vec[1]
chrom <- gntp_info_cols_vec[2]
position <- gntp_info_cols_vec[3]

gntp_info_dtfm <- fread(gntp_info_file, data.table = TRUE)
qlts_info_dtfm <- merge(qtls, gntp_info_dtfm, by.x = "SNP", by.y = gntp_info_idx_col)

mhtn_fig_p_thrd <- opts$mhtn_fig_p_thrd
for (pntp in target_pntp_vec) {
    pntp_qtls <- qlts_info_dtfm[qlts_info_dtfm$gene == pntp, c("SNP", chrom, position, "p-value")]

    if (nrow(pntp_qtls) > 0) {
        colnames(pntp_qtls) <- c("SNP", "CHR", "BP", "P")

        if (pv_opt_thrd < 0.9999) {
            lmb <- "NULL"
        } else {
            lmb <- lambda(pntp_qtls$P, "PVAL")
        }

        cat(str_glue("Inflation factor for {pntp}: {lmb}"), "\n")
        pdf(str_glue("MahattanPlot_{pntp}_{run_flag}.pdf"), width = 16, height = 9)
        manhattan(
            pntp_qtls[pntp_qtls[, "P"] <= mhtn_fig_p_thrd, ],
            main = str_glue("Manhattan plot for {pntp}"), ylab = "p-value(-log10)", annotateTop = TRUE
        )
        dev.off()

        pdf(str_glue("QQPlot_{pntp}_{run_flag}.pdf"), width = 16, height = 16)
        qq(pntp_qtls$P, main = str_glue("Q-Q plot {pntp} (lambda={lmb})"), pch = 18, cex = 1, las = 1)
        dev.off()
    } else {
        warning(str_glue("No QTL is found under p value ({pv_opt_thrd}) for {pntp} (trait or gene)"))
    }
}
