#!/usr/bin/env Rscript

# Create date: 2020-Jan-29
# Last update: Wed 29 Jan 2020 10:54:16 AM CET
# Version    : 0.1.0
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
#    5. For permutation test, the genotype should be shuffled simouteniously
#    with covariates.
#    6. Because R limits the each demision of a metrices up to 2^31-1, it's not
#    possbile to do permutation using all SNPS, if you want to exploit pbdMPI.
#    7. This script is meant to do permutation test, if your results are good
#    enough, you don't need it

# TODO:
#    1. A README.md to descript this shit.


library(pbdMPI)
library(stringr)
library(MatrixEQTL)

get_args <- function() {
    parser <- OptionParser(description = "Permutation using MatrixEQTL for QTL mapping.")
    parser <- add_option(parser, c("-w", "--work-dir"), action = "store", dest = "work_dir", type = "character", default = "./qtlmapping_opd", help = "Output direcotry.")
    parser <- add_option(parser, "--run-flag", action = "store", dest = "run_flag", type = "character", default = "test", help = "Running flag which help to discrimnate different runs.")

    # Phenotypes related parameters
    parser <- add_option(parser, c("-p", "--pntp-file"), action = "store", dest = "pntp_file", type = "character", help = "Phenotype input file.")
    parser <- add_option(parser, "--padding", action = "store", dest = "padding", type = "integer", default = 4, help = "The times of standard deviation as the boundary of outliers.")
    parser <- add_option(parser, "--target-pntp", action = "store", dest = "target_pntp", type = "character", help = "The phenotypes will be used, if more than one, using comma as delimiter.")
    parser <- add_option(parser, "--trps-pntp-dtfm", action = "store_true", dest = "trps_pntp_dtfm", help = "Whether should transpose the data.frame of phenotypes to make the columns as sample ID.")
    parser <- add_option(parser, "--pntp-idx-col", action = "store", dest = "pntp_idx_col", type = "character", default = "id", help = "The id column in phenotype file")

    # Covariates related parameters
    parser <- add_option(parser, c("-c", "--cvrt-file"), action = "store", dest = "cvrt_file", type = "character", help = "Covariates file.")
    parser <- add_option(parser, "--target-cvrt", action = "store", dest = "target_cvrt", type = "character", help = "Covariates will be used, if more than one, using comma as delimiter.")
    parser <- add_option(parser, "--trps-cvrt-dtfm", action = "store_true", dest = "trps_cvrt_dtfm", help = "Whether should transpose the data.frame of covariates.")
    parser <- add_option(parser, "--cvrt-idx-col", action = "store", dest = "cvrt_idx_col", type = "character", default = "id", help = "The id column in covariates file")

    # Genotypes related parameters
    parser <- add_option(parser, c("-d", "--gntp-dosage-file"), action = "store", dest = "gntp_dosage_file", type = "character", help = "The genotype dosage file (could be compressed).")
    parser <- add_option(parser, "--genotype-dosage-idx-col", action = "store", dest = "gntp_dosage_idx_col", type = "character", default = "id", help = "The id column in genotype file, usually its the name of column of SNP id. Default: id")
    parser <- add_option(parser, c("-i", "--gntp-info-file"), action = "store", dest = "gntp_info_file", type = "character", help = "The genotype information file (could be compressed).")
    parser <- add_option(parser, "--gntp-info-cols", action = "store", dest = "gntp_info_cols", type = "character", default = "rsID,SequenceName,Position,EffectAllele,AlternativeAllele", help = paste("The columns will be used.", "Default: rsID,SequenceName,Position,EffectAllele,AlternativeAllele"))
    parser <- add_option(parser, "--maf-thrd", action = "store", dest = "maf_thrd", type = "double", default = 0.05, help = "Minor allele frequency.")

    # Permutations
    parser <- add_option(parser, "--pm-times", action = "store", dest = "pm_times", type = "integer", default = 1000, help = "How many times of permutations should be done. If it's less than 1, no permutation will be performed but only 'raw' data will be used in the mapping. Default: 0")
    parser <- add_option(parser, "--pm-seed", action = "store", dest = "pm_seed", type = "integer", default = 31415, help = "The random seed for permutation. Defautl: 31415")

    return(parse_args2(parser))
}

jump_apply <- function(dtfm, func=min, win_size = 2) {
    n_cell <- base::nrow(dtfm)
    n_chunk <- n_cell / win_size
    base_idx <- 1 + ((1:n_chunk - 1) * win_size)

    optdtfm <- NULL
    for (x in 1: win_size) {
        tmpdtfm <- apply(dtfm[base_idx, ], 2, func)
        if(is.null(optdtfm)) {
            optdtfm <- tmpdtfm
        } else {
            optdtfm <- rbind(optdtfm, tmpdtfm)
        }

        base_idx <- base_idx + 1
    }

    rownames(optdtfm) <- 1:win_size

    return(optdtfm)
}

init()

comm_size <- comm.size()
comm_rank <- comm.rank()

if (file.exists("./utils.R")) {
    source("./utils.R")
}

chunk_data_to_scatter_pool <- NULL

if (comm_rank == 0) {
    library(optparse)
    #
    ## Parsing CLI arguments
    #
    opts <- get_args()$options

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
        dir.create(work_dir, recursive = T)
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
        pntp_chosen <- smtread(pntp_file, idxc = pntp_idx_col, kpr = target_pntp_vec, trps = T)
    }

    # Read covariates.
    cvrt_file <- opts$cvrt_file
    if (is.null(cvrt_file)) {
        with_cvrt <- FALSE
        warning("Not covariates is given, there won't be covariates in the model")
    } else {
        with_cvrt <- T

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
            cvrt_chosen <- smtread(cvrt_file, idxc = cvrt_idx_col, kpr = target_cvrt_vec, trps = T)
        }
    }

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
    pntp_4me <- as.data.frame(t(pntp_chosen))
    rownames(pntp_4me) <- colnames(pntp_chosen)
    colnames(pntp_4me) <- rownames(pntp_chosen)
    cat("Dim of phenotype data.frame will be piped into MatrixEQTL:", dim(pntp_4me), "\n")

    ## Covariates
    if (with_cvrt) {
        cvrt_4me <- as.data.frame(t(cvrt_chosen))
        rownames(cvrt_4me) <- colnames(cvrt_chosen)
        colnames(cvrt_4me) <- rownames(cvrt_chosen)
        cat("Dim of covariate data.frame will be piped into MatrixEQTL:", dim(cvrt_4me), "\n")
    } else {
        cat("[INFO] No covariate is supplied. Skipping covariates")
    }

    ## Genotypes
    gntp_dosage_idx_col <- opts$gntp_dosage_idx_col
    gntp_dosage <- smtread(gntp_dosage_file, idxc = gntp_dosage_idx_col)

    maf_thrd <- opts$maf_thrd
    gntp_round <- round(gntp_dosage) # For the phenotype level per genotype plot
    gntp_af <- apply(gntp_round, 1, sum) / (2 * dim(gntp_dosage)[[2]])
    nmvar_idx <- which((gntp_af >= maf_thrd) & (gntp_af <= 1 - maf_thrd))

    gntp_4me <- gntp_dosage[nmvar_idx, ]
    cat("Dim of genotype data.frame will be piped into MatrixEQTL:", dim(gntp_4me), "\n")

    pntp_samples <- colnames(pntp_4me)
    gntp_samples <- colnames(gntp_4me)

    if (with_cvrt) {
        cvrt_samples <- colnames(cvrt_4me)
        shared_samples <- intersect(pntp_samples, intersect(cvrt_samples, gntp_samples))
        cvrt_mtrx_4me <- as.matrix(cvrt_4me[, shared_samples])
    } else {
        shared_samples <- intersect(pntp_samples, gntp_samples)
        cvrt_mtrx_4me <- NULL
    }

    n_shared_samples <- length(shared_samples)
    cat("Samples (", n_shared_samples, ") will be fed to MatrixEQTL:", shared_samples, "\n")

    pntp_mtrx_4me <- as.matrix(pntp_4me[, shared_samples])
    gntp_mtrx_4me <- as.matrix(gntp_4me[, shared_samples])

    n_row <- base::nrow(gntp_mtrx_4me)
    chunk_size <- as.integer(n_row/comm_size) + 1
    chunk_bounds <- seq(1, n_row, chunk_size)

    chunk_data_to_scatter_pool <- list()
    for(x in chunk_bounds) {
        from <- x
        stop <- x + chunk_size - 1
        stop <- ifelse(stop >= n_row, n_row, stop)

        chunk_id <- paste0("chunck_", from, "_", stop)
        if(from == stop) {
            gntp_mtrx_4me_chunk <- as.matrix(gntp_mtrx_4me[from: stop, ])    
        } else {
            gntp_mtrx_4me_chunk <- gntp_mtrx_4me[from: stop, ]
        }

        chunk_data_to_scatter_pool[[chunk_id]] <- list(
            gntp=gntp_mtrx_4me_chunk, pntp=pntp_mtrx_4me, cvrt=cvrt_mtrx_4me,
            shared_samples=shared_samples, opts=opts, chunk_id=chunk_id
        )
    }
}

chunk_data <- scatter(chunk_data_to_scatter_pool)

#
## QTL mapping
#
gntp_mtrx_4me_chunk <- chunk_data$gntp
pntp_mtrx_4me <- chunk_data$pntp
cvrt_mtrx_4me <- chunk_data$cvrt
shared_samples <- chunk_data$shared_samples
opts <- chunk_data$opts
chunk_id <- chunk_data$chunk_id

cvrt <- NULL
if (!is.null(cvrt_mtrx_4me)) {
    cvrt <- SlicedData$new()
    cvrt$fileOmitCharacters <- "NA" # denote missing values;
    cvrt$fileSliceSize <- 2000;
    cvrt$CreateFromMatrix(cvrt_mtrx_4me)
}

snps <- SlicedData$new()
snps$fileOmitCharacters <- "NA" # denote missing values;
snps$fileSliceSize <- 2000;
snps$CreateFromMatrix(gntp_mtrx_4me_chunk)

gene <- SlicedData$new()
gene$fileOmitCharacters <- "NA" # denote missing values;
gene$fileSliceSize <- 2000

pv_opt_thrd <- 0.9999  # TODO: add a CLI option for it

set.seed(opts$pm_seed)

for (pm in 0: opts$pm_times) {
    if (pm) {
        shuffled_samples <- sample(x = shared_samples, replace = FALSE)
        pntp_mtrx_4me_4pm <- pntp_mtrx_4me[, shuffled_samples]
        gene$CreateFromMatrix(pntp_mtrx_4me_4pm)

        if (!is.null(cvrt)) {
            cvrt_mtrx_4me_4pm <- cvrt_mtrx_4me[, shuffled_samples]
            cvrt$CreateFromMatrix(cvrt_mtrx_4me_4pm)
        }
    } else {
        gene$CreateFromMatrix(pntp_mtrx_4me)
    }

    me <- MatrixEQTL::Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = NULL,
        pvOutputThreshold = pv_opt_thrd,
        useModel = MatrixEQTL::modelLINEAR,
        errorCovariance = numeric(),
        verbose = FALSE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = T,
        noFDRsaveMemory = FALSE
    )

    if (pm == 0) {
        min_gene_pv_dtfm <- as.data.frame(t(me$all$min.pv.gene))
        rownames(min_gene_pv_dtfm) <- str_glue("raw_rk{comm_rank}")
    } else {
        new_row_names <- c(rownames(min_gene_pv_dtfm), str_glue("pm{pm}_rk{comm_rank}"))
        min_gene_pv_dtfm <- rbind(min_gene_pv_dtfm, as.data.frame(t(me$all$min.pv.gene)))
        rownames(min_gene_pv_dtfm) <- new_row_names
    }
}

chunk_data_to_gather <- list(min_gene_pv=min_gene_pv_dtfm, chunk_id=chunk_id)
chunk_data_be_gather <- gather(chunk_data_to_gather)

#
## Dump permutation results.
#
if (comm_rank == 0) {
    library(data.table)

    min_gene_pv_dtfm <- NULL
    for (chunk in chunk_data_be_gather) {
        if (is.null(min_gene_pv_dtfm)) {
            min_gene_pv_dtfm <- chunk$min_gene_pv
        } else {
            min_gene_pv_dtfm <- rbind(min_gene_pv_dtfm, chunk$min_gene_pv)
        }
    }

    pm_pval_dtfm <- jump_apply(min_gene_pv_dtfm, func=min, win_size=opts$pm_times + 1)
    pm_pval_opt_file <- as.character(str_glue("permutation_pval_{run_flag}.tsv"))
    fwrite(pm_pval_dtfm, file = pm_pval_opt_file, sep = "\t", col.names = T, row.names = T)
}

finalize()
