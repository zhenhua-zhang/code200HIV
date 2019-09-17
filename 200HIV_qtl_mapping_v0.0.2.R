#!/usr/bin/env Rscript
# Using R 3.3.3
################################################################################
#' @version: 0.2.0
#' @author: Zhenhua Zhang <zhenhua.zhang217@gmail.com>
#' @date: 4th June, 2019
#' @reference: qtl_mapping.R (Xiaojing Chu)
################################################################################


# TODO: 1 When doing inverse-rank transform, use mean to replace NA, which is
#         dangerous. Need a better way to deal with NAs
#       2 Make the script rescuable at broken-points with additional command
#         arguments
#       3 Make the LocusZoom analysis local

# NOTE: 1 Settings in config file will be overrided by command line arguments
#       2 Some settings is only modifiable by configuration files
#       3 Keep the row as measurements, column as subjects
#       4 Please try to avoid using comma in each field, for the saperator is
#         comman in .csv files


# Used packages
# TODO: 1. Add version for each library
# library(docstring)
library(optparse, warn.conflicts = FALSE, quietly = TRUE)


# Processing the configuration file
prcs_cfg <- function(iptfl = NULL) {
  if (is.null(iptfl)) stop("iptfl is required but found NULL")

  cfglst <- list()
  cfgs <- paste(readLines(iptfl), sep = "\n")
  for (line in cfgs) {
    if (!startsWith(line, "#") && !stri_isempty(line)) {
      plst <- strsplit(line, "=")
      if (length(plst) != 1) stop("The length of plst should be 1")

      key <- plst[[1]][1]
      val <- plst[[1]][2]
      cfglst[[key]] <- val
    }
  }
  return(cfglst)
}


# Read file using fread() from data.table. Dealing with zipped files
smt_fread <- function(flnm, zipped = FALSE, data.table = FALSE, header = TRUE, stringsAsFactors = FALSE, verbose = FALSE, showProgress = FALSE) {
  if (zipped) {
    return(fread(cmd = paste0("zcat < ", flnm), data.table = data.table,
      header = header, stringsAsFactors = stringsAsFactors,
      verbose = verbose, showProgress = showProgress
    ))
  } else {
    return(fread(flnm, data.table = data.table, header = header,
      stringsAsFactors = stringsAsFactors, verbose = verbose,
      showProgress = showProgress
    ))
  }
}


# Preprocess
pprcs <- function(flnm, dscd_cols = NULL, dscd_rows = NULL, kept_cols = NULL,
  kept_rows = NULL, trps = TRUE, as_idx = "id", dscd_idx = FALSE, zipped = FALSE,
  swhd = 10) {
  cat("[INFO]  Processing:", flnm, "...\n")

  dtfm <- smt_fread(flnm, zipped = zipped)

  if (!is.null(as_idx)) {
    base::row.names(dtfm) <- dtfm[, as_idx]
    if (dscd_idx) dtfm[as_idx] <- NULL
    }

  col_names <- base::colnames(dtfm)
  dtfm <- dtfm[, base::sort(col_names)]
  col_names <- base::colnames(dtfm)

  if (is.null(kept_cols)) kept_cols <- col_names[which(!col_names %in% dscd_cols)]

  cat("[INFO]  Kept cols (first", min(swhd, length(kept_cols)), "of", length(kept_cols), "):", head(kept_cols, n = swhd), "...\n")
  dtfm <- dtfm[, kept_cols]

  row_names <- base::rownames(dtfm)
  if (is.null(kept_rows)) kept_rows <- row_names[which(!row_names %in% dscd_rows)]

  cat("[INFO]  Kept rows (first", min(swhd, length(kept_rows)), "of", length(kept_rows), "):", head(kept_rows, n = swhd), "...\n")
  dtfm <- dtfm[kept_rows, ]

  if (trps) {
    if (!is.null(as_idx)) {
      cat("[WARN]  Will transpose dataframe taking", as_idx, "as character\n")
      if (dscd_idx || (as_idx %in% dscd_cols)) {
        dtfm <- as.data.frame(t(dtfm), stringsAsFactors = FALSE)
      } else {
        idx_val <- dtfm[as_idx]
        dtfm <- dtfm[, which(!base::colnames(dtfm) %in% c(as_idx))]
        dtfm <- as.data.frame(t(dtfm), stringsAsFactors = FALSE)
        dtfm[as_idx, ] <- idx_val
        dtfm <- dtfm[c(as_idx, base::rownames(dtfm))]
      }
    } else {
      dtfm <- as.data.frame(t(dtfm), stringsAsFactors = FALSE)
    }
  }

  return(dtfm)
}


# Check the normality of each measuements
chck_nmlt <- function(dtfm, hst_opt = "Normality_check.pdf", bins = 30, idv = "id", stclnm = TRUE, fgsz = 3) {
  #' @title
  #' @param stclnm 以后得记得抓紧时间写上注释, 现在可好, 已经忘了为啥这么写了
  cat("[INFO]  Check the normality...\n")

  if (!is.data.frame(dtfm)) stop("`dtfm` should be a dataframe ...")

  if (stclnm) dtfm <- dtfm[, base::sort(base::colnames(dtfm))]

  exp_mlt <- reshape2::melt(dtfm, id.vars = idv)

  nvars <- length(unique(exp_mlt["variable"])[[1]])
  ncol <- as.integer(sqrt(nvars))

  p <- ggplot(data = exp_mlt)
  p <- p + geom_histogram(aes(x = value), na.rm = TRUE, bins = bins)
  p <- p + facet_wrap(. ~ variable, ncol = ncol, scales = "free")

  ggsave(hst_opt, width = ncol * fgsz, height = ncol * fgsz)
}


# Transform
ivsrk <- function(x, lude.mode = TRUE) {
  # Inverse rank normalization
  is_na <- is.na(x)
  x[is_na] <- mean(x, na.rm = TRUE)
  if (lude.mode) {
    res <- base::rank(x)
    res <- qnorm((0.5 + res - 1) / (length(res)))
    res <- (res * sd(x)) + mean(x)
  } else {
    res <- base::rank(x)
    res <- qnorm(res / (length(res) + 0.5))
    res <- (res * sd(x)) + mean(x)
  }
  res[is_na] <- NA

  return(res)
}


# Trans form log10 values into log2 values
lg10tolg2 <- function(x) { return(x * log2(10)) }

# Trans form log2 values into log10 values
lg2tolg10 <- function(x) { return(x * log10(2)) }


# NOTE: Make sure only the id column is characters, while the rest are numeric,
#       otherwise, the transpose would fail
# Transform non-normal distributed phnotype column into
trfm_nmlt <- function(dtfm, dscd_cols = NULL, trfm_cols = NULL, trfm_mthd = "log2", idcol = "id", stclnm = TRUE) {
  #' @title
  #' @param stclnm 忘了这是用来干啥的了
  cat("[INFO]  Transform data by", trfm_mthd, "\n")

  if (!is.null(idcol)) {
    id_val <- dtfm[idcol]
    trfm_cols <- trfm_cols[which(!trfm_cols %in% c(idcol))]
    dtfm[idcol] <- NULL
  }

  row_names <- base::rownames(dtfm)
  if (stclnm) dtfm <- dtfm[base::sort(base::colnames(dtfm))]
  col_names <- colnames(dtfm)

  if (!is.null(dscd_cols)) {
    cat("[WARN]  Will discard: \n         ", dscd_cols, "\n")
    dtfm <- dtfm[, which(!col_names %in% dscd_cols)]
    col_names <- base::colnames(dtfm)
  }

  if (!is.null(trfm_cols) && !is.na(trfm_cols)) {
    if (!is.null(dscd_cols)) {
      ovlp <- trfm_cols[which(dscd_cols %in% trfm_cols)]
      if (length(ovlp) != 0) stop("Overlapped between dscd_cols and trfm_cols \n", ovlp)
      }

    if (is.null(trfm_mthd)) stop("When trfm_cols is not NULL, trfm_mthd is required.")

    if (length(trfm_cols) == 1) dtfm[, trfm_cols] <- get(trfm_mthd)(dtfm[, trfm_cols])
    else dtfm[, trfm_cols] <- sapply(dtfm[, trfm_cols], trfm_mthd)

    base::colnames(dtfm) <- sapply(col_names,
      function(x) {
        if (x %in% trfm_cols) return(paste(x, trfm_mthd, sep = "_"))
        else return(x)
        }
    )
  } else {
    cat("[WARN]  trfm_cols is either NULL or NA\n")
  }

  dtfm <- as.data.frame(scale(dtfm), stringsAsFactors = FALSE)
  if (!is.null(idcol)) dtfm[idcol] <- id_val
  return(dtfm)
}


# Draw PCA
plt_pca <- function(dtfm, opt_flnm = "Outlier_check_PCA.pdf", idcol = "id", fgwd = 10, fght = 10, fgunt = "in", n_sd = 3) {
  #' @title Draw PCA plot
  #' @description Draw a PCA plot to check the outliers.
  #' @param dtfm A dataframe.

  cat("[INFO]  Plotting PCA... \n")

  if (!is.data.frame(dtfm)) stop("`dtfm` should be a dataframe ...")

  if (!is.null(idcol)) {
    dtfm_id <- dtfm[idcol]
    dtfm[idcol] <- NULL
  }

  dtfm[is.na(dtfm)] <- mean(base::rowMeans(dtfm), na.rm = TRUE)
  pca <- prcomp(dtfm)
  pcax <- as.data.frame(pca$x, stringsAsFactors = FALSE)

  if (!is.null(idcol)) pcax["id"] <- dtfm_id

  pc1 <- pcax[, "PC1"]
  pc1mn <- mean(pc1)
  pc1sd <- sd(pc1)

  pc2 <- pcax[, "PC2"]
  pc2mn <- mean(pc2)
  pc2sd <- sd(pc2)

  pc1otl <- ((pcax$PC1 > pc1mn + n_sd * pc1sd) | (pcax$PC1 < pc1mn - n_sd * pc1sd))
  pc2otl <- ((pcax$PC2 > pc2mn + n_sd * pc2sd) | (pcax$PC2 < pc2mn - n_sd * pc2sd))

  pcac <- pca$center

  p <- ggplot() + theme_bw()
  p <- p + geom_point(data = pcax, mapping = aes(x = PC1, y = PC2))
  p <- p + geom_text(data = pcax[which(pc1otl | pc2otl), ], mapping = aes(x = PC1, y = PC2, label = id))

  ggsave(opt_flnm, width = fgwd, height = fght, units = fgunt)
}


# Remove outliers
prcs_otls <- function(dtfm, n_sd = 3, idcol = "id") {
  cat("[INFO]  Removing outliers...\n")

  if (!is.null(idcol)) {
    id_val <- dtfm[idcol]
    dtfm[idcol] <- NULL
  }

  row_names <- base::rownames(dtfm)
  col_names <- base::colnames(dtfm)

  dtfm <- sapply(dtfm,
    function(x) {
      sd <- sd(x, na.rm = TRUE)
      mu <- mean(x, na.rm = TRUE)
      x[which((x <= mu - n_sd * sd) | (x >= mu + n_sd * sd))] <- NA
      return(x)
    }
  )

  dtfm <- as.data.frame(dtfm, row.names = row_names, stringsAsFactors = FALSE)
  cat("[INFO]  Number of masked subjects in each measurments:\n")
  nacnt <- sapply(dtfm, function(x) sum(is.na(x)))

  col_names <- base::colnames(dtfm)

  dtfm[idcol] <- id_val
  dtfm <- dtfm[c(idcol, col_names)]
  return(dtfm)
}


# Create SliceData object
mkslcdt <- function(ipt, dlmt = ",", omchr = "NA", skprw = 0, skpcl = 0, slice_size = 2000) {
  # cat("[INFO]  Creating SlicedData ...\n")
  slcdt <- SlicedData$new()
  slcdt$fileDelimiter <- dlmt
  slcdt$fileOmitCharacters <- omchr
  slcdt$fileSkipRows <- skprw
  slcdt$fileSkipColumns <- skpcl
  slcdt$fileSliceSize <- slice_size

  if (is.null(ipt)) stop("Input `ipt` should not be NULL!!")
  else if (is.matrix(ipt)) slcdt$CreateFromMatrix(ipt)
  else if (is.data.frame(ipt)) slcdt$CreateFromMatrix(as.matrix(ipt))
  else slcdt$LoadFile(ipt)

  return(slcdt)
}


# Main method to do QTL mapping
mkme <- function(pntp, gntp, cvrt = NULL, opt_file = NULL, pv_opt_thr = 5e-2) {
  cat("[INFO]  Mapping QTL ...\n")
  useModel <- modelLINEAR
  errorCovariance <- numeric()

  gene <- mkslcdt(pntp)
  snps <- mkslcdt(gntp, dlmt = "\t")
  cvrt <- mkslcdt(cvrt)

  me <- Matrix_eQTL_engine(
    snps = snps, gene = gene, cvrt = cvrt,
    output_file_name = opt_file, pvOutputThreshold = pv_opt_thr,
    useModel = useModel, errorCovariance = errorCovariance, verbose = FALSE,
    pvalue.hist = TRUE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE
  )

  return(me)
}


# Intersect Matrix eQTL with SNP information
itsct_qtlmt_ifmt <- function(qtls, snp_ifmt, mgby_x = "snps", mgby_y = "rsID") {
  cat("[INFO]  Intersect QTL matrix with SNP information ...\n")
  if (is.character(qtls))
    qtls <- smt_fread(qtls)
  if (is.character(snp_ifmt))
    snp_ifmt <- smt_fread(snp_ifmt)
  return(merge(qtls, snp_ifmt, by.x = mgby_x, by.y = mgby_y))
}


# Draw Manhattan plot and Q-Q plot
plt_mht <- function(qtls_lst, use_cols = NULL, opt_prfx = "qtl_mapping", svfmt = "pdf", wdth = 20, hght = 10, drwzsc = FALSE, drwqq = FALSE) {
  cat("[INFO]  Plotting Manhattan plot and Q-Q plot for", opt_prfx, "...\n")
  # use_cols should follow the order: snp_id, chromosome, position, p-value
  if (is.null(use_cols)) use_cols <- c("snps", "chr", "pos", "pvalue")
  else if (length(use_cols) != 4) stop("The length of use_cols should be 4 ...\n")

  gwas_rslt <- qtls_lst[, use_cols]
  base::colnames(gwas_rslt) <- c("SNP", "CHR", "BP", "P")
  opt_flnm <- paste0(opt_prfx, "_mhtn_pvalue.", svfmt)

  if (svfmt == "png") png(opt_flnm, width = wdth, height = hght)
  else pdf(opt_flnm, width = wdth, height = hght)

  manhattan(
    gwas_rslt, main = "Manhattan Plot (p-value)", ylab = "p-value(-log10)",
    suggestiveline = FALSE, annotatePval = 5e-8, annotateTop = FALSE
  )
  dev.off()

  if (drwzsc) {
    gwas_rslt <- transform(gwas_rslt, zscore = qnorm(P / 2, lower.tail = FALSE))
    opt_flnm <- paste0(opt_prfx, "_mhtn_zscore.", svfmt)

    if (svfmt == "png") png(opt_flnm, width = wdth, height = hght)
    else pdf(opt_flnm, width = wdth, height = hght)

    manhattan(gwas_rslt, main = "Manhattan Plot (z-score)", ylab = "Z-score", suggestiveline = FALSE, annotatePval = 0.00000005)
    dev.off()
  } else {
    cat("[INFO]  Skipping Manhattan plot for z-score ...\n")
  }

  if (drwqq) {
    opt_flnm <- paste0(opt_prfx, "_qq.", svfmt)

    if (svfmt == "png") png(opt_flnm, width = wdth, height = hght)
    else pdf(opt_flnm, width = wdth, height = hght)

    qq(gwas_rslt$P, main = "Q-Q plot of p-values", pch = 18, cex = 1.5, las = 1)
    dev.off()
  } else {
    cat("[INFO]  Skipping Q-Q plot ...\n")
  }
}


# Filter out by checking the number of mode
ckmd <- function(dtfm, thrshld = 0.5) {
  if (is.data.frame(dtfm)) {
    res <- sapply(dtfm,
      function(x) {
        cnt <- table(x)
        blvct <- ifelse(max(cnt) / sum(cnt) >= thrshld, TRUE, FALSE)
        return(blvct)
      })
    return(res)
  } else if (is.vector(dtfm)) {
    cnt <- table(dtfm)
    return(ifelse(max(cnt) / sum(cnt) >= thrshld, TRUE, FALSE))
  } else {
    stop("Can only handle vector or data.frame!!!")
  }
}

# Encode genotype
enc_gntp <- function(x, use_gntp, use_pntp) {
  #' @title enc_gntp
  #' @description Encode dosage level into genotypes
  #' @param x vector. A vector including SNP information
  #' @param use_gntp vector. A vector including SNPs will be processes
  #' @param use_pntp vector.

  snpid <- x[1]
  msmnt <- x[2]
  ref <- x[6]
  alt <- x[7]
  gntp_dsg <- use_gntp[snpid, ]
  gntp_code <- sapply(gntp_dsg,
    function(d) {
      if (length(d) != 1) return(NA)
      if (!is.numeric(d)) return(d)
      else if (round(d) == 0) return(paste0(ref, ref))
      else if (round(d) == 1) return(paste0(ref, alt))
      else if (round(d) == 2) return(paste0(alt, alt))
      else return(NA)
      })
  pntp_msmnt <- use_pntp[msmnt, ]
  gnpn <- rbind(pntp_msmnt, gntp_dsg, gntp_code)
  base::rownames(gnpn) <- c(msmnt, snpid, "code")
  return(gnpn)
}


# REPORT section
# Make VCF param for a readVcf
mkvcfprm <- function(rgstr = NULL, chrom = NULL, start = NULL, stop = NULL, onebased = TRUE) {
  # FIXME: unused parameter onebased
  #' @title mkvcfprm
  #' @version 0.1.0
  #' @description Create ScanVcfParam object
  #' @param rgstr Character.
  #' @param chrom Character.
  #' @param start Integer.
  #' @param stop Integer.
  #' @param onebased Boolean.
  #' @return ScanVcfParam.

  if (is.null(c(rgstr, chrom))) {
    stop("rgstr and chrom could not be NULL at the same time.")
  } else if (is.null(chrom)) {
    if (!is.vector(rgstr)) { rgstr <- c(rgstr) }
    for (sub_rgstr in rgstr) {
      rglst <- stri_split(sub_rgstr, regex = ":|-", simplify = TRUE)
      if (length(rglst) <= 3) {
        chrom <- c(chrom, rglst[1])
        start <- c(start, ifelse(length(rglst) < 2, 1, as.numeric(rglst[2])))
        stop <- c(stop, ifelse(length(rglst) < 3, start[length(start)], as.numeric(rglst[3])))
      } else {
        stop("Error while processing region string: ", rgstr)
      }
    }
  } else {
    chrom <- c(as.character(chrom))
    if (is.null(start)) start <- rep(1, length(chrom))
    if (is.null(stop)) stop <- start
  }

  if (length(chrom) != length(start) || length(start) != length(stop))
    stop("start and stop should have the same length ...")
  return(ScanVcfParam(which = GRanges(chrom, IRanges(start, stop))))
}


# Replicates the GW signicant SNPs in another independent GWAS
# rplct <- function(topsnp, smrst, pvt=5e-8) {
#' @title rplct(1)
#' @description Replicates genome wide significant SNPs from current GWAS in another independent GWAS
#'
#' @param topsnp dataframe. Top SNPs in current GWAS
#' @param smrst dataframe. Summary statistics from another independent GWAS
#' @return NULL
# }


# DuMP Top Snp information in VCF format
dmptsvcf <- function(tspc, vcfdb, vcfrfg, opt_dir = "reports", kept_info = NULL, dscd_info = NULL) {
  # TODO: A parameters in configuration file to decide add which INFO into the report.csv
  #' @title Dump top SNP information as a VCF file to disk
  #' @param tspc dataframe. Top SNPs per chromosome

  vcfprm <- mkvcfprm(apply(tspc, 1, function(x) paste0(x["chr"], ":", as.character(x["pos"]))))

  tspc_vcf <- readVcf(TabixFile(vcfdb), vcfrfg, vcfprm)
  vcfrpt_path <- file.path(opt_dir, "topSNPReport.vcf.gz")
  writeVcf(tspc_vcf, vcfrpt_path, index = TRUE)
}


# Fetch Top Snp Per Chromosome
ftspc <- function(qtl_dtfm, pvt = 5e-8, qtl_col = NULL) {
  #'
  #' @param qtl_col vector. Columns will be used, should follow <snps, gene, pvalue, chr, pos, ref, alt>
  #'
  if (is.null(qtl_col))
    qtl_col <- c("snps", "gene", "pvalue", "chr", "pos", "ref", "alt")

  qtl_dtfm <- qtl_dtfm[, qtl_col]
  qtl_col <- c("snps", "gene", "pvalue", "chr", "pos", "ref", "alt")
  base::colnames(qtl_dtfm) <- qtl_col

  qtl_gb <- c("snps", "chr", "gene")
  # Filter QTL datafram
  if (is.vector(qtl_col) && length(qtl_col) == 7) {
    tspc <- qtl_dtfm %>%
      filter(pvalue <= pvt) %>%
      group_by(chr, gene) %>%
      filter(pvalue == min(pvalue)) %>%
      ungroup() %>%
      as.data.frame()
  } else {
    stop("qtl_col should be a vector with length 7 ...")
  }

  return(tspc)
}


# LoaD genotype DoSaGe files
ldgntp <- function(gntp_ipt, gntp_ptn = NULL, gntp_zip = FALSE, gntp_idx_col = NULL,
  gntp_dscd_idx = NULL, dscd_rows = NULL, dscd_cols = NULL, kept_rows = NULL,
  kept_cols = NULL) {
  #' @title
  #' @description

  if (file_test("-d", gntp_ipt)) {
    if (is.null(gntp_ptn)) stop("When gntp_ipt is a dir, gntp_ptn is required")

    flnmlst <- list.files(gntp_ipt, pattern = gntp_ptn)
    if (length(flnmlst) == 0)
      stop("gntp_ipt is a dir, but NO file following pattern: ", gntp_ptn)

    gntp_dtfm <- data.frame()
    for (flnm in flnmlst) {
      gntp_file <- file.path(gntp_ipt, flnm)
      tmp_dtfm <- pprcs(gntp_file, trps = FALSE, zipped = gntp_zip,
        as_idx = gntp_idx_col, dscd_idx = gntp_dscd_idx, dscd_rows = dscd_rows,
        dscd_cols = dscd_cols, kept_rows = kept_rows, kept_cols = kept_cols
      )
      gntp_dtfm <- rbind(gntp_dtfm, tmp_dtfm)
    }
  } else if (file_test("-f", gntp_ipt)) {
    gntp_file <- gntp_ipt
    gntp_dtfm <- pprcs(gntp_file, trps = FALSE, zipped = gntp_zip,
      as_idx = gntp_idx_col, dscd_idx = gntp_dscd_idx, dscd_rows = dscd_rows,
      dscd_cols = dscd_cols, kept_rows = kept_rows, kept_cols = kept_cols
    )
  } else {
    stop("gntp_ipt should be a directory, normal file, or zipped file")
  }

  return(gntp_dtfm)
}


# Draw dosage boxplot
plt_dsg <- function(gnpnlst, svfmt = "pdf", opt_dir = "reports/") {
  #' @title plot dosage boxplot
  #' @description Function to draw the dosage level for each genotype

  for (gnpn in gnpnlst) {
    gnpn_t <- as.data.frame(t(gnpn))

    col_names <- base::colnames(gnpn_t)
    lvl <- col_names[1]
    snp <- col_names[2]
    opt_flnm <- file.path(opt_dir, paste0(lvl, "_", snp, ".", svfmt))

    base::colnames(gnpn_t) <- c("lvl", "dsg", "code")

    gnpn_t["lvl"] <- as.numeric(gnpn_t[, "lvl"])
    gnpn_t["dsg"] <- as.numeric(gnpn_t[, "dsg"])

    p <- ggplot(gnpn_t) + theme_bw()
    p <- p + geom_boxplot(aes(x = code, y = lvl, color = code))
    p <- p + geom_point(aes(x = code, y = lvl))
    ggsave(opt_flnm)
  }
}

#
# END REPORT section
#

#
# Utils starts
#
# Split arguments with coma or semi-colon in to a vector
hdmtopt <- function(arg) {
  if (is.na(arg)) return(cat("[WARN]  Keyword argument arg is NA, return NULL\n"))
  if (is.null(arg)) return(cat("[WARN]  Keyword argument arg is NULL\n"))
  if (length(arg) == 0) return(cat("[WARN]  Keyword argument arg is empty, return NULL\n"))

  blnks <- c(" ", "\t", "\n") # Check blanks
  if (any(stri_detect_fixed(arg, blnks))) stop("Blank or tab is not allowed!!")

  arg <- stri_split_fixed(arg, ";")[[1]]
  arg <- stri_split_fixed(arg, ",")
  if (length(arg) == 1) arg <- arg[[1]]

  return(arg)
}


# Parsing command line arguments
prsarg <- function() {
  parser <- OptionParser()

  parser <- add_option(parser, c("-G", "--genotypes"), type = "character", help = "Genotypes file or directory")
  parser <- add_option(parser, c("-C", "--covariates"), type = "character", help = "Covariates")
  parser <- add_option(parser, c("-P", "--phenotypes"), type = "character", help = "Phenotypes file")
  parser <- add_option(parser, c("-S", "--config-file"), type = "character", default = "config.cfg", help = "File including settings, overrided by command arguments")
  parser <- add_option(parser, c("-T", "--tmp-dir"), type = "character", help = "File for temporary files, will cleanup if no error is catched")
  parser <- add_option(parser, c("-o", "--output-dir"), type = "character", help = "Output directory for the QTL mapping results")
  agmnts <- parse_args(parser, positional_arguments = 1, convert_hyphens_to_underscores = TRUE)

  # cknmlt: check normality; trfm: Transform data;
  # qltmp: QTL mapping; qtlrpt: produce QTL report.
  subcmd <- agmnts$args
  subcmdlst <- c("cknmlt", "trfm", "qtlmp", "qtlrpt")
  if (!subcmd %in% subcmdlst) stop("Unknown: ", subcmd, ". Opts: cknmlt, trfm, qtlmp, qtlrpt")

  return(agmnts)
}


# Setup working space
stwkspc <- function(tmp_dir = NULL, opt_dir = NULL, mode = "0755") {
  # OPT_DIR
  # |-- logs
  # |-- reports
  # |   |-- GenotypeLevelPlot
  # |   |-- QtlInformation
  # |   |-- LocusZoomPlot
  # |   |-- ManhattanPlot
  # |   |-- Preprocessing
  # |   |-- Pathway
  # |-- report.tsv  TODO: 还没决定是不是要用这个文件名

  if (dir.exists(tmp_dir)) cat("[WARN] ", tmp_dir, "exists\n")
  else dir.create(tmp_dir, mode = mode)

  if (dir.exists(opt_dir)) {
    cat("[WARN] ", opt_dir, "exists\n")
  } else {
    dir.create(opt_dir, mode = mode)
    dir.create(file.path(opt_dir, "logs"), mode = mode)
    dir.create(file.path(opt_dir, "reports"), mode = mode)
    dirs <- c("GenotypeLevelPlot", "QtlInformation", "LocusZoomPlot", "ManhattanPlot", "Preprocessing", "Pathway")
    for (dir in dirs)
      dir.create(file.path(opt_dir, "reports", dir), mode = mode)
    }
}


# Merge the settings from configuration file and command line arguments
mgcfgcmd <- function(agmntlst, cfglst) {
  if (!is.list(agmntlst) || !is.list(cfglst)) stop("Keyword arguments both agmntlst and cfglst should be list")

  cmdopts <- agmntlst$options
  cmdopts["subcmd"] <- agmntlst$args
  cfglst[names(cmdopts)] <- cmdopts
  return(cfglst)
}

#
# Utils ends
#


#
# Start of Sub-commands
#
# sub-command `chnmlt`
cknmlt <- function(optlst) {
  tmp_dir <- optlst$tmp_dir
  opt_dir <- optlst$opt_dir

  # Arguments
  pntp_file <- optlst$phenotypes
  pntp_zip <- as.logical(optlst$pntp_zip)
  pntp_idx_col <- optlst$pntp_idx_col
  pntp_dscd_idx <- as.logical(optlst$pntp_dscd_idx)
  pntp_dtfm_trp <- as.logical(optlst$pntp_dtfm_trp)
  pntp_dscd_cols <- hdmtopt(optlst$pntp_dscd_cols)
  pntp_dscd_rows <- c(hdmtopt(optlst$pntp_dscd_rows), optlst$glbl_blck_lst)

  # Preprocessing
  pntp_dtfm <- pprcs(pntp_file, dscd_cols = pntp_dscd_cols, dscd_rows = pntp_dscd_rows, trps = FALSE, as_idx = pntp_idx_col, dscd_idx = pntp_dscd_idx)

  # Check normality
  pntp_raw_nmltCk_flnm <- file.path(optlst$opt_pprc_dir, optlst$pntp_raw_nmltCk_flnm)
  chck_nmlt(pntp_dtfm, hst_opt = pntp_raw_nmltCk_flnm)
  cat("[INFO]  Check", pntp_raw_nmltCk_flnm, "for distributions.\n")

  # Save the preprocessed dataframe for the next step.
  cknmlt_tmp_dtfm <- file.path(tmp_dir, "cknmlt_tmp.csv")
  fwrite(pntp_dtfm, cknmlt_tmp_dtfm, showProgress = FALSE, verbose = FALSE)
  cat("[INFO]  Check", cknmlt_tmp_dtfm, "for preprocessed input data.\n")
}


# sub-command `trfm`
trfm <- function(optlst) {
  tmp_dir <- optlst$tmp_dir
  opt_dir <- optlst$opt_dir

  # Read temp-file from the last step
  ipt_tmp_flnm <- file.path(tmp_dir, "cknmlt_tmp.csv")
  if (file.exists(ipt_tmp_flnm)) pntp_dtfm <- pprcs(ipt_tmp_flnm, trps = FALSE)
  else stop("Failed to find temporary file ", ipt_tmp_flnm)

  # Transform the raw data by specific method: log2, log10, inverse-rank
  pntp_trfm_mthd <- hdmtopt(optlst$pntp_trfm_mthd)
  pntp_trfm_cols <- hdmtopt(optlst$pntp_trfm_cols)
  pntp_trfm_excols <- hdmtopt(optlst$pntp_trfm_excols)

  if (pntp_trfm_cols == "*") {
    pntp_trfm_cols <- colnames(pntp_dtfm)
    excld <- which(!pntp_trfm_cols %in% pntp_trfm_excols)
    pntp_trfm_cols <- pntp_trfm_cols[excld]
  }

  pntp_dtfm <- trfm_nmlt(pntp_dtfm, trfm_cols = pntp_trfm_cols, trfm_mthd = pntp_trfm_mthd)

  ## Check the normality of the transformed raw data
  pntp_trfm_nmltCk_flnm <- file.path(optlst$opt_pprc_dir, optlst$pntp_trfm_nmltCk_flnm)
  chck_nmlt(pntp_dtfm, hst_opt = pntp_trfm_nmltCk_flnm)
  cat("[INFO]  Check", pntp_trfm_nmltCk_flnm, "for distributions after transforming.\n")

  ## Check outliers of transformed data using PCA
  pntp_trfm_otlCk_flnm <- file.path(optlst$opt_pprc_dir, optlst$pntp_trfm_otlCk_flnm)
  plt_pca(pntp_dtfm, pntp_trfm_otlCk_flnm, idcol = "id")
  cat("[INFO]  Check", pntp_trfm_otlCk_flnm, "for PCA plots.\n")

  ## Exclude outliers using 3 or 4 times of sd
  pntp_otl_n_sd <- as.numeric(optlst$pntp_otl_n_sd)
  pntp_dtfm <- prcs_otls(pntp_dtfm, n_sd = pntp_otl_n_sd)

  ## Check normality of dataframe without outliers
  pntp_exotl_nmltCk_flnm <- file.path(optlst$opt_pprc_dir, optlst$pntp_exotl_nmltCk_flnm)
  chck_nmlt(pntp_dtfm, hst_opt = pntp_exotl_nmltCk_flnm)
  cat("[INFO]  Check", pntp_exotl_nmltCk_flnm, "for distributions after excluding outliers.\n")

  ## Save processed data frame into temporary file
  opt_tmp_file <- file.path(optlst$tmp_dir, "trfm_tmp.csv")
  fwrite(pntp_dtfm, opt_tmp_file, showProgress = FALSE, verbose = FALSE)
  cat("[INFO]  Check", opt_tmp_file, " for preprocessed input data.\n")
}


# sub-command `qtlmp`
qtlmp <- function(optlst) {
  tmp_dir <- optlst$tmp_dir
  opt_dir <- optlst$opt_dir
  glbl_blck_lst <- optlst$glbl_blck_lst

  ipt_tmp_flnm <- file.path(tmp_dir, "trfm_tmp.csv")
  if (file.exists(ipt_tmp_flnm))
    pntp_dtfm <- pprcs(ipt_tmp_flnm, dscd_cols = c("id"), trps = TRUE)
  else
    stop("Failed to fined temporary file ", ipt_tmp_flnm)

  # Covariates
  cvrt_ipt <- optlst$covariates
  cvrt_trp <- as.logical(optlst$cvrt_trp)
  cvrt_zip <- as.logical(optlst$cvrt_zip)
  cvrt_idx_col <- optlst$cvrt_idx_col
  cvrt_dscd_idx <- as.logical(optlst$cvrt_dscd_idx)
  cvrt_dtfm_trp <- as.logical(optlst$cvrt_dtfm_trp)
  cvrt_dscd_cols <- hdmtopt(optlst$cvrt_dscd_cols)
  cvrt_kept_cols <- hdmtopt(optlst$cvrt_kept_cols)
  cvrt_dscd_rows <- c(glbl_blck_lst, hdmtopt(optlst$cvrt_dscd_rows))
  cvrt_kept_rows <- hdmtopt(optlst$cvrt_kept_rows)

  cvrt_dtfm <- pprcs(
    cvrt_ipt, dscd_cols = cvrt_dscd_cols, dscd_rows = cvrt_dscd_rows,
    zipped = cvrt_zip, as_idx = cvrt_idx_col, dscd_idx = cvrt_dscd_idx,
    kept_cols = cvrt_kept_cols, kept_rows = cvrt_kept_rows
  )

  # Genotypes
  gntp_ipt <- optlst$genotypes
  gntp_ptn <- optlst$gntp_ptn
  gntp_zip <- optlst$gntp_zip
  gntp_idx_col <- optlst$gntp_idx_col
  gntp_dscd_idx <- as.logical(optlst$gntp_dscd_idx)
  gntp_dtfm_trp <- as.logical(optlst$gntp_dtfm_trp)
  gntp_kept_rows <- hdmtopt(optlst$gntp_kept_rows)
  gntp_dscd_rows <- hdmtopt(optlst$gntp_dscd_rows)
  gntp_kept_cols <- hdmtopt(optlst$gntp_kept_cols)
  gntp_dscd_cols <- c(hdmtopt(optlst$gntp_dscd_cols), glbl_blck_lst)

  gntp_dtfm <- ldgntp(
    gntp_ipt = gntp_ipt, gntp_ptn = gntp_ptn, gntp_zip = gntp_zip,
    gntp_idx_col = gntp_idx_col, gntp_dscd_idx = gntp_dscd_idx,
    dscd_cols = gntp_dscd_cols
  )

  qtl_dtfm <- mkme(pntp_dtfm, gntp_dtfm, cvrt_dtfm)$all$eqtls

  # Save QTL mapping result into temporary file
  opt_tmp_file <- file.path(optlst$tmp_dir, "qtlmp_tmp.csv")
  fwrite(qtl_dtfm, opt_tmp_file, showProgress = FALSE, verbose = FALSE)
  cat("[INFO]  Please check", opt_tmp_file, " for temporary QTL list.\n")
}


# sub-command `qtlrpt`
qtlrpt <- function(optlst) {
  tmp_dir <- optlst$tmp_dir
  opt_dir <- optlst$opt_dir
  glbl_blck_lst <- optlst$glbl_blck_lst

  # Read temp-file from last step
  ipt_qtl_flnm <- file.path(tmp_dir, "qtlmp_tmp.csv")
  if (file.exists(ipt_qtl_flnm))
    qtl_dtfm <- pprcs(ipt_qtl_flnm, trps = FALSE, as_idx = NULL)
  else
    stop("Failed to find temporary file ", ipt_qtl_flnm)

  # QTLs with genotype information
  gntp_ifmt_flnm <- optlst$gntp_ifmt_flnm
  gntp_ifmt_zip <- optlst$gntp_ifmt_zip

  qtl_dtfm <- itsct_qtlmt_ifmt(qtl_dtfm, gntp_ifmt_flnm) # overwrite qtl_dtfm

  ## Save QTLs for each measurements and draw Manhattan plots
  qtlrpt_mhtnplt_cols <- hdmtopt(optlst$qtlrpt_mhtnplt_cols)
  qtlrpt_mhtnplt_svfmt <- optlst$qtlrpt_mhtnplt_svfmt
  qtlrpt_svrpt <- as.logical(optlst$qtlrpt_svrpt)
  qtlrpt_svpvt <- as.numeric(optlst$qtlrpt_svpvt)
  pvalue <- qtlrpt_mhtnplt_cols[4]

  msmnts <- unique(qtl_dtfm$gene) # msmnts: measurements a.k.a traits or phenotypes
  for (msmnt in msmnts[1:3]) {
    msmnt_dtfm <- qtl_dtfm[which(qtl_dtfm$gene == msmnt), ]

    opt_prfx <- file.path(optlst$opt_mhtn_dir, msmnt)
    plt_mht(msmnt_dtfm, opt_prfx = opt_prfx, use_cols = qtlrpt_mhtnplt_cols)
    cat("[INFO] Check", paste0(optlst$opt_mhtn_dir, opt_prfx), "for Manhattan plot\n")

    # whether save QTLs
    if (qtlrpt_svrpt) {
      opt_flnm <- file.path(optlst$opt_qtl_dir, paste0(msmnt, ".csv"))
      msmnt_dtfm <- msmnt_dtfm[which(msmnt_dtfm[pvalue] < qtlrpt_svpvt), ]
      fwrite(msmnt_dtfm, opt_flnm, row.names = FALSE, quote = FALSE, sep = ",")
      cat("[INFO] Check", opt_flnm, "for QTL information\n")
    }
  }

  # Fetch top SNPs per chrom
  qtlrpt_tspc_pvt <- as.numeric(optlst$qtlrpt_tspc_pvt) # Top SNP per chrom p-value threshold
  qtlrpt_tspc_cols <- hdmtopt(optlst$qtlrpt_tspc_cols) # Cols will be used from the qtl_dtfm

  tspc <- ftspc(qtl_dtfm = qtl_dtfm, pvt = qtlrpt_tspc_pvt, qtl_col = qtlrpt_tspc_cols)

  if (dim(tspc)[1] == 0)
    return(cat("[WARN]  tspc(dataframe top snp per chrom) is empty. Exit ...\n"))

  #  Plot genotype level for top SNPs
  ## Fetch phenotypes by read temp-files. A dataframe with column as subjects, rows as measurements
  ipt_tmp_flnm <- file.path(tmp_dir, "trfm_tmp.csv")
  pntp_idx_col <- optlst$pntp_idx_col
  if (!file.exists(ipt_tmp_flnm))
    stop("Failed to find temporary file ", ipt_tmp_flnm)
  pntp_dtfm <- pprcs(ipt_tmp_flnm, dscd_cols = pntp_idx_col, trps = TRUE)
  use_pntp <- pntp_dtfm[which(base::rownames(pntp_dtfm) %in% unique(tspc$gene)), ]

  # Dump VCF records for the top SNPs per chromosome to drive
  # FIXME: need a lot dependencies, could be not efficient.
  qtlrpt_vcf_db <- optlst$qtlrpt_vcf_db
  qtlrpt_vcf_rfg <- optlst$qtlrpt_vcf_rfg
  qtlrpt_vcf_cols <- hdmtopt(optlst$qtlrpt_vcf_cols) # FIXME: Perhaps useless

  dmptsvcf(tspc, vcfdb = qtlrpt_vcf_db, vcfrfg = qtlrpt_vcf_rfg, opt_dir = opt_dir)

  ## working on top snps
  ## Fetch dosages of genotype
  gntp_ipt <- optlst$genotypes
  gntp_ptn <- optlst$gntp_ptn
  gntp_zip <- as.logical(optlst$gntp_zip)
  gntp_idx_col <- optlst$gntp_idx_col
  gntp_dscd_idx <- as.logical(optlst$gntp_dscd_idx)
  gntp_dscd_cols <- c(hdmtopt(optlst$gntp_dscd_cols), glbl_blck_lst)

  gntp_dtfm <- ldgntp(
    gntp_ipt = gntp_ipt, gntp_ptn = gntp_ptn, gntp_zip = gntp_zip,
    gntp_idx_col = gntp_idx_col, gntp_dscd_idx = gntp_dscd_idx,
    dscd_cols = gntp_dscd_cols
  )
  use_gntp <- gntp_dtfm[which(base::rownames(gntp_dtfm) %in% unique(tspc$snps)), ]

  # dsg_dtfm, to draw dosage ~ genotype boxplot by plt_dsg()
  qtlrpt_dsgplt_svfmt <- optlst$qtlrpt_dsgplt_svfmt
  qtlrpt_dsgplt_optdir <- optlst$opt_dsg_dir

  dsg_dtfm <- apply(tspc, 1, enc_gntp, use_gntp, use_pntp)
  plt_dsg(dsg_dtfm, svfmt = qtlrpt_dsgplt_svfmt, opt_dir = qtlrpt_dsgplt_optdir)

  cat("[INFO]  Finished dump of reports")
}
#
# Sub-commands ends
#


#
# Main function
#
main <- function() {
  tryCatch({
    # Command line arguments
    agmntlst <- prsarg()

    suppressWarnings(require(dplyr, quietly = TRUE, warn.conflicts = FALSE))
    suppressWarnings(require(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
    suppressWarnings(require(stringi, quietly = TRUE, warn.conflicts = FALSE))
    suppressWarnings(require(reshape2, quietly = TRUE, warn.conflicts = FALSE))
    suppressWarnings(require(data.table, quietly = TRUE, warn.conflicts = FALSE))

    opts <- agmntlst$options
    cfglst <- prcs_cfg(opts$config_file)
    optlst <- mgcfgcmd(agmntlst, cfglst)

    tmp_dir <- optlst$tmp_dir
    opt_dir <- optlst$opt_dir

    stwkspc(tmp_dir, opt_dir)

    # Add internal options to option list
    opt_lg_dir <- file.path(opt_dir, "logs")
    opt_rpt_dir <- file.path(opt_dir, "reports")
    optlst$opt_lg_dir <- opt_lg_dir
    optlst$opt_rpt_dir <- opt_rpt_dir
    optlst$opt_dsg_dir <- file.path(opt_rpt_dir, "GenotypeLevelPlot")
    optlst$opt_qtl_dir <- file.path(opt_rpt_dir, "QtlInformation")
    optlst$opt_mhtn_dir <- file.path(opt_rpt_dir, "ManhattanPlot")
    optlst$opt_pprc_dir <- file.path(opt_rpt_dir, "Preprocessing")

    # Sub-command
    subcmd <- optlst$subcmd

    # Err files
    # errobj <- file.path(opt_lg_dir, paste0(subcmd, ".err"))
    # sink(errobj, append=TRUE, type="message")

    # Log files
    # logobj <- file.path(opt_lg_dir, paste0(subcmd, ".log"))
    # sink(logobj, append=TRUE, type="output")

    # Global black list
    optlst$glbl_blck_lst <- hdmtopt(optlst$glbl_blck_lst)

    ## Check the normality of raw data
    if (subcmd == "cknmlt") {
      cknmlt(optlst)
    } else if (subcmd == "trfm") {
      trfm(optlst)
    } else if (subcmd == "qtlmp") {
      suppressWarnings(require(MatrixEQTL, quietly = TRUE, warn.conflicts = FALSE))
      qtlmp(optlst)
    } else if (subcmd == "qtlrpt") {
      suppressWarnings(require(qqman, quietly = TRUE, warn.conflicts = FALSE))
      suppressWarnings(require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE))
      suppressWarnings(require(VariantAnnotation, quietly = TRUE, warn.conflicts = FALSE))
      qtlrpt(optlst)
    } else {
      stop("Unknow error while processing sub-command: ", subcmd)
    }

    cat("[INFO]  Success!!!\n\n------------------------------------\n")
  }, error = function(e) {
    cat("[ERROR] ", e$message, "\n")
    print(e)
  })
}

#
# Main function
#
main()


#
# [TEST]
#
test_fetch_vcf <- function() {
  tspc <- data.frame(
    snp = c(
      "rs200709595", "rs190867312", "rs201488854", "rs183605470", "rs188652299",
      "rs191297051", "rs183209871", "rs187855973", "rs191379015", "rs201721682",
      "rs183470350", "rs187802690"
    ),
    gene = c(paste0("gene", as.character(1:12))),
    pvalue = c(rep(5e-8, 12)),
    chr = as.character(rep(1, 12)),
    pos = c(84030, 84079, 84133, 84139, 84156, 84244, 84295, 84346, 84453, 84683, 84705, 85063),
    ref = c("G", "T", "A", "A", "A", "A", "G", "T", "C", "A", "T", "T"),
    alt = c("A", "C", "T", "T", "C", "C", "A", "C", "G", "G", "G", "C")
  )

  vcfprm <- mkvcfprm(rgstr = apply(tspc, 1, function(x) {
    paste0(x["chr"], ":", as.character(x["pos"]))
  }))

  vcfrfg <- "hg19"
  vcfdb <- "test.vcf.gz"
  tspc_vcf <- readVcf(TabixFile(vcfdb), vcfrfg, vcfprm)

  # Remove elements of INFO
  tspc_info <- info(tspc_vcf)
  tspc_info <- tspc_info[, which(!colnames(tspc_info) %in% c("ASS", "ASP"))]
  info(tspc_vcf) <- tspc_info

  tspc_header_info <- info(header(tspc_vcf))

  tspc_header_info <- tspc_header_info[which(!rownames(tspc_header_info) %in% c("ASP", "ASS")), ]
  info(header(tspc_vcf)) <- tspc_header_info
  info(header(tspc_vcf))
}

# test_fetch_vcf()
