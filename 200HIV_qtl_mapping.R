#!/usr/bin/env Rscript
# -*- utf-8 -*-

# @version: 0.0.1
# @author: Zhenhua Zhang
# @eamil: zhenhua.zhang217@gmail.com
# @refer: qtl_mapping.R (Xiaojing Chu)
# @date: 4th June, 2019

library(ggplot2)
library(reshape2)
library(outliers)
library(data.table)
library(MatrixEQTL)
library(qqman)

rm(list=ls())

# Preprocess
pprcs <- function(cvrt_file, dscd=NULL, blck=NULL, trps=TRUE, as_idx="id",
                  zipped=FALSE) {
  
  if(zipped) {
    dtfm <- fread(
      cmd=paste0("zcat < ", cvrt_file), data.table=FALSE, header=TRUE,
      stringsAsFactors=FALSE, verbose=FALSE, showProgress=FALSE
    )
  } else {
    dtfm <- fread(
      cvrt_file, data.table=FALSE, header=TRUE, stringsAsFactors=FALSE,
      verbose=FALSE, showProgress=FALSE
    )
  }
  base::row.names(dtfm) <- dtfm[, as_idx]
  
  col_names <- base::colnames(dtfm)
  row_names <- base::rownames(dtfm)
  dtfm <- dtfm[, sort(col_names)]
  
  if (is.null(dscd)) {
    print("[INFO]: discarding list is empty ...")
  } else {
    kept_cols <- col_names[which( !col_names %in% dscd)]
    dtfm <- dtfm[, kept_cols]
  }
  
  if (is.null(blck)) {
    print("[INFO]: blck list in empty ...")
  } else {
    kept_rows <- row_names[which(!row_names %in% blck)]
    dtfm <- dtfm[kept_rows, ]
  }
  
  dtfm <- dtfm[, !(base::colnames(dtfm) %in% c(as_idx))]
  
  if (trps == TRUE) return(t(dtfm))
  else return(dtfm)
}

# Inverse rank normalization
ivsrk <- function(x, lude.mode =TRUE){
  if(lude.mode){
    res <- rank(x)
    res <- qnorm((0.5+res-1)/(length(res)))
    res <- (res*sd(x))+ mean(x)
  }
  else{
    res <- rank(x)
    res <- qnorm(res/(length(res)+0.5))
    res <- (res*sd(x))+ mean(x)
  }  
  return(res)
}

# Transform non-normal distributed phnotype column into
trfm_nmlt <- function(dtfm, dscd_cols=NULL, trfm_cols=NULL, trfm_mthd="log2",
                      trps=FALSE) {

  row_names <- base::rownames(dtfm)
  cols_names <- sort(base::colnames(dtfm))
  dtfm <- dtfm[cols_names]

  if(!is.null(dscd_cols)) {
    cat("[WARN]  will discard:", paste(dscd_cols, sep=", "), "\n")
    dscd_idx <- which(cols_names %in% dscd_cols)
    dtfm <- dtfm[, -dscd_idx]
    cols_names <- base::colnames(dtfm)
  }

  if(!is.null(trfm_cols)) {
    if(!is.null(dscd_cols)){
      ovlp_cols <- trfm_cols[which(dscd_cols %in% trfm_cols)]
      if(length(ovlp_cols) != 0) {
        stop("[ERROR]  overlapped between dscd_cols and trfm_cols: ", ovlp_cols)
      }
    }

    if (is.null(trfm_mthd)) {
      stop("[ERROR]  When trfm_cols is not NULL, trfm_mthd is required.")
    }
    
    dtfm[, trfm_cols] <- sapply(dtfm[, trfm_cols], trfm_mthd)
    base::colnames(dtfm) <- sapply(cols_names, function(x) {ifelse(x%in%trfm_cols, paste(x, trfm_mthd, sep="_"), x)})
  }
  
  if(trps) {
    cat("[INFO]  Will transpose dataframe: dtfrm")
    dtfm <- t(dtfm)
  }
  
  dtfm <- as.data.frame(scale(dtfm))
  return(dtfm)
}


# Check nomality
chck_nmlt <- function(dtfm, hst_opt="Normality_check.pdf") {
  row_names <- base::rownames(dtfm)
  cols_names <- sort(base::colnames(dtfm))
  dtfm <- dtfm[cols_names]
  dtfm <- as.data.frame(scale(dtfm, center=FALSE))

  exp_mlt <- reshape2::melt(dtfm)

  nvars <- length(unique(exp_mlt["variable"])[[1]])
  ncol <- as.integer(sqrt(nvars))
  p <- ggplot(data=exp_mlt)
  p <- p + geom_histogram(aes(x=value))
  p <- p + facet_wrap(.~variable, ncol=ncol, scales="free")
  ggsave(hst_opt, width=ncol*3, height=ncol*3)
}


# Draw PCA
plt_pca <- function(dtfm, opt_flnm="PCA.pdf") {
  if(!is.data.frame(dtfm)) {
    stop("[ERROR]  `dtfm` should be a dataframe ...")
  }
  
  dtfm[is.na(dtfm)] <- mean(base::rowMeans(dtfm), na.rm=TRUE)
  pca <- prcomp(dtfm)
  pcax <- pca$x
  pcac <- pca$center
  
  p <- ggplot() + theme_bw()
  p <- p + geom_point(data=as.data.frame(pcax), mapping=aes(x=PC1, y=PC2))
  
  ggsave(opt_flnm, width=10, height=10, units="in")
}


# A PCA analysis to remove outliers (out of 3*sd or 4*sd)
prcs_otls <- function(dtfm, n_sd=3, trps=FALSE) {
  row_names <- rownames(dtfm)
  masked <- sapply(
    dtfm, function(x) {
      sd=sd(x, na.rm=TRUE)
      mu=mean(x, na.rm=TRUE)
      flt <- which((x<=mu-n_sd*sd) | (x>=mu+n_sd*sd))
      x[flt] <- NA
      return(x)
    }
  )
  masked <- as.data.frame(masked, row.names=row_names)
  cat("[INFO]  Number of masked subjects in each measurments: \n")
  print(sapply(masked, function(x) {sum(is.na(x))}))
  if(trps) {
    return(t(masked))
  }
  return(masked)
}


crt_slcdt <- function(ipt, dlmt=',', omchr="NA", skprw=0, skpcl=0,
                      slice_size=2000) {
  slcdt <- SlicedData$new()
  slcdt$fileDelimiter <- dlmt
  slcdt$fileOmitCharacters <- omchr
  slcdt$fileSkipRows <- skprw
  slcdt$fileSkipColumns <- skpcl
  slcdt$fileSliceSize <- slice_size
  
  if (is.null(ipt)){
    stop("[ERROR]  input `ipt` should not be NULL ...")
  } else if(is.matrix(ipt)) {
    slcdt$CreateFromMatrix(ipt)
  } else if(is.data.frame(ipt)) {
    slcdt$CreateFromMatrix(as.matrix(ipt))
  } else {
    slcdt$LoadFile(ipt)
  }

  return(slcdt)
}

# Main method to do QTL mapping
mkme <- function(pntp, gntp, cvrt=NULL, otpt_file="qtl_output.txt",
                 pvOutputThreshold=1e-2) {
  
  useModel = modelLINEAR
  errorCovariance <- numeric()
  
  snps <- crt_slcdt(pntp)
  gene <- crt_slcdt(gntp, dlmt="\t", skprw=1, skpcl=1)
  cvrt <- crt_slcdt(cvrt)

  me <- Matrix_eQTL_engine(
    snps=snps,
    gene=gene,
    cvrt=cvrt,
    output_file_name=otpt_file,
    pvOutputThreshold=pvOutputThreshold,
    useModel=useModel,
    errorCovariance=errorCovariance,
    verbose=TRUE,
    pvalue.hist=TRUE,
    min.pv.by.genesnp=FALSE,
    noFDRsaveMemory=FALSE
  )
  
  return(me)
}


# Covariates
blck_lst <- c("X1012",  "X1053", "X1109",  "X1128", "X1129",  "X1142", "X1150", "X1193", "X1009", "X1101", "X1125", "X1126")
cvrt_file <- "/home/umcg-zzhang/Documents/projects/200HIV/inputs/datasets/metaData_2019-03-20.csv"
cvrt_dtfm <- pprcs(cvrt_file, blck=blck_lst)

# Phenotypes
pntp_file <- "/home/umcg-zzhang/Documents/projects/200HIV/inputs/datasets/M_chemokineCytokineData200Hiv_boxCorData_log10_2019-03-20.csv"
exp_dtfm <- pprcs(pntp_file, dscd=c("cohort"), blck=blck_lst, trps=FALSE)

chck_nmlt(exp_dtfm, hst_opt="Normality_check_raw.pdf")

trans <- colnames(exp_dtfm) # c("AAT", "Beta_TG", "CCL5", "CRP", "CXCL7", "d.dimer", "IFABP", "IL-10", "IL-1A", "IL-1RA", "IL6_LOD", "Leptin",  "PF4", "Resistin", "sCD14",  "sCD163", "TAT2", "vWf")
discard <- NULL # c("AAT", "IL-1A", "IL1b_LLOD")

exp_nmlz_dtfm <- trfm_nmlt(exp_dtfm, dscd_cols=discard, trfm_cols=trans)
chck_nmlt(exp_nmlz_dtfm, hst_opt="Normality_check_trans.pdf")

plt_pca(exp_nmlz_dtfm)
exp_nmlz_rmotl_dtfm <- prcs_otls(exp_nmlz_dtfm, trps=TRUE, n_sd=4)

# QTL mapping
gntp_drct <- "/home/umcg-zzhang/Documents/projects/200HIV/inputs/dosage/200HIV_dosages"
for (gntp_file in list.files(gntp_drct, pattern="*.gz")) {
  qtl_opt <- sub(".gz", "_qtl_mapping.txt", gntp_file)
  gntp_file <- file.path(gntp_drct, gntp_file)
  
  gntp_dtfm <- pprcs(gntp_file, dscd=blck_lst, trps=FALSE, zipped=TRUE)  # Load genotypes
  tmp <- mkme(exp_nmlz_rmotl_dtfm, gntp_dtfm, cvrt_dtfm, otpt_file=qtl_opt) # QTL mapping
}

