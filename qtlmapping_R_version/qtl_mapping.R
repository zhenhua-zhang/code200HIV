#!/usr/bin/env Rscript
#
# File Name  : qtl_mapping.R
# Author     : Xiaojing Chu
# E-mail     : chuxj1128@gmail.com
# Created on : Thu 30 May 2019 12:29:22 PM CEST
# Version    : v0.0.1
# License    : 
#

library(MatrixEQTL);

useModel = modelLINEAR;
SNP_file_name="snp_for_list4_mapping.txt";
expression_file_name="exp_for_list4_mapping.txt";
covariates_file_name="cov_for_list4_mapping.txt";
output_file_name="qtl_out_list4.txt";

pvOutputThreshold = 1e-2;
errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

