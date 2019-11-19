#!/usr/bin/env Rscript

#
## mewrapper.R
#

require(optparse, quietly=T, warn.conflicts=F)
require(data.table, quietly=T, warn.conflicts=F)
require(MatrixEQTL, quietly=T, warn.conflicts=F)


# Function to create SliceData from given input file
mk_sliceddata <- function(ipt, dlm="\t", omch="NA", skrw=1, skcl=1, ssize=2000)
{
  # @description Create SliceData object from given matrix, data.frame or file
  # @param ipt string/matrix/data.frame. Input
  # @param dlm string. Delimiter in given input. Default: "\t"
  # @param omch string. Omitted the char. Default: "NA"
  # @param skrw integer. Skipping rows. Default: 0
  # @param skcl integer. Skipping columns. Default: 0
  # @param ssize integer. The size of slice. Default: 2000

  slcdt <- SlicedData$new()
  slcdt$fileDelimiter <- dlm
  slcdt$fileOmitCharacters <- omch
  slcdt$fileSkipRows <- skrw
  slcdt$fileSkipColumns <- skcl
  slcdt$fileSliceSize <- ssize

  if (is.null(ipt))
    stop("Input `ipt` should not be NULL!!")
  else if (is.matrix(ipt))
    slcdt$CreateFromMatrix(ipt)
  else if (is.data.frame(ipt))
    slcdt$CreateFromMatrix(as.matrix(ipt))
  else
    slcdt$LoadFile(ipt)

  return(slcdt)
}


mk_matrixeqtl <- function(phtp, gntp, cvrt=NULL, opt_file=NULL, pv_thr=1)
{
  # @description Run Matrix_eQTL_engine
  # @param phtp data.frame. A data.frame including phenotype measurements.
  # @param gntp data.frame. A data.frame including genotype measurements.
  # @param cvrt data.frame. A data.frame including covariates. Default: NULL
  # @param opt_file string. A string indicate the output file. Default: NULL
  # @param pv_thr float. P-value output threshold. Default: 5e-2

  useModel <- modelLINEAR  # For qtl mapping it should be linear model
  errorCovariance <- numeric()

  phtp <- mk_sliceddata(phtp)
  snps <- mk_sliceddata(gntp)
  cvrt <- mk_sliceddata(cvrt)

  me <- Matrix_eQTL_engine(
    snps=snps, gene=phtp, cvrt=cvrt,
	output_file_name=opt_file,
	pvOutputThreshold=pv_thr,
	useModel=useModel,
	errorCovariance=errorCovariance,
	verbose=F,
	pvalue.hist=T,
	noFDRsaveMemory=F,
	min.pv.by.genesnp=F
  )

  return(me)
}


# Read file using fread() from data.table. Dealing with zipped files
fr_wrapper <- function(flnm, zipped=F, data.table=F, header=T, stringsAsFactors=F, showProgress=F, verbose=F)
{
	print(flnm)
  if (zipped) {
	return(fread(cmd=paste0("zcat < ", flnm), data.table=data.table, header=header, stringsAsFactors=stringsAsFactors, verbose=verbose, showProgress=showProgress))
  } else {
	return(fread(flnm, data.table=data.table, header=header, stringsAsFactors=stringsAsFactors, verbose=verbose, showProgress=showProgress))
  }
}


# Parsing command line options
get_arg <- function()
{
	parser <- OptionParser()

	parser <- add_option(parser, c("-p", "--phenotypes"), dest="phenotype", type="character", help="Phenotypes")
	parser <- add_option(parser, c("-c", "--covariates"), dest="covariate", type="character", help="Covariates")
	parser <- add_option(parser, c("-g", "--genotypes"), dest="genotype", type="character", help="Genotypes")
	parser <- add_option(parser, c("-o", "--output"), dest="output", type="character", help="Output")

	agmnts <- parse_args(parser)
}


main <- function()
{
	opt_list = get_arg()
	gntp_inpt = opt_list$genotype
	pntp_inpt = opt_list$phenotype
	cvrt_inpt = opt_list$covariate
	output = opt_list$output

	me <- mk_matrixeqtl(pntp_inpt, gntp_inpt, cvrt_inpt)
	qtls_dtfm <- me$all$eqtls
	fwrite(qtls_dtfm, output, showProgress=F, verbose=F)
}

main()
