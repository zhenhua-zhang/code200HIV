#!/usr/bin/env Rscript

# @version: 0.0.1
# @author: Zhenhua Zhang
# @eamil: zhenhua.zhang217@gmail.com
# @refer: qtl_mapping.R (Xiaojing Chu)
# @date: 4th June, 2019


# TODO: 1. When doing inverse-rank transform, use mean to replace NA, which is
#           dangerous. Need a better way to deal with NAs
#       2. Make the script rescuable at broken-points with additional command
#           arguments


# NOTE: 1. Settings in config file will be overrided by command line arguments
#       2. Some settings is only modifiable by configuration files
#		3. Keep the row as measurements, column as subjects


library(qqman)
library(dplyr)
library(ggplot2)
library(methods)
library(stringi)
library(optparse)
library(reshape2)
library(data.table)
library(MatrixEQTL)

rm(list=ls())


# Processing the configuration file
prcs_cfg <- function(iptfl) {
	cfglst <- list()
	cfgs <- paste(readLines(iptfl), sep="\n")
	for(line in cfgs) {
		if(!startsWith(line, "#") && !stri_isempty(line)){
			plst <- strsplit(line, "=")
			if(length(plst) != 1) {
				stop("[ERROR]  The length of plst should be 1")
			}
			key <- plst[[1]][1]
			val <- plst[[1]][2]
			cfglst[[key]] <- val
		}
	}
	return(cfglst)
}


# Read file using fread() from data.table. Dealing with zipped files
smt_fread <- function(flnm, zipped=FALSE, data.table=FALSE, header=TRUE,
					  stringsAsFactors=FALSE, verbose=FALSE,
					  showProgress=FALSE) {
	#	cat("[INFO]  Reading\n          ", flnm, "\n")
	if(zipped) {
		return(fread(cmd=paste0("zcat < ", flnm), data.table=data.table,
					 header=header, stringsAsFactors=stringsAsFactors,
					 verbose=verbose, showProgress=showProgress 
			 ))
	} else {
		return(fread(flnm, data.table=data.table, header=header,
					 stringsAsFactors=stringsAsFactors, verbose=verbose,
					 showProgress=showProgress 
			 ))
	}
}

# Preprocess
# NOTE: 1. dscd working on columns, blck working on rows
pprcs <- function(flnm, dscd=NULL, blck=NULL, trps=TRUE, as_idx="id",
				  zipped=FALSE) {
	cat("[INFO]  Processing\n          ", flnm, "\n")

	dtfm <- smt_fread(flnm, zipped=zipped)

	if(!is.null(as_idx)){ base::row.names(dtfm) <- dtfm[, as_idx] }

	col_names <- base::colnames(dtfm)
	row_names <- base::rownames(dtfm)
	dtfm <- dtfm[, sort(col_names)]
	col_names <- base::colnames(dtfm)

	if (is.null(dscd)) {
		cat("[INFO]  Discarding list is empty ...\n")
	} else {
		cat("[WARN]  Will discard column(s)\n          ", dscd, "...\n")
		kept_cols <- col_names[which( !col_names %in% dscd)]
		dtfm <- dtfm[, kept_cols]
	}

	if (is.null(blck)) {
		cat("[INFO]  Black list is empty ...\n")
	} else {
		cat("[WARN]  Will black-out rows(s)\n          ", blck, "...\n")
		kept_rows <- row_names[which(!row_names %in% blck)]
		dtfm <- dtfm[kept_rows, ]
	}

	# dtfm <- dtfm[, !(base::colnames(dtfm) %in% c(as_idx))]

	if(trps) {
		cat("[WARN]  Will transpose dataframe: dtfm\n")
		return(as.data.frame(t(dtfm), stringsAsFactors=FALSE))
	} else {
		return(dtfm)
	}
}


###############################################################################
##       Check the normality of each measuements                             ##
###############################################################################
chck_nmlt <- function(dtfm, hst_opt="Normality_check.pdf", bins=30) {
	cat("[INFO]  Check the normality...\n")

	if(!is.data.frame(dtfm)) { stop("[ERROR]  `dtfm` should be a dataframe ...") }

	col_names <- sort(base::colnames(dtfm))
	dtfm <- dtfm[col_names]
	exp_mlt <- reshape2::melt(dtfm, id.vars='id')

	nvars <- length(unique(exp_mlt["variable"])[[1]])
	ncol <- as.integer(sqrt(nvars))

	p <- ggplot(data=exp_mlt)
	p <- p + geom_histogram(aes(x=value), na.rm=TRUE, bins=bins)
	p <- p + facet_wrap(.~variable, ncol=ncol, scales="free")

	ggsave(hst_opt, width=ncol*3, height=ncol*3)
}


###############################################################################
##       Transform                                                           ##
###############################################################################
# Inverse rank normalization
ivsrk <- function(x, lude.mode=TRUE){
	is_na <- is.na(x)
	x[is_na] <- mean(x, na.rm=TRUE)
	if(lude.mode){
		res <- base::rank(x)
		res <- qnorm((0.5+res-1)/(length(res)))
		res <- (res*sd(x)) + mean(x)
	}
	else{
		res <- base::rank(x)
		res <- qnorm(res/(length(res)+0.5))
		res <- (res*sd(x))+ mean(x)
	}  
	res[is_na] <- NA

	return(res)
}

# Trans form log10 values into log2 values
lg10tolg2 <- function(x) { return(x * log2(10)) }
lg2tolg10 <- function(x) { return(x * log10(2)) }


# NOTE: Make sure only the id column is characters, while the rest are numeric,
#       otherwise, the transpose would fail
# Transform non-normal distributed phnotype column into
trfm_nmlt <- function(dtfm, dscd_cols=NULL, trfm_cols=NULL, trfm_mthd="log2",
					  trps=FALSE, idcol="id") {
	cat("[INFO]  Transform data by", trfm_mthd, "\n")
	row_names <- base::rownames(dtfm)
	dtfm[idcol] <- NULL
	col_names <- sort(base::colnames(dtfm))
	dtfm <- dtfm[col_names]

	if(!is.null(dscd_cols)) {
		cat("[WARN]  Will discard: \n         ", dscd_cols, "\n")
		dscd_idx <- which(col_names %in% dscd_cols)
		dtfm <- dtfm[, -dscd_idx]
		col_names <- base::colnames(dtfm)
	}

	if(!is.null(trfm_cols) && !is.na(trfm_cols)) {
		if(!is.null(dscd_cols)){
			ovlp_cols <- trfm_cols[which(dscd_cols %in% trfm_cols)]
			if(length(ovlp_cols) != 0) {
				stop("[ERROR]  Overlapped between dscd_cols and trfm_cols \n",
					 ovlp_cols)
			}
		}

		if (is.null(trfm_mthd)) {
			stop("[ERROR]  When trfm_cols is not NULL, trfm_mthd is required.")
		}

		if(length(trfm_cols)==1){
			dtfm[, trfm_cols] <- get(trfm_mthd)(dtfm[, trfm_cols])
		} else {
			dtfm[, trfm_cols] <- sapply(dtfm[, trfm_cols], trfm_mthd)
		}
		base::colnames(dtfm) <- sapply(col_names, 
			function(x) {ifelse(x%in%trfm_cols, paste(x, trfm_mthd, sep="_"), x)}
		)
	} else {
		cat("[WARN]  trfm_cols is either NULL or NA\n")
	}

	if(trps) {
		cat("[INFO]  Will transpose dataframe: dtfm \n")
		dtfm <- t(dtfm)
	}

	dtfm <- as.data.frame(scale(dtfm), stringsAsFactors=FALSE)
	dtfm[idcol] <- row_names
	return(dtfm)
}


# Draw PCA
pltpca <- function(dtfm, opt_flnm="Outlier_check_PCA.pdf", idcol="id") {
	cat("[INFO]  Plotting PCA... \n")
	if(!is.data.frame(dtfm)) {
		stop("[ERROR]  `dtfm` should be a dataframe ...")
	}

	dtfm[idcol] <- NULL
	dtfm[is.na(dtfm)] <- mean(base::rowMeans(dtfm), na.rm=TRUE)
	pca <- prcomp(dtfm)
	pcax <- pca$x
	pcac <- pca$center

	p <- ggplot() + theme_bw()
	p <- p + geom_point(data=as.data.frame(pcax, stringsAsFactors=FALSE), mapping=aes(x=PC1, y=PC2))

	ggsave(opt_flnm, width=10, height=10, units="in")
}


# A PCA analysis to remove outliers (out of 3*sd or 4*sd)
prcs_otls <- function(dtfm, n_sd=3, trps=FALSE, idcol="id") {
	cat("[INFO]  Removing outliers...\n")
	row_names <- base::rownames(dtfm)
	dtfm[idcol] <- NULL

	masked <- sapply(dtfm,
									 function(x) {
										 sd=sd(x, na.rm=TRUE)
										 mu=mean(x, na.rm=TRUE)
										 flt <- which((x<=mu-n_sd*sd) | (x>=mu+n_sd*sd))
										 x[flt] <- NA
										 return(x)
									 })

	masked <- as.data.frame(masked, row.names=row_names, stringsAsFactors=FALSE)
	cat("[INFO]  Number of masked subjects in each measurments:\n")
	nacnt <- sapply(masked, function(x) {sum(is.na(x))})

	col_names <- base::colnames(masked)
	masked[idcol] <- row_names
	masked <- masked[c(idcol, col_names)]

	if(trps) {
		cat("[INFO]  Will transpose dataframe: masked\n")
		return(t(masked))
	} else {
		return(masked) 
	}
}


crt_slcdt <- function(ipt, dlmt=",", omchr="NA", skprw=0, skpcl=0,
											slice_size=2000) {
	#	cat("[INFO]  Creating SlicedData ...\n")
	slcdt <- SlicedData$new()
	slcdt$fileDelimiter <- dlmt
	slcdt$fileOmitCharacters <- omchr
	slcdt$fileSkipRows <- skprw
	slcdt$fileSkipColumns <- skpcl
	slcdt$fileSliceSize <- slice_size

	if (is.null(ipt)){
		stop("[ERROR]  input `ipt` should not be NULL!!")
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
mkme <- function(pntp, gntp, cvrt=NULL, otpt_file=NULL, pv_opt_thr=5e-2) {
	cat("[INFO]  Mapping QTL ...\n")
	useModel = modelLINEAR
	errorCovariance <- numeric()
	print(cvrt[1:5, 1:5])

	gene <- crt_slcdt(pntp)
	snps <- crt_slcdt(gntp, dlmt="\t")
	cvrt <- crt_slcdt(cvrt)

	me <- Matrix_eQTL_engine(
		 snps=snps, gene=gene, cvrt=cvrt,
		 output_file_name=otpt_file, pvOutputThreshold=pv_opt_thr,
		 useModel=useModel, errorCovariance=errorCovariance, verbose=FALSE,
		 pvalue.hist=TRUE, min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE
	)

	return(me)
}


# Intersect Matrix eQTL with SNP information
intsct_qtlmt_ifmt <- function(qtls, snp_ifmt, mgby_x="snps", mgby_y="rsID") {
	cat("[INFO]  Intersect QTL matrix with SNP information ...\n")
	if(is.character(qtls)) { qtls <- smt_fread(qtls) }
	if(is.character(snp_ifmt)) { snp_ifmt <- smt_fread(snp_ifmt) }
	return(merge(qtls, snp_ifmt, by.x=mgby_x, by.y=mgby_y))
}


# Draw Manhattan plot and Q-Q plot
pltme <- function(qtls_lst, use_cols=NULL, otpt_prfx="qtl_mapping",
									svfmt="pdf", wdth=20, hght=10, drwzsc=FALSE, drwqq=FALSE) {
	cat("[INFO]  Plotting Manhattan plot and Q-Q plot for", otpt_prfx, "...\n")
	# use_cols should follow the order: snp_id, chromosome, position, p-value
	if(is.null(use_cols)) {
		use_cols <- c("snps", "chr", "pos", "pvalue") 
	} else if(length(use_cols) != 4) {
		stop("[ERROR]  The length of use_cols should be 4 ...\n")
	}

	gwas_rslt <- qtls_lst[, use_cols]
	base::colnames(gwas_rslt) <- c("SNP", "CHR", "BP", "P")
	otpt_flnm <- paste0(otpt_prfx, "_mhtn_pvalue.", svfmt)

	if(svfmt=="png") { png(otpt_flnm, width=wdth, height=hght) }
	else { pdf(otpt_flnm, width=wdth, height=hght) }

	manhattan(gwas_rslt, main="Manhattan Plot (p-value)", ylan="p-value(-log10)",
						suggestiveline=FALSE, annotatePval=5e-8, annotateTop=FALSE
						)
	dev.off()

	if(drwzsc) {
		gwas_rslt <- transform(gwas_rslt, zscore=qnorm(P/2, lower.tail=FALSE))
		otpt_flnm <- paste0(otpt_prfx, "_mhtn_zscore.", svfmt)
		if(svfmt=="png") { png(otpt_flnm, width=wdth, height=hght) }
		else { pdf(otpt_flnm, width=wdth, height=hght) }
		manhattan(gwas_rslt, main="Manhattan Plot (z-score)", ylab="Z-score",
							suggestiveline=FALSE, annotatePval=0.00000005
							)
		dev.off()
	} else {
		cat("[INFO]  Skipping Manhattan plot for z-score ...\n")
	}

	if(drwqq) {
		otpt_flnm <- paste0(otpt_prfx, "_qq.", svfmt)
		if(svfmt=="png") { png(otpt_flnm, width=wdth, height=hght) }
		else { pdf(otpt_flnm, width=wdth, height=hght) }
		qq(gwas_rslt$P, main="Q-Q plot of p-values", pch=18, cex=1.5, las=1)
		dev.off()
	} else {
		cat("[INFO]  Skipping Q-Q plot ...\n")
	}
}


# Filter out check the number of mode
ckmd <- function(dtfm, thrshld=0.5) {
	if(is.data.frame(dtfm)) {
		return(sapply(dtfm, 
									function(x) {
										cnt <- table(x)
										ifelse(max(cnt)/sum(cnt) >= thrshld, TRUE, FALSE)
									}))
	} else if (is.vector(dtfm)) {
		cnt <- table(dtfm)
		return(ifelse(max(cnt)/sum(cnt) >= thrshld, TRUE, FALSE))
	} else {
		stop("[ERROR]  I can only handle vector or data.frame!!!")
	}
}


# Draw dosage boxplot
# use_col should be as the following order: 
# TRY NOT BE COMPLICATED
plt_dsg <- function(qtllst, pntp, gntp, gntp_ptn="*.dosage", gntp_zipped=TRUE,
					svfmt="pdf", pvt=5e-8, qtl_col=NULL, gntp_blck_lst=NULL,
					opt_dir="qtl_mapping_output_dir/plots") {

	# Fetch top snp on each chr, filted by pvalue
	if(is.null(qtl_col)) { 
		qtl_col <- c( "snps", "gene", "pvalue", "chr", "pos", "ref", "alt") 
	}

	qtllst <- qtllst[, qtl_col]
	qtl_col <- c("snps", "gene", "pvalue", "chr", "pos", "ref", "alt") 
	base::colnames(qtllst) <- qtl_col

	qtl_gb <- c("snps", "chr", "gene")

	if(is.vector(qtl_col) && length(qtl_col)==7){
		top_snps_per_chr <- qtllst %>%
			filter(pvalue <= pvt) %>%
			group_by(chr, gene) %>%
			mutate(minpval=min(pvalue)) %>%
			filter(pvalue == minpval) %>%
			ungroup() %>%
			as.data.frame()

		use_msmnt <- unique(top_snps_per_chr$gene)
		use_snps <- unique(top_snps_per_chr$snps)
	} else {
		stop("qtl_col should be a vector with length 7")
	}
	
	# Fetch phenotypes
	if(is.data.frame(pntp)) {  # A dataframe with column as subjects, rows as measurements
		pntp_rwnm <- base::rownames(pntp)
		use_pntp <- pntp[which(pntp_rwnm %in% use_msmnt), ]
	} else {
		stop("pntp should be a dataframe")
	}

	# Fetch genotype
	# Repeated code, should rewrite into a function
	if(file_test("-d", gntp)) {
		if(is.null(gntp_ptn)) { stop("When gntp is a directory, gntp_ptn shouldn't be NULL or empty") }

		flnmlst <- list.files(gntp, pattern=gntp_ptn)
		if(length(flnmlst) == 0) { stop("gntp is a directory, but NO file following pattern (", gntp_ptn, ")") }

		gntp_dtfm <- data.frame()
		for (flnm in flnmlst) {
			gntp_file <- file.path(gntp, flnm)
			tmp_dtfm <- pprcs(gntp_file, dscd=gntp_blck_lst, trps=FALSE)
			gntp_dtfm <- rbind(gntp_dtfm, tmp_dtfm)
		}
	} else if(file_test("-f", gntp)) {
		gntp_dtfm <- pprcs(gntp_file, dscd=gntp_blck_lst, trps=FALSE)
	} else {
		stop("gntp should be a directory, normal file, or zipped file")
	}
	gntp_rwnm <- base::rownames(gntp_dtfm)
	use_gntp <- gntp_dtfm[which(gntp_rwnm %in% use_snps), ]

	# print(rbind(use_gntp, use_pntp))
	# Encode genotype
	tmp_func <- function(x, use_gntp, use_pntp) {
		snpid <- x[1]
		msmnt <- x[2]
		ref <- x[6]
		alt <- x[7]
		gntp_dsg <- use_gntp[snpid, ]
		gntp_code <- sapply(gntp_dsg, function(d) {
				if(length(d) != 1) { return(NA) }
				if(round(d) == 0){ return(paste0(ref, ref)) }
				else if(round(d) == 1){ return(paste0(ref, alt)) }
				else if(round(d) == 2){ return(paste0(alt, alt)) }
				else{ return(NA) }
			}
		)
		pntp_msmnt <- use_pntp[msmnt, ]
		gnpn <- rbind(pntp_msmnt, gntp_dsg, gntp_code)
		base::rownames(gnpn) <- c(msmnt, snpid, "code")
		return(gnpn)
	}

	lgt <- dim(top_snps_per_chr)[1]
	if(lgt != 0) {
		gnpnlst <- apply(top_snps_per_chr, 1, tmp_func, use_gntp, use_pntp)
	} else {
		cat("[WARN]  top_snps_per_chr is empty. Exit ...\n")
		return(NULL)
	}

	for(gnpn in gnpnlst) {
		gnpn_t <- as.data.frame(t(gnpn))

		col_names <- base::colnames(gnpn_t)
		lvl <- col_names[1]
		snp <- col_names[2]
		opt_flnm <- paste0(lvl, "_", snp, ".", svfmt)
		opt_flnm <- file.path(opt_dir, opt_flnm)

		base::colnames(gnpn_t) <- c("lvl", "dsg", "code")

		gnpn_t["lvl"] <- as.numeric(gnpn_t[, "lvl"])
		gnpn_t["dsg"] <- as.numeric(gnpn_t[, "dsg"])

		p <- ggplot(gnpn_t) + theme_bw()
		p <- p + geom_boxplot(aes(x=code, y=lvl, color=code))
		p <- p + geom_point(aes(x=code, y=lvl))
		ggsave(opt_flnm)
	}
}


# Split arguments with coma or semi-colon in to a vector
hndl_mltp_opt <- function(arg) {
	if(!is.character(arg)) {
		stop("Keyword argument arg should be character") 
	}
	if(length(arg) == 0) {
		cat("[WARN]  Keyword argument arg is empty, return \"\"\n")
		return("")
	}
	if(is.na(arg)) { 
		cat("[WARN]  Keyword argument arg is NA, return NULL\n")
		return(NULL) 
	}
	if(stri_detect_fixed(arg, ",")) {
		return(stri_split_fixed(arg, ",")[[1]])
	} else if(stri_detect_fixed(arg, ";")) {
		return(stri_split_fixed(arg, ";")[[1]])
	} else {
		return(arg)
	}
}
# Parsing command line arguments
# parser <- add_option(parser, c("-", "--"), type="character", default="", help="")
# TODO: -g is only for --gui ??
prsarg <- function() {
	tryCatch(
					 {
						 parser <- OptionParser()

						 parser <- add_option(parser, c("-G", "--genotypes"), type="character", help="Genotypes file or directory")
						 parser <- add_option(parser, c("-C", "--covariates"), type="character", help="Covariates")
						 parser <- add_option(parser, c("-P", "--phenotypes"), type="character", help="Phenotypes file")
						 parser <- add_option(parser, c("-S", "--config-file"), type="character", default="config.cfg", help="File including settings, can be overrided by command arguments")
						 parser <- add_option(parser, c("-T", "--tmp-dir"), type="character", default=NULL, help="File including temporary files")
						 parser <- add_option(parser, c("-o", "--output-dir"), type="character", default=NULL, help="Output directory for the QTL mapping results")
						 parser <- add_option(parser, c("-I", "--interactive"), type="store_true", help="Excute the script interactively")
						 parser <- add_option(parser, c("", "--gntp_ptn"), type="character", default=NULL, help="Genptype files pattern")
						 parser <- add_option(parser, c("", "--gntp_zip"), type="store_true", help="The genotype input file is zipped")

						 agmnts <- parse_args(parser, positional_arguments=1, convert_hyphens_to_underscores=TRUE)

						 subcmd <- agmnts$args
						 subcmdlst <- c("cknmlt", "trfm", "qtlmp", "pltmht", "pltdsg")  # pprcs: preprocessing. qltmp: QTL mapping. pstmp: post QLT mapping
						 if(!subcmd %in% subcmdlst) { stop("Unknown subcmd: ", subcmd, ". Options: [cknmlt, trfm, qtlmp, pltmht, pltdsg]\n") }

						 return(agmnts)
					 }, error=function(e) {
						 cat("[ERROR] ", e$message, "\n")
					 }, finally=function(msg) {
						 cat("[INFO]  Finished tryCatch scope", msg, "\n")
					 }
	)
}


# Setup working space
stwkspc <- function(tmp_dir=NULL, opt_dir=NULL, mode="0755") {
	# output
	# |-- logs
	# |-- plots
	# |-- reports
	# |-- README.md

	if(dir.exists(tmp_dir)) { cat("[WARN] ", tmp_dir, "exists\n") }
	else { dir.create(tmp_dir, mode=mode) }

	if(dir.exists(opt_dir)) { cat("[WARN] ", opt_dir, "exists\n") }
	else {
		dir.create(opt_dir, mode=mode)
		dir.create(file.path(opt_dir, "logs"), mode=mode)
		dir.create(file.path(opt_dir, "plots"), mode=mode)
		dir.create(file.path(opt_dir, "reports"), mode=mode)
	}
}


# Merge the settings from configuration file and command line arguments
mgcfgcmd <- function(agmntlst, cfglst) {
	if(!is.list(agmntlst) || !is.list(cfglst)) { stop("Keyword arguments both agmntlst and cfglst should be list") }
	cmdopts <- agmntlst$options
	cmdopts["subcmd"] <- agmntlst$args
	cfglst[names(cmdopts)] <- cmdopts
	return(cfglst)
}


main <- function() {
	# Trycatch clean up tmp_dir
	tryCatch(
		{
			# blck_lst <- c("X1012",  "X1053", "X1109",  "X1128", "X1129",  "X1142", "X1150", "X1193", "X1009", "X1101", "X1125", "X1126")
			# Command line arguments
			agmntlst <- prsarg()
			opts <- agmntlst$options
			cfglst <- prcs_cfg(opts$config_file)
			optlst <- mgcfgcmd(agmntlst, cfglst)

			tmp_dir <- optlst$tmp_dir
			output_dir <- optlst$output_dir

			stwkspc(tmp_dir, output_dir)
			opt_logs_dir <- file.path(output_dir, "logs")
			opt_plots_dir <- file.path(output_dir, "plots")
			opt_reports_dir <- file.path(output_dir, "reports")

			# Positional argument
			subcmd <- optlst$subcmd

			# sink(file(file.path(opt_logs_dir, paste0(subcmd, ".err"))), type="message")
			# sink(file(file.path(opt_logs_dir, paste0(subcmd, ".log"))), type="output")

			glbl_blck_lst <- hndl_mltp_opt(optlst$glbl_blck_lst)  # Global black list

			phntp_file <- optlst$phenotypes
			phntp_blckRows <- hndl_mltp_opt(optlst$phntp_blckRows)
			phntp_dscdCols <- hndl_mltp_opt(optlst$phntp_dscdCols)

			## Check the normality of raw data
			if(subcmd == "cknmlt") {
				# Phenotypes: ~/Documents/projects/200HIV/inputs/datasets/200HIV_Plasmamarkers_20190528_YANG.csv
				# Preprocessing on phenotypes
				exp_dtfm <- pprcs(phntp_file, dscd=phntp_dscdCols, blck=glbl_blck_lst, trps=FALSE)

				phntp_raw_nmltCk_flnm <- file.path(opt_plots_dir, optlst$phntp_raw_nmltCk_flnm)
				chck_nmlt(exp_dtfm, hst_opt=phntp_raw_nmltCk_flnm)

				# Save the preprocessed dataframe for the next step.
				# NOTE: Please try to avoid using comma in each field, for the saperator is comman in .csv files
				tmp_file <- file.path(tmp_dir, "cknmlt_tmp_dtfm.csv")
				fwrite(exp_dtfm, tmp_file, showProgress=FALSE, verbose=FALSE) 
				cat("[INFO]  Please check", phntp_raw_nmltCk_flnm, " for distribution of each measurement.\n")
				cat("[INFO]  Please check", tmp_file, " for preprocessed input data.\n")
			} else if(subcmd == "trfm") {
				# This step could be ran several times, therefore first check if there
				# is an processed file by this step (file name trfm_tmp_dtfm.tsv), if
				# not then go to check if there's a temporary file clled
				# cknmlt_tmp_dtfm.csv.
				# TODO: Perhaps an argument option to decide overlap the last proprecess or keep all
				## Transform the raw data by specific method: log2, log10, inverse-rank
				opt_tmp_file <- file.path(tmp_dir, "trfm_exotl_tmp.csv")
				ipt_tmp_file <- file.path(tmp_dir, "cknmlt_tmp_dtfm.csv")
				if(file.exists(opt_tmp_file)) {
					exp_dtfm <- pprcs(opt_tmp_file, trps=FALSE)
					cat("[WARN]  Will overwrite existing", opt_tmp_file, "\n")
				} else if (file.exists(ipt_tmp_file)) {
					exp_dtfm <- pprcs(ipt_tmp_file, trps=FALSE)
				} else {
					stop("Failed to find temporary file(cknmlt_tmp_dtfm.csv) from last step [cknmlt] or file(trfm_exotl_tmp.csv) from current step [trfm]")
				}

				phntp_trfm_mthd <- optlst$phntp_trfm_mthd 
				phntp_trfm_cols <- optlst$phntp_trfm_cols 
				if(phntp_trfm_cols == "*") {
					phntp_trfm_cols <- colnames(exp_dtfm) 
					phntp_trfm_cols <- phntp_trfm_cols[which(!phntp_trfm_cols %in% c("id"))]
				}
				exp_dtfm <- trfm_nmlt(exp_dtfm, trfm_cols=phntp_trfm_cols, trfm_mthd=phntp_trfm_mthd)

				## Check the normality of the transformed raw data
				phntp_trfm_nmltCk_flnm <- file.path(opt_plots_dir, optlst$phntp_trfm_nmltCk_flnm)
				chck_nmlt(exp_dtfm, hst_opt=phntp_trfm_nmltCk_flnm)
				cat("[INFO]  Please check", phntp_trfm_nmltCk_flnm, " for distribution of each measurement after transforming.\n")

				## Check outliers of transformed data using PCA
				phntp_trfm_otlCk_flnm <- file.path(opt_plots_dir, optlst$phntp_trfm_otlCk_flnm)
				pltpca(exp_dtfm, phntp_trfm_otlCk_flnm, idcol="id")  # Check the outliers using PCA
				cat("[INFO]  Please check", phntp_trfm_otlCk_flnm, " for PCA plots.\n")

				## Exclude outliers using 3 or 4 times of sd
				exp_dtfm <- prcs_otls(exp_dtfm, n_sd=4)
				phntp_exotl_nmltCk_flnm <- file.path(opt_plots_dir, optlst$phntp_exotl_nmltCk_flnm)
				chck_nmlt(exp_dtfm, hst_opt=phntp_exotl_nmltCk_flnm)
				cat("[INFO]  Please check", phntp_exotl_nmltCk_flnm, " for distribution of each measurement after excluding outliers.\n")

				fwrite(exp_dtfm, opt_tmp_file, showProgress=FALSE, verbose=FALSE) 
				cat("[INFO]  Please check", opt_tmp_file, " for preprocessed input data.\n")
			} else if(subcmd == "qtlmp") {
				phntp_ipt <- file.path(tmp_dir, "trfm_exotl_tmp.csv")
				if(file.exists(phntp_ipt)) {
					exp_dtfm <- pprcs(phntp_ipt, dscd=c("id"), trps=TRUE)
				} else {
					stop("Failed to find temporary file(trfm_exotl_tmp.csv) from last step [trfm]")
				}

				# Covariates
				# cvrt_file <- file.path("~", "Documents", "projects", "200HIV", "inputs", "datasets", "metaData_pcntgMnct_ssnlt_CD4CD8TC_2019Mar20.csv")
				cvrt_ipt <- optlst$covariates 
				cvrt_idxCol <- optlst$cvrt_idxCol 
				cvrt_dtfmTrp <- optlst$cvrt_dtfmTrp 
				cvrt_blckRows <- hndl_mltp_opt(optlst$cvrt_blckRows)
				cvrt_dscdCols <- hndl_mltp_opt(optlst$cvrt_dscdCols)

				cvrt_dtfm <- pprcs(cvrt_ipt, dscd=cvrt_dscdCols, blck=glbl_blck_lst)

				# QTL mapping
				# gntp_drct <- file.path("~", "Documents", "projects", "200HIV", "inputs", "dosage", "200HIV_dosages")
				gntp_ipt <- optlst$genotypes 
				gntp_ptn <- optlst$gntp_ptn 
				gntp_zip <- optlst$gntp_zip 
				gntp_blck_rows <- optlst$gntp_blck_rows 
				gntp_dscd_cols <- optlst$gntp_dscd_cols 

				qtl_dtfm <- data.frame()
				if(file_test("-d", gntp_ipt)) {
					if(is.null(gntp_ptn)) { stop("When --genotypes is a directory, gntp_ptn shouldn't be NULL or empty") }

					flnmlst <- list.files(gntp_ipt, pattern=gntp_ptn)
					if(length(flnmlst) == 0) { stop("--genotypes is a directory, but NO file following pattern (", gntp_ptn, ")") }

					for (flnm in flnmlst) {
						gntp_file <- file.path(gntp_ipt, flnm)
						gntp_dtfm <- pprcs(gntp_file, dscd=c("id", glbl_blck_lst), trps=FALSE, zipped=TRUE)
						qtl_dtfm <- rbind(qtl_dtfm, mkme(exp_dtfm, gntp_dtfm, cvrt_dtfm)$all$eqtls)
					}
				} else if(file_test("-f", gntp_ipt)) {
					gntp_dtfm <- pprcs(gntp_file, dscd=c("id", glbl_blck_lst), trps=FALSE, zipped=TRUE)
					qtl_dtfm <- mkme(exp_dtfm, gntp_dtfm, cvrt_dtfm)$all$eqtls
				} else {
					stop("--genotypes should be a directory, normal file, or zipped file")
				}

				qtl_tmp_file <- file.path(tmp_dir, "qtlmp_tmp_file.csv")
				fwrite(qtl_dtfm, qtl_tmp_file, showProgress=FALSE, verbose=FALSE) 
				cat("[INFO]  Please check", qtl_tmp_file, " for QTL list.\n")
			} else if(subcmd == "pltmht") {

				qtl_dtfm_flnm <- file.path(tmp_dir, "qtlmp_tmp_file.csv")
				if(file.exists(qtl_dtfm_flnm)) {
					qtl_dtfm <- pprcs(qtl_dtfm_flnm, trps=FALSE, as_idx=NULL)
				} else {
					stop("Failed to find temporary file", qtl_dtfm_flnm)
				}

				# Manhattan plots
				#	snp_ifmt <- file.path("~", "Documents", "projects", "200HIV", "inputs", "dosage", "200HIV_dosages", "variantInfo.gz")
				gntp_ifmt_flnm <- optlst$gntp_ifmt_flnm 
				gntp_ifmt_zip <- optlst$gntp_ifmt_zip 

				qtls <- intsct_qtlmt_ifmt(qtl_dtfm, gntp_ifmt_flnm)

				# Save some data to make life easy
				fwrite(qtls, file.path(tmp_dir, "qtls_gntpInfo_tmp_file.csv"),
							 row.names=FALSE, quote=FALSE, sep=",")

				sbjcts <- unique(qtls$gene)

				use_cols <- c("snps", "SequenceName", "Position", "pvalue")
				opt_pval_thrld <- 5e-2

				for(sbjct in sbjcts) {
					sbjct_qtls <- qtls[which(qtls$gene == sbjct), ]

					prfx <- file.path(opt_plots_dir, sbjct)
					pltme(sbjct_qtls, otpt_prfx=prfx, use_cols=use_cols)

					otpt_flnm <- file.path(opt_reports_dir, paste0(sbjct, ".tsv"))
					sbjct_qtls <- sbjct_qtls[which(sbjct_qtls[, 'pvalue'] < opt_pval_thrld), ]
					fwrite(sbjct_qtls, otpt_flnm, row.names=FALSE, quote=FALSE, sep="\t")
				}
			} else if(subcmd == "pltdsg") {
				# Dosage plots
				exp_dtfm_flnm <- file.path(tmp_dir, "trfm_exotl_tmp.csv")
				if(file.exists(exp_dtfm_flnm)) {
					exp_dtfm <- pprcs(exp_dtfm_flnm, trps=TRUE, as_idx="id")
				} else {
					stop("Failed to find temporary file", exp_dtfm_flnm)
				}

				qtls_dtfm_flnm <- file.path(tmp_dir, "qtls_gntpInfo_tmp_file.csv")
				if(file.exists(qtls_dtfm_flnm)) {
					qtls_dtfm <- pprcs(qtls_dtfm_flnm, trps=FALSE, as_idx=NULL)
				} else {
					stop("Failed to find temporary file", exp_dtfm_flnm)
				}

				gntp_ipt <- optlst$genotypes 
				gntp_ptn <- optlst$gntp_ptn 
				qtl_col <- c("snps", "gene", "pvalue", "SequenceName", "Position", "EffectAllele", "AlternativeAllele")
				plt_dsg(
					qtls_dtfm, exp_dtfm, gntp_ipt, qtl_col=qtl_col, pvt=5e-8,
					gntp_blck_lst=c("id", glbl_blck_lst), opt_dir=opt_plots_dir
				)
			} else {
				cat("[]")
			}
		}, error=function(e) {
			cat("[ERROR]  ", e$message, "\n")
			print(e)
			cat("[WARN]  NO cleaning-up ...\n")
		}, finally=function(m) {
			cat("[INFO]  Cleaning up ...\n")
			unlink(tmp_dir, TRUE, TRUE)
		}
	)
}

main()

# /vim:ts=8:ft=R/
