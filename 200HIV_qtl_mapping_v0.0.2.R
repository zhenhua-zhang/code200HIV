#!/usr/bin/env Rscript
###############################################################################
## @version: 0.0.1
## @author: Zhenhua Zhang
## @eamil: zhenhua.zhang217@gmail.com
## @refer: qtl_mapping.R (Xiaojing Chu)
## @date: 4th June, 2019
###############################################################################


# TODO: 1. When doing inverse-rank transform, use mean to replace NA, which is
#		dangerous. Need a better way to deal with NAs
#       2. Make the script rescuable at broken-points with additional command
#       arguments

# NOTE: 1. Settings in config file will be overrided by command line arguments
#       2. Some settings is only modifiable by configuration files
#		3. Keep the row as measurements, column as subjects
# 		4. Please try to avoid using comma in each field, for the saperator is
# 		comman in .csv files


# Used packages
# TODO: 1. Add version. 2. Suppress warnings(??)
library(qqman)
library(dplyr)
library(ggplot2)
library(methods)
library(stringi)
library(optparse)
library(reshape2)
library(data.table)
library(MatrixEQTL)


# Processing the configuration file
prcs_cfg <- function(iptfl) {
	cfglst <- list()
	cfgs <- paste(readLines(iptfl), sep="\n")
	for(line in cfgs) {
		if(!startsWith(line, "#") && !stri_isempty(line)) {
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


################################################################################  
##    Preprocess                                                              ##
################################################################################  
pprcs <- function(flnm, dscd_cols=NULL, dscd_rows=NULL, kept_cols=NULL,
				  kept_rows=NULL, trps=TRUE, as_idx="id", zipped=FALSE) {
	cat("[INFO]  Processing: ", flnm, "\n")

	dtfm <- smt_fread(flnm, zipped=zipped)

	if(!is.null(as_idx)){
		base::row.names(dtfm) <- dtfm[, as_idx] 
	}

	col_names <- base::colnames(dtfm)
	row_names <- base::rownames(dtfm)

	dtfm <- dtfm[, base::sort(col_names)]
	col_names <- base::colnames(dtfm)

	if (is.null(dscd)) {
		cat("[INFO]  Columns discarding list is empty ...\n")
	} else {
		if(is.null(kept_cols)) {
			cat("[INFO]  Columns keeping list is empty ...")
		} else {
			ovlps <- dscd_cols[which(dscd_cols %in% kept_cols)]
			if(is.null(ovlp) || length(ovlp)!=0) {
				stop("Overlap between dscd_cols and kept_cols: ", ovlps)
			}
		}
		cat("[WARN]  Will discard column(s):", dscd, "...\n")
		kept_cols <- c(kept_cols, col_names[which(!col_names %in% dscd_rows)])
		dtfm <- dtfm[, kept_cols]
	}

	if (is.null(blck)) {
		cat("[INFO]  Row discarding list is empty ...\n")
	} else {
		if(is.null(kept_rows)) {
			cat("[INFO]  Rows keeping list is empty ...")
		} else {
			ovlps <- dscd_cols[which(dscd_rows %in% kept_rows)]
			if(is.null(ovlp) || length(ovlp)!=0) {
				stop("Overlap between dscd_rows and kept_rows: ", ovlps)
			}
		}
		cat("[WARN]  Will discard rows(s):", blck, "...\n")
		kept_rows <- c(kept_rows, col_names[which(!row_names %in% dscd_rows)])
		dtfm <- dtfm[kept_rows, ]
	}


	if(trps) {
		cat("[WARN]  Will transpose dataframe taking", as_idx, "as character\n")
		idx_val <- dtfm[, as_idx]
		dtfm <- dtfm[, !(base::colnames(dtfm) %in% c(as_idx))]
		dtfm <- as.data.frame(t(dtfm), stringsAsFactors=FALSE)
		dtfm[as_idx, ] <- idx_val
		return(dtfm)
	} else {
		return(dtfm)
	}
}


###############################################################################
##    Check the normality of each measuements                                ##
###############################################################################
chck_nmlt <- function(dtfm, hst_opt="Normality_check.pdf", bins=30, idv='id',
					  stclnm=TRUE, fgsz=3) {
	cat("[INFO]  Check the normality...\n")

	if(!is.data.frame(dtfm)) {
		stop(" `dtfm` should be a dataframe ...")
	}

	if(stclnm) {
		dtfm <- dtfm[, base::sort(base::colnames(dtfm))]
	}

	exp_mlt <- reshape2::melt(dtfm, id.vars=idv)

	nvars <- length(unique(exp_mlt["variable"])[[1]])
	ncol <- as.integer(sqrt(nvars))

	p <- ggplot(data=exp_mlt)
	p <- p + geom_histogram(aes(x=value), na.rm=TRUE, bins=bins)
	p <- p + facet_wrap(.~variable, ncol=ncol, scales="free")

	ggsave(hst_opt, width=ncol*fgsz, height=ncol*fgsz)
}


###############################################################################
##    Transform                                                              ##
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

# Trans form log2 values into log10 values
lg2tolg10 <- function(x) { return(x * log10(2)) }


# NOTE: Make sure only the id column is characters, while the rest are numeric,
#       otherwise, the transpose would fail
# Transform non-normal distributed phnotype column into
trfm_nmlt <- function(dtfm, dscd_cols=NULL, trfm_cols=NULL, trfm_mthd="log2",
					  trps=FALSE, idcol="id", stclnm=TRUE) {
	cat("[INFO]  Transform data by", trfm_mthd, "\n")

	id_val <- dtfm[idcol]
	dtfm[idcol] <- NULL

	row_names <- base::rownames(dtfm)
	if(stclnm) {
		dtfm <- dtfm[base::sort(base::colnames(dtfm))]
	}
	col_names <- colnames(dtfm)

	if(!is.null(dscd_cols)) {
		cat("[WARN]  Will discard: \n         ", dscd_cols, "\n")
		dtfm <- dtfm[, which(!col_names %in% dscd_cols)]
		col_names <- base::colnames(dtfm)
	}

	if(!is.null(trfm_cols) && !is.na(trfm_cols)) {
		if(!is.null(dscd_cols)){
			ovlp <- trfm_cols[which(dscd_cols %in% trfm_cols)]
			if(length(ovlp) != 0) {
				stop("Overlapped between dscd_cols and trfm_cols \n", ovlp)
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
			function(x) {
				if(x %in% trfm_cols) {return(paste(x, trfm_mthd, sep="_"))}
				else {return(x)}
			}
		)
	} else {
		cat("[WARN]  trfm_cols is either NULL or NA\n")
	}

	if(trps) {
		cat("[WARN]  Will transpose dataframe taking", idcol, "as character\n")
		dtfm <- as.data.frame(t(dtfm), stringsAsFactors=FALSE)
	}

	dtfm <- as.data.frame(scale(dtfm), stringsAsFactors=FALSE)
	dtfm[idcol] <- id_val
	return(dtfm)
}


# Draw PCA
plt_pca <- function(dtfm, opt_flnm="Outlier_check_PCA.pdf", idcol="id",
				   fgwd=10, fght=10, fgunt="in") {
	cat("[INFO]  Plotting PCA... \n")

	if(!is.data.frame(dtfm)) {
		stop("[ERROR]  `dtfm` should be a dataframe ...")
	}

	dtfm[idcol] <- NULL
	dtfm[is.na(dtfm)] <- mean(base::rowMeans(dtfm), na.rm=TRUE)
	pca <- prcomp(dtfm)
	pcax <- as.data.frame(pca$x, stringsAsFactors=FALSE)
	pcac <- pca$center

	p <- ggplot() + theme_bw()
	p <- p + geom_point(data=pcax, mapping=aes(x=PC1, y=PC2))

	ggsave(opt_flnm, width=fgwd, height=fght, units=fgunt)
}


# Remove outliers (out of 3*sd or 4*sd)
prcs_otls <- function(dtfm, n_sd=3, trps=FALSE, idcol="id") {
	cat("[INFO]  Removing outliers...\n")

	if(!is.null(idcol)) {
		id_val <- dtfm[idcol]
		dtfm[idcol] <- NULL
	}

	row_names <- base::rownames(dtfm)
	col_names <- base::colnames(dtfm)

	dtfm <- sapply(dtfm,
		function(x) {
			sd <- sd(x, na.rm=TRUE)
			mu <- mean(x, na.rm=TRUE)
			x[which((x <= mu - n_sd * sd) | (x >= mu + n_sd * sd))] <- NA
			return(x)
		}
	)

	dtfm <- as.data.frame(dtfm, row.names=row_names, stringsAsFactors=FALSE)
	cat("[INFO]  Number of masked subjects in each measurments:\n")
	nacnt <- sapply(dtfm, function(x) {sum(is.na(x))})

	col_names <- base::colnames(dtfm)
	if(trps) {
		cat("[WARN]  Will transpose dataframe taking", idcol, "as character\n")
		dtfm <- as.data.frame(t(dtfm), stringsAsFactors=FALSE)
	}

	dtfm[idcol] <- id_val
	dtfm <- dtfm[c(idcol, col_names)]
	return(dtfm)
}


# Create SliceDate object
mkslcdt <- function(ipt, dlmt=",", omchr="NA", skprw=0, skpcl=0,
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
mkme <- function(pntp, gntp, cvrt=NULL, opt_file=NULL, pv_opt_thr=5e-2) {
	cat("[INFO]  Mapping QTL ...\n")
	useModel = modelLINEAR
	errorCovariance <- numeric()
	print(cvrt[1:5, 1:5])

	gene <- mkslcdt(pntp)
	snps <- mkslcdt(gntp, dlmt="\t")
	cvrt <- mkslcdt(cvrt)

	me <- Matrix_eQTL_engine(
		 snps=snps, gene=gene, cvrt=cvrt,
		 output_file_name=opt_file, pvOutputThreshold=pv_opt_thr,
		 useModel=useModel, errorCovariance=errorCovariance, verbose=FALSE,
		 pvalue.hist=TRUE, min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE
	)

	return(me)
}


# Intersect Matrix eQTL with SNP information
# TODO: A better function name
itsct_qtlmt_ifmt <- function(qtls, snp_ifmt, mgby_x="snps", mgby_y="rsID") {
	cat("[INFO]  Intersect QTL matrix with SNP information ...\n")
	if(is.character(qtls)) { qtls <- smt_fread(qtls) }
	if(is.character(snp_ifmt)) { snp_ifmt <- smt_fread(snp_ifmt) }
	return(merge(qtls, snp_ifmt, by.x=mgby_x, by.y=mgby_y))
}


# Draw Manhattan plot and Q-Q plot
plt_mht <- function(qtls_lst, use_cols=NULL, opt_prfx="qtl_mapping",
					 svfmt="pdf", wdth=20, hght=10, drwzsc=FALSE, drwqq=FALSE) {
	cat("[INFO]  Plotting Manhattan plot and Q-Q plot for", opt_prfx, "...\n")
	# use_cols should follow the order: snp_id, chromosome, position, p-value
	if(is.null(use_cols)) {
		use_cols <- c("snps", "chr", "pos", "pvalue") 
	} else if(length(use_cols) != 4) {
		stop("[ERROR]  The length of use_cols should be 4 ...\n")
	}

	gwas_rslt <- qtls_lst[, use_cols]
	base::colnames(gwas_rslt) <- c("SNP", "CHR", "BP", "P")
	opt_flnm <- paste0(opt_prfx, "_mhtn_pvalue.", svfmt)

	if(svfmt=="png") { png(opt_flnm, width=wdth, height=hght) }
	else { pdf(opt_flnm, width=wdth, height=hght) }

	manhattan(gwas_rslt, main="Manhattan Plot (p-value)", ylan="p-value(-log10)",
						suggestiveline=FALSE, annotatePval=5e-8, annotateTop=FALSE
						)
	dev.off()

	if(drwzsc) {
		gwas_rslt <- transform(gwas_rslt, zscore=qnorm(P/2, lower.tail=FALSE))
		opt_flnm <- paste0(opt_prfx, "_mhtn_zscore.", svfmt)
		if(svfmt=="png") { png(opt_flnm, width=wdth, height=hght) }
		else { pdf(opt_flnm, width=wdth, height=hght) }
		manhattan(gwas_rslt, main="Manhattan Plot (z-score)", ylab="Z-score",
							suggestiveline=FALSE, annotatePval=0.00000005
							)
		dev.off()
	} else {
		cat("[INFO]  Skipping Manhattan plot for z-score ...\n")
	}

	if(drwqq) {
		opt_flnm <- paste0(opt_prfx, "_qq.", svfmt)
		if(svfmt=="png") { png(opt_flnm, width=wdth, height=hght) }
		else { pdf(opt_flnm, width=wdth, height=hght) }
		qq(gwas_rslt$P, main="Q-Q plot of p-values", pch=18, cex=1.5, las=1)
		dev.off()
	} else {
		cat("[INFO]  Skipping Q-Q plot ...\n")
	}
}


# Filter out check the number of mode
ckmd <- function(dtfm, thrshld=0.5) {
	if(is.data.frame(dtfm)) {
		res <- sapply(dtfm, 
			function(x) {
				cnt <- table(x)
				blvct <- ifelse(max(cnt)/sum(cnt) >= thrshld, TRUE, FALSE)
				return(blvct)
			}
		)
		return(res)
	} else if (is.vector(dtfm)) {
		cnt <- table(dtfm)
		return(ifelse(max(cnt)/sum(cnt) >= thrshld, TRUE, FALSE))
	} else {
		stop("[ERROR]  I can only handle vector or data.frame!!!")
	}
}

# Encode genotype
enc_gntp <- function(x, use_gntp, use_pntp) {
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


# Draw dosage boxplot
plt_dsg <- function(qtl_dtfm, pntp, gntp, gntp_ptn="*.dosage", gntp_zipped=TRUE,
					svfmt="pdf", pvt=5e-8, qtl_col=NULL, gntp_blck_lst=NULL,
					opt_dir="qtl_mapping_output_dir/plots") {

	# Fetch top snp on each chr, filted by pvalue
	if(is.null(qtl_col)) { 
		qtl_col <- c( "snps", "gene", "pvalue", "chr", "pos", "ref", "alt") 
	}

	qtl_dtfm <- qtl_dtfm[, qtl_col]
	qtl_col <- c("snps", "gene", "pvalue", "chr", "pos", "ref", "alt") 
	base::colnames(qtl_dtfm) <- qtl_col

	qtl_gb <- c("snps", "chr", "gene")

	if(is.vector(qtl_col) && length(qtl_col)==7){
		top_snps_per_chr <- qtl_dtfm %>%
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
	# A dataframe with column as subjects, rows as measurements
	if(is.data.frame(pntp)) {  
		pntp_rwnm <- base::rownames(pntp)
		use_pntp <- pntp[which(pntp_rwnm %in% use_msmnt), ]
	} else {
		stop("pntp should be a dataframe")
	}

	# Fetch genotype
	# Repeated code, should rewrite into a function
	if(file_test("-d", gntp)) {
		if(is.null(gntp_ptn)) {
			stop("When gntp is a dir, gntp_ptn is required")
		}

		flnmlst <- list.files(gntp, pattern=gntp_ptn)
		if(length(flnmlst) == 0) {
			stop("gntp is a dir, but NO file following pattern:", gntp_ptn)
		}

		gntp_dtfm <- data.frame()
		for (flnm in flnmlst) {
			gntp_file <- file.path(gntp, flnm)
			tmp_dtfm <- pprcs(gntp_file, dscd_cols=gntp_blck_lst, trps=FALSE)
			gntp_dtfm <- rbind(gntp_dtfm, tmp_dtfm)
		}
	} else if(file_test("-f", gntp)) {
		gntp_dtfm <- pprcs(gntp_file, dscd_cols=gntp_blck_lst, trps=FALSE)
	} else {
		stop("gntp should be a directory, normal file, or zipped file")
	}
	gntp_rwnm <- base::rownames(gntp_dtfm)
	use_gntp <- gntp_dtfm[which(gntp_rwnm %in% use_snps), ]

	# print(rbind(use_gntp, use_pntp))
	lgt <- dim(top_snps_per_chr)[1]
	if(lgt != 0) {
		gnpnlst <- apply(top_snps_per_chr, 1, enc_gntp, use_gntp, use_pntp)
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
hdmtopt <- function(arg) {
	if (!is.character(arg)) {
		stop("Keyword argument arg should be character!!") 
	}
	if (length(arg) == 0) {
		cat("[WARN]  Keyword argument arg is empty, return \"\"\n")
		return("")
	}
	if (is.na(arg)) { 
		cat("[WARN]  Keyword argument arg is NA, return NULL\n")
		return(NULL) 
	}

	blnks <- c(" ", "\t", "\n")  # Check blanks
	if (any(stri_detect_fixed(arg, blnks))) {
		stop("Blank or tab is not allowed!!")
	} 

	arg <- stri_split_fixed(arg, ";")[[1]]
	arg <- stri_split_fixed(arg, ",")

	return(arg)
}


# Parsing command line arguments
prsarg <- function() {
	tryCatch(
		{
			parser <- OptionParser()

			parser <- add_option(
				parser, c("-G", "--genotypes"), type="character",
				help="Genotypes file or directory"
			)
			parser <- add_option(
				parser, c("-C", "--covariates"), type="character",
				help="Covariates"
			)
			parser <- add_option(
				parser, c("-P", "--phenotypes"), type="character",
				help="Phenotypes file"
			)
			parser <- add_option(
				parser, c("-S", "--config-file"), type="character",
				default="config.cfg",
				help="File including settings, overrided by command arguments"
			)
			parser <- add_option(
				parser, c("-T", "--tmp-dir"), type="character",
				help="File for temporary files"
			)
			parser <- add_option(
				parser, c("-o", "--output-dir"), type="character",
				help="Output directory for the QTL mapping results"
			)
			agmnts <- parse_args(
				parser, positional_arguments=1,
				convert_hyphens_to_underscores=TRUE
			)

			# cknmlt: check normality; trfm: Transfomr data;
			# qltmp: QTL mapping; qtlrpt: produce QTL report
			subcmd <- agmnts$args
			subcmdlst <- c("cknmlt", "trfm", "qtlmp", "qtlrpt")
			if(!subcmd %in% subcmdlst) {
				stop("Unknown: ", subcmd, "Opts:", subcmdlst)
			}
		}, error=function(e) {
			cat("[ERROR] ", e$message, "\n")
			print(e)
		}, finally=function(msg) {
			cat("[INFO]  Finished tryCatch scope", msg, "\n")
			return(agmnts)
		}
	)
}


# Setup working space
stwkspc <- function(tmp_dir=NULL, opt_dir=NULL, mode="0755") {
	# OPT_DIR
	# |-- logs
	# |-- reports
	# |   |-- GenotypeLevelPlot
	# |   |-- QtlInformation
	# |   |-- LocusZoomPlot
	# |   |-- ManhattanPlot
	# |   |-- Preprocessing
	# |   |-- Pathway
	# |-- README.md

	if(dir.exists(tmp_dir)) {
		cat("[WARN] ", tmp_dir, "exists\n") 
	} else {
		dir.create(tmp_dir, mode=mode)
	}

	if(dir.exists(opt_dir)) { 
		cat("[WARN] ", opt_dir, "exists\n") 
	} else {
		dir.create(opt_dir, mode=mode)
		dir.create(file.path(opt_dir, "logs"), mode=mode)
		dir.create(file.path(opt_dir, "reports"), mode=mode)
		dirs <- c(
			"GenotypeLevelPlot", "QtlInformation", "LocusZoomPlot",
			"ManhattanPlot", "Preprocessing", "Pathway"
		)
		for(dir in dirs) {
			dir.create(file.path(opt_dir, dir), mode=mode)
		}
	}
}


# Merge the settings from configuration file and command line arguments
mgcfgcmd <- function(agmntlst, cfglst) {
	if(!is.list(agmntlst) || !is.list(cfglst)) {
		stop("Keyword arguments both agmntlst and cfglst should be list")
	}
	cmdopts <- agmntlst$options
	cmdopts["subcmd"] <- agmntlst$args
	cfglst[names(cmdopts)] <- cmdopts
	return(cfglst)
}


# sub-command `chnmlt`
cknmlt <- function(optlst) {
	# Arguments
	pntp_file <- optlst$phenotypes
	pntp_dtfm_trp <- optlst$pntp_dtfm_trp
	pntp_dscd_cols <- hdmtopt(optlst$pntp_dscd_cols)
	pntp_dscd_rows <- hdmtopt(optlst$pntp_dscd_rows)
	glbl_blck_lst <- c(pntp_dscd_rows, optlst$glbl_blck_lst)

	# Preprocessing
	pntp_dtfm <- pprcs(
		pntp_file, dscd_cols=pntp_dscd_cols, dscd_rows=pntp_dscd_rows,
		trps=FALSE
	)

	# Check normality
	pntp_raw_nmltCk_flnm <- file.path(
		optlst$opt_pprc_dir, optlst$pntp_raw_nmltCk_flnm
	)
	chck_nmlt(pntp_dtfm, hst_opt=pntp_raw_nmltCk_flnm)
	cat("[INFO]  Check", pntp_raw_nmltCk_flnm, "for distributions.\n")

	# Save the preprocessed dataframe for the next step.
	cknmlt_tmp_dtfm <- file.path(tmp_dir, "cknmlt_tmp.csv")
	fwrite(pntp_dtfm, tmp_file, showProgress=FALSE, verbose=FALSE) 
	cat("[INFO]  Check", cknmlt_tmp_dtfm, "for preprocessed input data.\n")
}


# sub-command `trfm`
trfm <- function() {
	# TODO: Perhaps an argument option to decide overlap the last proprecess or
	# keep all

	tmp_dir <- optlst$tmp_dir
	opt_dir <- optlst$opt_dir

	# Read temp-file from the last step
	ipt_tmp_flnm <- file.path(tmp_dir, "cknmlt_tmp.csv")
	if (file.exists(ipt_tmp_flnm)) {
		pntp_dtfm <- pprcs(ipt_tmp_flnm, trps=FALSE)
	} else {
		stop("Failed to find temporary file", ipt_tmp_flnm)
	}

	# Transform the raw data by specific method: log2, log10, inverse-rank
	pntp_trfm_mthd <- hdmtopt(optlst$pntp_trfm_mthd)
	pntp_trfm_cols <-  hdmtopt(optlst$pntp_trfm_cols)
	pntp_trfm_excols <- hdmtopt(optlst$pntp_trfm_excols)

	if(pntp_trfm_cols == "*") {
		pntp_trfm_cols <- colnames(pntp_dtfm) 
		excld <- which(!pntp_trfm_cols %in% pntp_trfm_excols)
		pntp_trfm_cols <- pntp_trfm_cols[excld]
	}

	pntp_dtfm <- trfm_nmlt(
		 pntp_dtfm, trfm_cols=pntp_trfm_cols, trfm_mthd=pntp_trfm_mthd
	)

	## Check the normality of the transformed raw data
	pntp_trfm_nmltCk_flnm <- file.path(
		optlst$opt_pprc_dir, optlst$pntp_trfm_nmltCk_flnm
	)
	chck_nmlt(pntp_dtfm, hst_opt=pntp_trfm_nmltCk_flnm)
	cat(
		"[INFO]  Check", pntp_trfm_nmltCk_flnm,
		"for distributions after transforming.\n"
	)

	## Check outliers of transformed data using PCA
	pntp_trfm_otlCk_flnm <- file.path(
		optlst$opt_pprc_dir, optlst$pntp_trfm_otlCk_flnm
	)
	plt_pca(pntp_dtfm, pntp_trfm_otlCk_flnm, idcol="id")  
	cat("[INFO]  Check", pntp_trfm_otlCk_flnm, "for PCA plots.\n")

	## Exclude outliers using 3 or 4 times of sd
	pntp_otl_n_sd <- as.numeric(optlst$pntp_otl_n_sd)
	pntp_dtfm <- prcs_otls(exp_dtfm, n_sd=pntp_otl_n_sd)

	## Check normality
	pntp_exotl_nmltCk_flnm <- file.path(
		optlst$opt_pprc_dir, optlst$pntp_exotl_nmltCk_flnm
	)
	chck_nmlt(pntp_dtfm, hst_opt=pntp_exotl_nmltCk_flnm)
	cat(
		"[INFO]  Check", pntp_exotl_nmltCk_flnm,
		"for distributions after excluding outliers.\n"
	)

	## Save processed data frame into temporary file
	opt_tmp_file <- file.path(optlst$tmp_dir, "trfm_tmp.csv")
	fwrite(pntp_dtfm, opt_tmp_file, showProgress=FALSE, verbose=FALSE) 
	cat("[INFO]  Check", opt_tmp_file, " for preprocessed input data.\n")
}


# sub-command `qtlmp`
qtlmp <- function(optlst) {
	tmp_dir <- optlst$tmp_dir
	opt_dir <- optlst$opt_dir
	glbl_blck_lst <- optlst$glbl_blck_lst

	ipt_tmp_flnm <- file.path(tmp_dir, "trfm_tmp.csv")
	if(file.exists(ipt_tmp_flnm)) {
		pntp_dtfm <- pprcs(ipt_tmp_flnm, dscd_cols=c("id"), trps=TRUE)
	} else {
		stop("Failed to fined temporary file", ipt_tmp_flnm)
	}

	# Covariates
	cvrt_ipt <- optlst$covariates
	cvrt_trp <- optlst$cvrt_trp
	cvrt_zip <- optlst$cvrt_zip
	cvrt_idx_col <- optlst$cvrt_idx_col
	cvrt_dtfm_trp <- optlst$cvrt_dtfm_trp
	cvrt_dscd_cols <- hdmtopt(optlst$cvrt_dscd_cols)
	cvrt_dscd_rows <- c(glbl_blck_lst, hdmtopt(optlst$cvrt_dscd_rows))

	cvrt_dtfm <- pprcs(
		cvrt_ipt, dscd_cols=cvrt_dscd_cols, dscd_rows=cvrt_dscd_rows
	)

	# Genotypes
	gntp_ipt <- optlst$genotypes
	gntp_ptn <- optlst$gntp_ptn
	gntp_zip <- optlst$gntp_zip  # Not piped into function, yet
	gntp_idx_col <- optls$gntp_idx_col
	gntp_dtfm_trp <- optlst$gntp_dtfm_trp  # Not piped into function, yet
	gntp_dscd_rows <- hdmtopt(optlst$gntp_dscd_rows)
	gntp_dscd_cols <- c(hdmtopt(optlst$gntp_dscd_cols), glbl_blck_lst)

	qtl_dtfm <- data.frame()
	if(file_test("-d", gntp_ipt)) {
		if(is.null(gntp_ptn)) {
			stop("--genotypes is a directory, gntp_ptn is required")
		}

		flnmlst <- list.files(gntp_ipt, pattern=gntp_ptn)
		if(length(flnmlst) == 0) {
			stop("--genotypes is a directory, but NO file by: ", gntp_ptn)
		}

		for (flnm in flnmlst) {
			gntp_file <- file.path(gntp_ipt, flnm)
			gntp_dtfm <- pprcs(
				gntp_file, dscd_cols=gntp_dscd_cols, trps=FALSE, zipped=TRUE
			)
			qtl_dtfm <- rbind(
				qtl_dtfm, mkme(pntp_dtfm, gntp_dtfm, cvrt_dtfm)$all$eqtls
			)
		}
	} else if(file_test("-f", gntp_ipt)) {
		gntp_file <- gntp_ipt
		gntp_dtfm <- pprcs(
			gntp_file, dscd_cols=gntp_dscd_cols, trps=FALSE, zipped=TRUE
		)
		qtl_dtfm <- mkme(pntp_dtfm, gntp_dtfm, cvrt_dtfm)$all$eqtls
	} else {
		stop("--genotypes should be a directory, normal file, or zipped file")
	}

	# Save QTL mapping result into temporary file
	opt_tmp_file <- file.path(optlst$tmp_dir, "qtlmp_tmp.csv")
	fwrite(qtl_dtfm, opt_tmp_file, showProgress=FALSE, verbose=FALSE) 
	cat("[INFO]  Please check", opt_tmp_file, " for temporary QTL list.\n")
}


# sub-command `qtlrpt`
qtlrpt <- function(optlst) {
	tmp_dir <- optlst$tmp_dir
	opt_dir <- optlst$opt_dir
	glbl_blck_lst <- optlst$glbl_blck_lst

	# Read temp-file from last step
	ipt_qtl_flnm <- file.path(tmp_dir, "qtlmp_tmp.csv")
	if(file.exists(ipt_qtl_flnm)) {
		qtl_dtfm <- pprcs(ipt_qtl_flnm, trps=FALSE, as_idx=NULL)
	} else {
		stop("Failed to find temporary file", ipt_qtl_flnm)
	}

	# QTLs with genotype information
	gntp_ifmt_flnm <- optlst$gntp_ifmt_flnm
	gntp_ifmt_zip <- optlst$gntp_ifmt_zip

	# qtl_dtfm is override
	qtl_dtfm <- itsct_qtlmt_ifmt(qtl_dtfm, gntp_ifmt_flnm)

	## Save QTLs for each measurements and draw Manhattan plots
	qtlrpt_mhtnplt_cols <- hdmtopt(optlst$qtlrpt_mhtnplt_cols)
	qtlrpt_mhtnplt_svfmt <- optlst$qtlrpt_mhtnplt_svfmt
	qtlrpt_opt_pvt <- as.numeric(optlst$qtlrpt_opt_pvt)
	pvalue <- qtlrpt_mhtnplt_cols[4]

	msmnts <- unique(qtl_dtfm$gene)
	for(msmnt in msmnts) {
		msmnt_dtfm <- qtl_dtfm[which(qtl_dtfm$gene == msmnt), ]

		opt_prfx <- file.path(opt_mhtn_dir, msmnt)
		plt_mht(msmnt_dtfm, opt_prfx=opt_prfx, use_cols=qtlrpt_mhtnplt_cols)
		cat("[INFO] Check", paste0(opt_prfx, opt_prfx), "for Manhattan plot\n")

		opt_flnm <- file.path(optlst$opt_qtl_dir, paste0(msmnt, ".csv"))
		msmnt_dtfm <- msmnt_dtfm[which(msmnt_dtfm[pvalue] < qtlrpt_opt_pvt), ]
		fwrite(msmnt_dtfm, opt_flnm, row.names=FALSE, quote=FALSE, sep=",")
		cat("[INFO] Check", opt_flnm, "for QTL informations\n")
	}

	#  Plot genotype level for top SNPs
	## Read temp-files
	pntp_idx_col <- optlst$pntp_idx_col
	ipt_pntp_flnm <- file.path(tmp_dir, "trfm_tmp.csv")
	if(file.exists(ipt_pntp_flnm)) {
		pntp_dtfm <- pprcs(ipt_pntp_flnm, trps=TRUE, as_idx=pntp_idx_col)
	} else {
		stop("Failed to find temporary file", ipt_pntp_flnm)
	}

	gntp_ipt <- optlst$genotypes 
	gntp_ptn <- optlst$gntp_ptn 
	gntp_zip <- optlst$gntp_zip
	gntp_dscd_cols <- c(hdmtopt(optlst$gntp_dscd_cols), glbl_blck_lst)
	qtlrpt_dsgplt_pvt <- as.numeric(optlst$qtlrpt_dsgplt_pvt)
	qtlrpt_dsgplt_cols <- hdmtopt(optlst$qtlrpt_dsgplt_cols)

	plt_dsg(
		qtl_dtfm, pntp_dtfm, gntp_ipt, qtl_col=qtlrpt_dsgplt_cols,
		gntp_ptn=gntp_ptn, pvt=gntp_dsg_pvt, gntp_blck_lst=gntp_dscd_cols,
		opt_dir=opt_dsg_dir
	)
}


# Create a README.md file to descript current
mkrdme <- function() {
}


# Main function
main <- function() {
	# Trycatch clean up tmp_dir
	tryCatch(
		{
			# Command line arguments
			agmntlst <- prsarg()
			opts <- agmntlst$options
			cfglst <- prcs_cfg(opts$config_file)
			optlst <- mgcfgcmd(agmntlst, cfglst)

			# Err files
			errobj <- file(file.path(opt_lg_dir, paste0(subcmd, ".err")))
			sink(errobj, append=TRUE, type="message")

			# Log files
			logobj <- file(file.path(opt_lg_dir, paste0(subcmd, ".log")))
			sink(logobj, append=TRUE, type="output")

			tmp_dir <- optlst$tmp_dir
			opt_dir <- optlst$opt_dir

			stwkspc(tmp_dir, opt_dir)

			optlst$opt_lg_dir <- file.path(opt_dir, "logs")
			optlst$opt_rpt_dir <- file.path(opt_dir, "reports")
			optlst$opt_dsg_dir <- file.path(opt_rpt_dir, "GenotypeLevelPlot")
			optlst$opt_qtl_dir <- file.path(opt_rpt_dir, "QtlInformation")
			optlst$opt_mhtn_dir <- file.path(opt_rpt_dir, "ManhattanPlot")
			optlst$opt_pprc_dir <- file.path(opt_rpt_dir, "Preprocessing")

			# Sub-command
			subcmd <- optlst$subcmd

			# Global black list
			optlst$glbl_blck_lst <- hdmtopt(optlst$glbl_blck_lst)

			## Check the normality of raw data
			if(subcmd == "cknmlt") {
				cknmlt(optlst)
			} else if(subcmd == "trfm") {
				trfm(optlst)
			} else if(subcmd == "qtlmp") {
				qtlmp(optlst)
			} else if(subcmd == "qtlrpt") {
				qtlrpt(optlst)
			} else {
				stop("Unknow error while processing sub-command: ", subcmd)
			}
		}, error=function(e) {
			cat("[ERROR]  ", e$message, "\n")
			print(e)
		}, finally=function(m) {
			cat("[INFO]  Success!!!\n\n------------------------------------")
		}
	)
}

main()

