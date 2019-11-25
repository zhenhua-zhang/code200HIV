#!/usr/bin/env Rscript

# NOTE: 
#    1. The id column is mandatory for each input file, i.e phenotypes, covariates
#    2. For the phenotype and the covariate file, the script supposes the columns are traits, while
#    each row represents one sample. Therefore, if the --phenotypeDataFrameTranspose or
#    --covariateDataFrameTranspose is given, it means the input file isn't in the formated as the
#    script supposes.

# TODO: 
#    1. The format of input files
#    2. Change all variable name etc. into snake_case instead of camlCase.
# 3. A README.md to descript this shit.
# 4. Create a global configuration file to control options of lintr.

# FIXME:
# 1. Issue with lintr. The gitter marks are in wrong place.

suppressPackageStartupMessages(library(qqman, quietly = TRUE))
library(stringr)
library(ggplot2)
library(optparse)
library(data.table)
library(MatrixEQTL)

#
## Plot phenotype level per genotype
#
# TODO:
# 	1. Need a smart implementation
#	2. Remove this funciton and pipe it into a seperated file.
pppg <- function(snpid, phenotype, qtlsInformation, phenotypeLevels, commonSamples, genotypeDosageRound) {
    myConditions <- qtlsInformation$SNP == snpid & qtlsInformation$gene == phenotype
    currentQtlInfo <- qtlsInformation[myConditions, ]

    print(currentQtlInfo)
    if (base::nrow(currentQtlInfo) != 0) {
        alternativeAllele <- currentQtlInfo[["AlternativeAllele"]]
        effectAllele <- currentQtlInfo[["EffectAllele"]]
        pValue <- currentQtlInfo[["p-value"]]
        beta <- currentQtlInfo[["beta"]]
        chr <- currentQtlInfo[["SequenceName"]]
        pos <- currentQtlInfo[["Position"]]

        topsnpDosageRound <- as.data.frame(genotypeDosageRound[snpid, commonSamples])
        topsnpDosageRoundPhenotype <- as.data.frame(t(rbind.data.frame(topsnpDosageRound, phenotypeLevels)))

        topsnpDosageRoundPhenotypes[, snpid] <- sapply(topsnpDosageRoundPhenotype[, snpid], function(e){
            if (e == 0){
                return(paste0(alternativeAllele, alternativeAllele))
            } else if(e == 1) {
                return(paste0(effectAllele, alternativeAllele))
            } else if(e == 2) {
                return(paste0(effectAllele, effectAllele))
            } else {
                warning("Bad genotype coding occured!")
                return(NA)
            }
        })
        # topsnpDosageRoundPhenotype[topsnpDosageRoundPhenotype[, snpid] == 0, snpid] <- paste0(alternativeAllele, alternativeAllele)
        # topsnpDosageRoundPhenotype[topsnpDosageRoundPhenotype[, snpid] == 1, snpid] <- paste0(effectAllele, alternativeAllele)
        # topsnpDosageRoundPhenotype[topsnpDosageRoundPhenotype[, snpid] == 2, snpid] <- paste0(effectAllele, effectAllele)

        genotypeCode <- topsnpDosageRoundPhenotype[, snpid]
        topsnpDosageRoundPhenotype[, snpid] <- factor(genotypeCode, levels = sort(unique(genotypeCode), decreasing = (effectAllele > alternativeAllele)))

        countPerGenotype <- table(topsnpDosageRoundPhenotype[, snpid])
        genotypes <- sort(names(countPerGenotype), decreasing = (effectAllele > alternativeAllele))
        xTickLables <- paste0(paste(genotypes, countPerGenotype, sep = "("), ")")

        x_lable <- str_glue("Genotype (P-val: {signif(pValue, 3)}; Beta: {signif(beta, 3)}; Pos: chr{chr}:{pos}; SNP: {effectAllele}>{alternativeAllele})")
        y_lable <- str_glue("Phenotype ({phenotype})")
        g <- ggplot(topsnpDosageRoundPhenotype) + theme_bw()
        g <- g + geom_boxplot(aes_string(x = snpid, y = phenotype, color = snpid, fill = snpid), alpha = 0.5)
        g <- g + geom_point(aes_string(x = snpid, y = phenotype))
        g <- g + scale_x_discrete(labels = xTickLables)
        g <- g + labs(title = snpid)
        g <- g + xlab(x_lable)
        g <- g + ylab(y_lable)

        ggsave(paste0(phenotype, "_phenotypeLevelPerGenotype_", snpid, runFlag, ".pdf"))
    } else {
        cat("Empty QTL information for", snpid)
    }
}

#
## A transpose function exploiting melt and dcast from reshape2, which doesn't ruin the data type.
#
# TODO:
# 1. More arguments to control the transpose
myTranspose <- function(dataframe) {
    newDataframe <- reshape2::dcast(reshape2::melt(dataframe), ...~id)

    rownames(newDataframe) <- newDataframe$variable
    newDataframeColNames <- colnames(newDataframe)
    colNamesToUse <- newDataframeColNames[which(! newDataframeColNames %in% c("variable"))]
    newDataframe <- newDataframe[, colNamesToUse]

    return(newDataframe)
}

#
## Discard some columns of some data.frame
#
dscdCol <- function(dtfm, colToDiscard) {
 dtfmColNames <- colnames(dtfm)
 return(dtfm[, which(! dtfmColNames %in% colToDiscard)])
}

#
## Parsing CLI arguments
#
parser <- OptionParser(description = "A QTL mapping script based on MatrixEQTL")

parser <- add_option(parser, c("-w", "--workDir"), help = "Output direcotry.", type = "character", default = "./QTL_mapping_optdir")
parser <- add_option(parser, c("--runFlag"), help = "Running flag which help to discrimnate different runs.", type = "character", default = "RUN")

# Phenotypes related parameters
parser <- add_option(parser, c("-p", "--phenotypeFile"), help = "Phenotype input file.")
parser <- add_option(parser, c("--padding"), help = "The times of standard deviation as the boundary of outliers.", type = "integer", default = 4)
parser <- add_option(parser, c("--phenotypesToUse"), help = "Phenotypes will be used, if more than one, using comma as delimiter.", type = "character")
parser <- add_option(parser, c("--phenotypeDataFrameToTranspose"), help = "Whether should transpose the data.frame of phenotypes.", type = "logical", action = "store_false")

# Covariates related parameters
parser <- add_option(parser, c("-c", "--covariateFile"), help = "Covariates file.", type = "character")
parser <- add_option(parser, c("--covariatesToUse"), help = "Covariates will be used, if more than one, using comma as delimiter.")
parser <- add_option(parser, c("--covariateDataFrameToTranspose"), help = "Whether should transpose the data.frame of covariates.", action = "store_false")

# Phenotypes and covariates correlation
parser <- add_option(parser, c("-t", "--traitsToPairCorrelate"), help = "Traits will be correlated with in paire-wised way.")

# Genotypes related parameters
parser <- add_option(parser, c("-d", "--genotypeDosageFile"), help = "Genotype dosage file.")
parser <- add_option(parser, c("-i", "--genotypeInformationFile"), help = "Genotype information file.")
parser <- add_option(parser, c("-m", "--minorAlleleFreq"), help = "Minor allele frequency.", default = 0.05)

parser_list <- parse_args2(parser)
args <- parser_list$args
opts <- parser_list$options

# setwd("~/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AGC4nHd")
# runFlag <- "_DNAandRNA_RNAvsDNA_padding4"
workDir <- opts$workDir
runFlag <- opts$runFlag
setwd(workDir)

#
## Visualization of the correlation of phenotypes and covariates
#
# phenotypeFile <- "~/Documents/projects/200HIV/inputs/datasets/20190524_HIVreservoir_GENT.csv"
phenotypeFile <- opts$phenotypFile
phenotypes <- fread(phenotypeFile, data.table = FALSE)

if (opts$phenotypeDataFrameToTranspose) {
 phenotypes <- myTranspose(phenotypes)
}

rownames(phenotypes) <- phenotypes$id

if(is.na(opts$phenotypesToUse)) {
    phenotypesToUse <- rownames(phenotypes)
} else {
    phenotypesToUse <- str_split(phenotypesToUse, ",")[[1]]
}

# phenotypesToUse <- c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG")
phenotypesChosen <- phenotypes[, c("id", phenotypesToUse)]


# covariateFile <- "~/Documents/projects/200HIV/inputs/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv"
covariateFile <- opts$covariateFile
covariates <- fread(covariateFile, data.table = F)

if (opts$covariateDataFrameToTranspose) {
 covariates <- myTranspose(covariates)
}

rownames(covariates) <- covariates$id

if(is.na(opts$covariatesToUse)) {
    covariatesToUse <- rownames(covariates)
} else {
    covariatesToUse <- str_split(covariatesToUse, ",")[[1]]
}

covariatesChosen <- covariates[, c("id", covariatesToUse)]

# FIXME: Stopped here at 19:30, Fri. Nov. 22 2019.
phenotypesCovariatesChosen <- merge(phenotypesChosen, covariatesChosen, by = c("id"), all = TRUE)
phenotypesCovariatesChosen <- dscdCol(phenotypesCovariatesChosen, c("id"))

# phenotypesCovariatesChosenColNames <- colnames(phenotypesCovariatesChosen)
# phenotypesCovariatesChosenTrimIdColNames <- phenotypesCovariatesChosenColNames[which(!phenotypesCovariatesChosenColNames %in% c("id"))]
# phenotypesCovariatesTrimId <- phenotypesCovariatesChosen[, phenotypesCovariatesChosenTrimIdColNames]
# phenotypesCovariatesTrimId$gender <- as.factor(phenotypesCovariatesTrimId$gender)
# phenotypesCovariatesTrimId$smoking <- as.factor(phenotypesCovariatesTrimId$smoking)
# phenotypesCovariatesTrimId$HIV_DURATION <- log10(phenotypesCovariatesTrimId$HIV_DURATION)
# phenotypesCovariatesTrimId$CD4_NADIR <- log10(phenotypesCovariatesTrimId$CD4_NADIR)

# c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "HIV_DURATION", "CD4_NADIR")
traitsToCorrelate <- opts$traitsToCorrelate
if (is.na(traitsToCorrelate)){
    traitsToCorrelateVec <- rownames(phenotypesCovariatesChosen)
} else {
    traitsToCorrelateVec <- str_split(traitsToCorrelate, ",")[[1]]
}

if (dim(phenotypesCovariatesChosen)[[2]] < 16) {
    pairwisePlot <- ggpairs(phenotypesCovariatesChosen[, traitsToCorrelateVec], cardinality_threshold = 16)
    ggsave(paste0("phenotypes_covariates_pairwise", runFlag, ".pdf"), plot = pairwisePlot, width = 25, height = 25)
} else {
    warning("More than 16 variables in the data.frame, exceeding cardinality_threshold")
}

# TODO: add a section to draw Spearman's rank correlation matrix???

# Remove outliers of phenotypes
# padding <- 4
padding <- opts$padding
phenotypesMeans <- sapply(as.data.frame(phenotypesChosen[, phenotypesToUse]), FUN = mean, na.rm = T)
phenotypesStdev <- sapply(as.data.frame(phenotypesChosen[, phenotypesToUse]), FUN = sd, na.rm = T)

upperBound <- phenotypesMeans + padding * phenotypesStdev
lowerBound <- phenotypesMeans - padding * phenotypesStdev

names(upperBound) <- phenotypesToUse
names(lowerBound) <- phenotypesToUse

phenotypeLogRownames <- rownames(phenotypesChosen)
for (phenotype in phenotypesToUse) {
    outlierRowIndex <- which((lowerBound[phenotype] > phenotypesChosen[, phenotype]) | (phenotypesChosen[, phenotype] > upperBound[phenotype]))
    phenotypesChosen[outlierRowIndex, phenotype] <- NA
    cat("For ", phenotype, " the outlier index out of ", padding, " * SD boundary [", lowerBound[phenotype], ", ", upperBound[phenotype], "]: \n ", sep = "")
    cat(phenotypeLogRownames[outlierRowIndex], "\n")
}

#
## Phenotypes
#
# Remove id to prevent the strings ruin the matrix transform
usedPhenotypeColNames <- phenotypeLogRownames[which(! phenotypeLogRownames %in% c("id"))]
phenotypeForMe <- as.data.frame(t(phenotypesChosen[, usedPhenotypeColNames]))

rownames(phenotypeForMe) <- phenotypesToUse
colnames(phenotypeForMe) <- rownames(phenotypesChosen)

cat("phenotypeForMe dim:", dim(phenotypeForMe), "\n")

# Preprocessing covariates
covariatesColNames <- colnames(covariates)
usedCovariateColNames <- covariatesColNames[which(! covariatesColNames %in% c("id"))]

covariatesTransposed <- as.data.frame(t(covariates[, usedCovariateColNames]))

rownames(covariatesTransposed) <- colnames(covariates)
colnames(covariatesTransposed) <- rownames(covariates)

#
## Covariates
#
covariateForME <- covariatesTransposed[c("age", "gender", "HIV_DURATION", "CD4_NADIR"), ]
cat("covariateForME dim:", dim(covariateForME), "\n")

# Preprocessing genotypes
# genotypeDosageFile <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/chr11_200HIV_dosage.gz" # For debugging
# genotypeDosageFile <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/chr8_200HIV_dosage.gz" # For debugging
# genotypeDosageFile <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/200HIV_dosage.gz"
genotypeDosageFile <- opts$genotypeDosageFile
genotypeDosage <- fread(genotypeDosageFile, header = TRUE, data.table = FALSE)

rownames(genotypeDosage) <- genotypeDosage$id

genotypeDosageColNames <- colnames(genotypeDosage)
usedGenotypeDosageColNames <- genotypeDosageColNames[which(! genotypeDosageColNames %in% c("id"))]
genotypeDosage <- genotypeDosage[, usedGenotypeDosageColNames] # remove id

minorAlleleFreq <- opts$minorAlleleFreq
nSamples <- dim(genotypeDosage)[[2]]
genotypeDosageRound <- round(genotypeDosage) # For the phenotype level per genotype plot
genotypeAlleleFrequency <- apply(genotypeDosageRound, 1, sum) / (2 * nSamples)
normalVariants <- which((genotypeAlleleFrequency >= minorAlleleFreq) & (genotypeAlleleFrequency <= 1 - minorAlleleFreq))

#
## Geontypes
#
genotypeForME <- genotypeDosage[normalVariants, ]
cat("genotypeForME dim:", dim(genotypeForME), "\n")

phenotypeSamples <- colnames(phenotypeForMe)
covariatesSamples <- colnames(covariateForME)
genotypesSamples <- colnames(genotypeForME)

commonSamples <- intersect(phenotypeSamples, intersect(covariates_samples, genotypes_samples))
nSharedSample <- length(commonSamples)
cat("Samples (", nSharedSample, ") will be fed to MatrixEQTL:", commonSamples, "\n")

phenotypeMatrixForME <- as.matrix(phenotypeForMe[, commonSamples])
fwrite(phenotypeForMe[, commonSamples], file = paste0("phenotypes_fedtoMatrixEQTL", runFlag, ".tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

covariateMatrixForME <- as.matrix(covariateForME[, commonSamples])
fwrite(covariateForME[, commonSamples], file = paste0("covariates_fedtoMatrixEQTL", runFlag, ".tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

genotypeMatrixForME <- as.matrix(genotypeForME[, commonSamples])

#
## QTL mapping
#
snps <- SlicedData$new()
snps$fileOmitCharacters <- "NA" # denote missing values;
snps$fileSliceSize <- 2000;
snps$CreateFromMatrix(genotypeMatrixForME)

cvrt <- SlicedData$new()
cvrt$fileOmitCharacters <- "NA" # denote missing values;
cvrt$fileSliceSize <- 2000;
cvrt$CreateFromMatrix(covariateMatrixForME)

phenotypeMatrixForMeForPermutation <- phenotypeMatrixForME # Matrix for the permutations operation
useModel <- modelLINEAR
pvOutputThreshold <- 0.999
errorCovariance <- numeric()

permutationTimes <- 0
seed <- 1234
set.seed(seed)

for (pm in 0:permutationTimes) {
    gene <- SlicedData$new()
    gene$fileOmitCharacters <- "NA" # denote missing values;
    gene$fileSliceSize <- 2000

    if (pm) {
        cat("Permuted genotype:", pm, "\n")
        outputFileName <- paste0("qtls_hivresorvior", runFlag, "_permut", pm, ".tsv")
        phenotypeMatrixForMeForPermutation <- t(apply(phenotypeMatrixForMeForPermutation, 1, sample, size = nSharedSample, replace=TRUE))

        gene$CreateFromMatrix(phenotypeMatrixForMeForPermutation)
    } else {
        cat("Raw phenotype\n")
        outputFileName <- paste0("qtls_hivresorvior", runFlag, ".tsv")
        gene$CreateFromMatrix(phenotypeMatrixForME)
    }

    me = Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = outputFileName,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    )
}

#
## Report
#
qtls <- fread(outputFileName, data.table = FALSE)

genotypeInformationFile <- opts$genotypeInformationFile
genotypeInformation <- fread(genotypeInformationFile, data.table = FALSE)

qtlsInformation <- merge(qtls, genotypeInformation, by.x = "SNP", by.y = "rsID")

phenotypeLevelPerGenotypeThreshold <- 5e-7 # The threshold of p-value for the plot of phenotype level per genotype
manhattanPlotThreshold <- 1e-2

# phenotypeLevelPerGenotypeSnpIds <- c() # Extra SNPs to draw the plot of phenotype level per genotype 
phenotypeLevelPerGenotypeSnpIds <- c("rs7113204", "rs2613996", "rs7817589") # Extra SNPs to draw the plot of phenotype level per genotype 

for (phenotype in phenotypesToUse) {
    qtlmappingResults <- qtlsInformation[qtlsInformation$gene == phenotype, c("SNP", "SequenceName", "Position", "p-value")]

    if (nrow(qtlmappingResults) == 0) {
        next
    }

    pdf(paste0(phenotype, "_MahattanPlot", runFlag, ".pdf"), width = 16, height = 9)
    colnames(qtlmappingResults) <- c("SNP", "CHR", "BP", "P")
    manhattan(
        qtlmappingResults[qtlmappingResults["P"] <= manhattanPlotThreshold, ], 
        main = paste("Manhattan plot", phenotype), ylab = "p-value(-log10)", annotateTop = TRUE
    )
    dev.off()

    pdf(paste0(phenotype, "_QQPlot", runFlag, ".pdf"), width = 16, height = 16)
    qq(qtlmappingResults$P, main = paste("Q-Q plot", phenotype), pch = 18, cex = 1, las = 1)
    dev.off()

    phenotypeLevels <- phenotypesForMe[phenotype, commonSamples]
    chroms <- c(unique(qtlmappingResults$CHR))

    for (chr in chroms) {
        qtlmappingResultsChr <- qtlmappingResults[qtlmappingResults$CHR == chr, ] 
        minP <- min(qtlmappingResultsChr[, "P"])
        if (minP > phenotypeLevelPerGenotypeThreshold) { next }
        phenotypeLevelPerGenotypeSnpIds <- c(phenotypeLevelPerGenotypeSnpIds, qtlmappingResultsChr[qtlmappingResultsChr$P == minP, "SNP"])
    }

    phenotypeLevelPerGenotypeSnpIds <- unique(phenotypeLevelPerGenotypeSnpIds)
    for (snpid in phenotypeLevelPerGenotypeSnpIds) {
        print(snpid)
        pppg(
             snpid = snpid,
             phenotype = phenotype,
             qtlsInformation = qtlsInformation,
             phenotypeLevels = phenotypeLevels,
             commonSamples = commonSamples,
             genotypeDosageRound = genotypeDosageRound
        )
    }
}
