#!/usr/bin/env Rscript
require(stringr)
require(ggplot2)
require(qqman)
require(data.table)
require(MatrixEQTL)

pppg <- function(snpid, phenotype, qtlsInformation, phenotypeLevels, commonSamples, genotypeDosageRound) {
    condition <- qtlsInformation$SNP == snpid & qtlsInformation$gene == phenotype
    currentQtlInfo <- qtlsInformation[condition, ]

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

        topsnpDosageRoundPhenotype[topsnpDosageRoundPhenotype[, snpid] == 0, snpid] <- paste0(alternativeAllele, alternativeAllele)
        topsnpDosageRoundPhenotype[topsnpDosageRoundPhenotype[, snpid] == 1, snpid] <- paste0(effectAllele, alternativeAllele)
        topsnpDosageRoundPhenotype[topsnpDosageRoundPhenotype[, snpid] == 2, snpid] <- paste0(effectAllele, effectAllele)

        genotypeCode <- topsnpDosageRoundPhenotype[, snpid]
        topsnpDosageRoundPhenotype[, snpid] <- factor(genotypeCode, levels = sort(unique(genotypeCode), decreasing = (effectAllele > alternativeAllele)))

        countPerGenotype <- table(topsnpDosageRoundPhenotype[, snpid])
        genotypes <- sort(names(countPerGenotype), decreasing = (effectAllele > alternativeAllele))
        xTickLables <- paste0(paste(genotypes, countPerGenotype, sep="("), ")")

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


setwd("~/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AGC4nHd")
runFlag <- "_DNAandRNA_RNAvsDNA_padding4"

# Visualization of the correlation of phenotypes and covariates
phenotypesFile <- "~/Documents/projects/200HIV/inputs/datasets/20190524_HIVreservoir_GENT.csv"
phenotypes <- fread(phenotypesFile, data.table = F)

phenotypes["RNAvsDNA_CD4LOG"] <- phenotypes["RNAHIV_CD4LOG"] - phenotypes["DNAHIV_CD4LOG"]
phenotypes["DNAvsRNA_CD4LOG"] <- phenotypes["DNAHIV_CD4LOG"] - phenotypes["RNAHIV_CD4LOG"]

phenotypeNames <- c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG") #, "DNAvsRNA_CD4LOG")
phenotypeLog <- phenotypes[, c("id", phenotypeNames)]

rownames(phenotypeLog) <- phenotypeLog$id

covariatesFile <- "~/Documents/projects/200HIV/inputs/datasets/metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv"
covariates <- fread(covariatesFile, data.table = F)
rownames(covariates) <- covariates$id

phenotypesCovariates <- merge(phenotypeLog, covariates, by=c("id"), all=T)
rownames(phenotypesCovariates) <- phenotypesCovariates$id
phenotypesCovariatesTrimId <- phenotypesCovariates[, -1]
phenotypesCovariatesTrimId$gender <- as.factor(phenotypesCovariatesTrimId$gender)
phenotypesCovariatesTrimId$smoking <- as.factor(phenotypesCovariatesTrimId$smoking)
phenotypesCovariatesTrimId$HIV_DURATION <- log10(phenotypesCovariatesTrimId$HIV_DURATION)
phenotypesCovariatesTrimId$CD4_NADIR <- log10(phenotypesCovariatesTrimId$CD4_NADIR)

pairwisePlot <- ggpairs(phenotypesCovariatesTrimId[, c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "HIV_DURATION", "CD4_NADIR")], cardinality_threshold = 16)
ggsave(paste0("phenotypes_covariates_pairwise", runFlag, ".pdf"), plot = pairwisePlot, width = 25, height = 25)

# Remove outliers of phenotypes
padding <- 4
phenotypesMeans <- sapply(as.data.frame(phenotypeLog[, phenotypeNames]), FUN = mean, na.rm = T)
phenotypesSd <- sapply(as.data.frame(phenotypeLog[, phenotypeNames]), FUN = sd, na.rm = T)

upperBound <- phenotypesMeans + padding * phenotypesSd
lowerBound <- phenotypesMeans - padding * phenotypesSd

names(upperBound) <- phenotypeNames
names(lowerBound) <- phenotypeNames

phenotypeLogRownames <- rownames(phenotypeLog)
for (targetPhenotype in phenotypeNames) {
    outlierRowIndex <- which((lowerBound[targetPhenotype] > phenotypeLog[, targetPhenotype]) | (phenotypeLog[, targetPhenotype] > upperBound[targetPhenotype]))
    phenotypeLog[outlierRowIndex, targetPhenotype] <- NA
    cat("For ", targetPhenotype, " the outlier index out of ", padding, " * SD boundary [", lowerBound[targetPhenotype], ", ", upperBound[targetPhenotype], "]: \n  ", sep = "")
    cat(phenotypeLogRownames[outlierRowIndex], "\n")
}

#
## Phenotypes
#
# Remove id to prevent the strings ruin the matrix transform
phenotypesForMe <- as.data.frame(t(phenotypeLog[, -1]))
cat("phenotypeForMe dim:", dim(phenotypeForMe), "\n")

rownames(phenotypeForMe) <- phenotypeNames
colnames(phenotypeForMe) <- rownames(phenotypeLog)

# Preprocessing covariates
covariatesTransposed <- as.data.frame(t(covariates[, -1]))
rownames(covariatesTransposed) <- colnames(covariates[, -1])
colnames(covariatesTransposed) <- rownames(covariates)

#
## Covariates
#
covariateForME <- covariatesTransposed[c("age", "gender", "HIV_DURATION", "CD4_NADIR"), ]
cat("covariateForME dim:", dim(covariateForME), "\n")

# Preprocessing genotypes
# genotypeFile <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/chr11_200HIV_dosage.gz"  # For debugging
# genotypeFile <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/chr18_200HIV_dosage.gz"  # For debugging
genotypeFile <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/200HIV_dosage.gz"
genotypeDosage <- fread(genotypeFile, header = TRUE, data.table = FALSE)

rownames(genotypeDosage) <- genotypeDosage$id
genotypeDosage <- genotypeDosage[, -1]  # remove id
nSamples <- dim(genotypeDosage)[[2]]

mafThreshold <- 0.05
genotypeDosageRound <- round(genotypeDosage)  # For the phenotype level per genotype plot
genotypeAlleleFrequency <- apply(genotypeDosageRound, 1, sum) / (2 * nSamples)
normalVariants <- which((genotypeAlleleFrequency >= mafThreshold) & (genotypeAlleleFrequency <= 1 - mafThreshold))

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
cat("Samples (", length(commonSamples), ") will be fed to MatrixEQTL:", commonSamples, "\n")

phenotypeMatrixForME <- as.matrix(phenotypeForMe[, commonSamples])
fwrite(phenotypeForMe[, commonSamples], file = paste0("phenotypes_fedtoMatrixEQTL", runFlag, ".tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

covariateMatrixForME <- as.matrix(covariateForME[, commonSamples])
fwrite(covariateForME[, commonSamples], file = paste0("covariates_fedtoMatrixEQTL", runFlag, ".tsv"), sep = "\t", row.names = TRUE, col.names = TRUE)

genotypeMatrixForME <- as.matrix(genotypeForME[, commonSamples])

#
## QTL mapping
#
snps <- SlicedData$new()
# snps$fileDelimiter <- "\t"    # the TAB character
snps$fileOmitCharacters <- "NA" # denote missing values;
# snps$fileSkipRows <- 1
# snps$fileSkipColumns <- 1
snps$fileSliceSize <- 2000;
snps$CreateFromMatrix(genotypeMatrixForME)

cvrt <- SlicedData$new()
# cvrt$fileDelimiter <- "\t"    # the TAB character
cvrt$fileOmitCharacters <- "NA" # denote missing values;
# cvrt$fileSkipRows <- 1      # one row of column labels
# cvrt$fileSkipColumns <- 1     # one column of row labels
cvrt$CreateFromMatrix(covariateMatrixForME)

phenotypeMatrixForMeForPermutation <- phenotypeMatrixForME  # Matrix for the permutations operation
useModel <- modelLINEAR
pvOutputThreshold <- 0.999
errorCovariance <- numeric()

permutationTimes <- 0
seed <- 1234
set.seed(seed)

for (pm in 0:permutationTimes) {
    gene <- SlicedData$new()
    # gene$fileDelimiter <- "\t"
    gene$fileOmitCharacters <- "NA" # denote missing values;
    # gene$fileSkipRows <- 1
    # gene$fileSkipColumns <- 1
    gene$fileSliceSize <- 2000

    if (pm) {
        cat("Permuted genotype:", pm, "\n")
        outputFileName <- paste0("qtls_hivresorvior", runFlag, "_permut", pm, ".tsv")
        phenotypeMatrixForMeForPermutation = t(apply(phenotypeMatrixForMeForPermutation, 1, sample, size = nSharedSample, replace=TRUE))

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
qtlsFile <- paste0("~/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AGC4nHd/qtls_hivresorvior", runFlag, ".tsv")
qtls <- fread(qtlsFile, sep = "\t", data.table = FALSE)

information_file <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/variantInfo.gz"
information <- fread(information_file, sep = " ", data.table = FALSE)

qtlsInformation <- merge(qtls, information, by.x = "SNP", by.y = "rsID")

phenotypeLevelPerGenotypeThreshold <- 5e-7  # The threshold of p-value for the plot of phenotype level per genotype
manhattanPlotThreshold <- 1e-2

# phenotypeLevelPerGenotypeSnpIds <- c()  # Extra SNPs to draw the plot of phenotype level per genotype 
phenotypeLevelPerGenotypeSnpIds <- c("rs7113204", "rs2613996", "rs7817589")  # Extra SNPs to draw the plot of phenotype level per genotype 

for (phenotype in phenotypeNames) {
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
