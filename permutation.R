#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(gtools)

# The effect allele is encoded as 1

# TODO: gene-gene intercation analaysis https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001338
#
# TODO: a project


snpid <- "rs7113204"
snpid <- "rs2613996"
snpid <- "rs7817589"

# snpid <- "rs7926172"
# snpid <- "rs7936906"
#
## Dosages
#
genotype_file <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/chr8_200HIV_dosage.gz"  # For debugging
genotypes <- fread(genotype_file, data.table = FALSE)
rownames(genotypes) <- genotypes$id
genotypes <- genotypes[, -1]

snp_dosage <- as.data.frame(t(genotypes[snpid, ]))
names(snp_dosage) <- "dosage"

snp_dosage_round <- round(snp_dosage)
names(snp_dosage_round) <- "genotype"
snp_genotype <- cbind.data.frame(snp_dosage, snp_dosage_round)


#
## Genotype information
#
genotype_info_file <- "~/Documents/projects/200HIV/inputs/dosage/200HIV_dosages/chr8.variantInfo.gz"
genotype_info <- fread(genotype_info_file, data.table = FALSE, sep=" ")
rownames(genotype_info) <- genotype_info$rsID
genotype_info <- genotype_info[, -1]

snp_info <- genotype_info[snpid, ]
ea <- snp_info$EffectAllele  # Effect allele encode as 1
aa <- snp_info$AlternativeAllele  # Alternative allele as 0

#
## Phenotypes
#
phenotypes_file <- "~/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AGC4nHd/phenotypes_fedtoMatrixEQTL_DNAandRNA_RNAvsDNA_padding4.tsv"
phenotypes <- fread(phenotypes_file, data.table = FALSE)
rownames(phenotypes) <- phenotypes$V1
phenotypes <- as.data.frame(t(phenotypes[, -1]))

#
## Covariates
#
covariates_file <- "~/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AGC4nHd/covariates_fedtoMatrixEQTL_DNAandRNA_RNAvsDNA_padding4.tsv"
covariates <- fread(covariates_file, data.table = FALSE)
rownames(covariates) <- covariates$V1
covariates <- as.data.frame(t(covariates[, -1]))

common_samples <- intersect(intersect(rownames(phenotypes), rownames(covariates)), rownames(snp_dosage))
data_dataframe <- cbind.data.frame(phenotypes[common_samples, ], snp_genotype[common_samples, ], covariates[common_samples, ])

eg <- paste0(ea, ea) # Effect genotype, AA
hg <- paste0(ea, aa) # Heterzygous genotype, Aa
ag <- paste0(aa, aa) # Alternative genotype, aa

data_dataframe[data_dataframe$genotype == 0, "genotype"] <- ag
data_dataframe[data_dataframe$genotype == 1, "genotype"] <- hg
data_dataframe[data_dataframe$genotype == 2, "genotype"] <- eg

gc <- as.data.frame(table(data_dataframe$genotype))  # Number per genotype, genotype counts
rownames(gc) <- gc$Var1

obaan <- gc[ag, "Freq"]  # Observed alternative genotype number
obean <- gc[eg, "Freq"]  # observed effect genotype number
obhtn <- gc[hg, "Freq"]  # observed heter-genotype number

gt <- matrix(c(ifelse(is.na(obaan), 0, obaan), obhtn / 2, obhtn / 2, obean), ncol = 2)
ct <- chisq.test(gt)
ft <- fisher.test(gt)


m1 <- lm(RNAvsDNA_CD4LOG~dosage + age + gender + HIV_DURATION + CD4_NADIR, data=data_dataframe)
summary(m1)

m2 <- glm(RNAHIV_CD4LOG~genotype + age + gender + HIV_DURATION + CD4_NADIR, data=data_dataframe)
summary(m2)

m3 <- glm(DNAHIV_CD4LOG~genotype + age + gender + HIV_DURATION + CD4_NADIR, data=data_dataframe)
summary(m3)

snp_statistic <- cor(data_dataframe$RNAvsDNA_CD4LOG, data_dataframe$dosage, use = "complete.obs", method = "spearman")
snp_statistic <- summary(lm(RNAvsDNA_CD4LOG~dosage + age + gender + HIV_DURATION + CD4_NADIR, data=data_dataframe))$coefficients["dosage", "t value"]

statistic_pool <- c()
size <- length(common_samples)
permut_dtfm <- data_dataframe

phenotype_names <- c("DNAHIV_CD4LOG", "DNAHIV_CD4LOG", "RNAvsDNA_CD4LOG")
seed <- 1234
set.seed(seed)

for (x in 1:1000) {
  # permut_dosage <- sample(data_dataframe[, "dosage"], size, replace = TRUE)
  # statistic_pool <- c(statistic_pool, cor(data_dataframe$RNAvsDNA_CD4LOG, permut_dosage, use = "complete.obs", method = "spearman"))
  for (phtp in phenotype_names) {
    permut_dtfm[, phtp] <- sample(data_dataframe[, phtp], size, replace = TRUE)
  }
  
  statistic_pool <- c(statistic_pool, summary(lm(RNAvsDNA_CD4LOG~dosage + age + gender + HIV_DURATION + CD4_NADIR, data=permut_dtfm))$coefficients["dosage", "t value"])
}


g <- ggplot() + theme_bw()
g <- g + geom_histogram(aes(x=statistic_pool))
g <- g + geom_vline(xintercept = snp_statistic, color = "red")
g <- g + xlab("T-statistic") + ylab("Counts")
g

pnorm(snp_statistic, mean=mean(statistic_pool), sd=sd(statistic_pool))
