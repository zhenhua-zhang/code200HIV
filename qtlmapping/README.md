# `qtlmap` instructions

## Objective
Scripts to do QTL mapping using `MatrixEQTL` package in R.

## Introduction
Data is preprocessed using `numpy` + `pandas`. Mapping results are visualized
using `assocplots` which was published in 2017 on Bioinformatics.

## Usage


## I/O
Input:
    - Phenotypes
    - Genotypes
    - Covariates (e.g. age / gender)
    - Gene features (e.g. GFF)

Output:
    - Manhattan plot
    - Q-Q plot
    - Report for QTL mapping results

## Functions
Preprocessing
    - Check the normality of the phenotype data.
    - Transformation using log10 / inverse-rank / log2 etc. if not normal distributed.
    - Remove outliers using PCA or $mean -/+ n * sd$

QTL-mapping
    - Logistic regression
    - P-value adjustment??

Visualization
    - Manhattan plot
    - Q-Q plot
    - Dosage boxplot (optionally adjusted by age)
    - *Normality check*
    - *Normality (histgrams) after transformation*

Generate mapping report
    - Locus-zoom plot (manually)
    - Genes in LD (linkage disequilibrium)

Bonus analysis
- Mendelian randomization
- Pathway enrichment analysis

## Misc
### Phenotype file format
    | Phenotype_1 | Phenotype_2 | ... | Phenotype_n |
 id1|             |             |     |             |
 id2|             |             |     |             |
 ...|             |             |     |             |
 idn|             |             |     |             |