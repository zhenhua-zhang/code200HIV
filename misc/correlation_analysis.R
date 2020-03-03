#!/usr/bin/env Rscript
# A script to analyze the correlation between phenotypes and covariates

#TODO:
#NOTE:

library(ggplot2)
library(stringr)
library(corrplot)
library(data.table)

setwd("/home/umcg-zzhang/Documents/projects/200HIV/outputs/HIVReservior/HIVReservior_AGC4nHd/correlation_analysis/")
#
## Utils
#
read_merge <- function(file_root, file_path_vec = NULL, merge_on = "id") {
  if (is.null(file_path_vec)) {
    file_path_vec <- c(file_root)
  } else {
    file_path_vec <- paste0(file_root, "/", file_path_vec)
  }

  dtfm <- NULL
  for(file_path in file_path_vec) {
    tmp_dtfm <- NULL
    
    if (is.null(dtfm)){
      dtfm <- fread(file_path, verbose = F, data.table = F, stringsAsFactors = F)
    } else {
      tmp_dtfm <- fread(file_path, verbose = F, data.table = F, stringsAsFactors = F)
      dtfm <- merge(dtfm, tmp_dtfm, by = "id")
    }
  }

  return(dtfm)
}

#
## Correlation matrix
#
cor_matrix <- function(data, candidates, new_name, method = "spearman", plot_name=NULL) {
  candicate_dtfm <- data[candidates] 
  colnames(candicate_dtfm) <- new_name

  # Correlation matrix
  if (is.null(plot_name)) {
    plot_name <- paste("corr_matrix", method, "heatmap.png", sep = "_")
  }

  corr <- cor(candicate_dtfm, use = "complete.obs", method = method)

  png(plot_name, width = 1080, height = 720)
  corrplot(corr, type = "upper", order = "original", method = "number")
  dev.off()
}

#
## Correlation analysis between DNAHIV_CD4lOG and cell counts
#
model_fit <- function(response, target_predictor, data, cvrt_var_vec = NULL, plot_name = NULL, report_name = NULL) {
  if (!is.null(cvrt_var_vec)) {
    predictor_vec <- c(target_predictor, cvrt_var_vec)
  } else {
    predictor_vec <- target_predictor
  }
  
  fml <- formula(paste(response, "~", paste(predictor_vec, collapse = " + ")))
  m <- glm(formula = fml, data = data)
  
  dev_total <- with(summary(m), 1 - deviance / null.deviance)
  
  m_anova <- anova(m, test = "F")
  dev_trgt_var_per <- m_anova[target_predictor, "Deviance"] / max(m_anova$`Resid. Dev`, na.rm = T)
  names(dev_trgt_var_per) <- target_predictor
  
  if (!is.null(cvrt_var_vec)) {
    dev_cvrt_var_per <- m_anova[cvrt_var_vec, "Deviance"] / max(m_anova$`Resid. Dev`, na.rm = T)
    names(dev_cvrt_var_per) <- cvrt_var_vec
  }
  
  g <- ggplot(data = data, aes_string(y = response, x = target_predictor)) + theme_bw()
  g <- g + geom_point()
  g <- g + geom_smooth(method="lm")

  if (is.null(plot_name)) {
    plot_name <- paste0(response, "_vs_", target_predictor, ".png") 
  }
  ggsave(plot_name, g)
  
  if (is.null(report_name)) {
    report_name <- paste0(response, "_vs_", target_predictor, ".txt") 
  }
  capture.output(list(model = m, fromula = fml, anova = m_anova, dev_total = dev_total, dev_trgt_var_per = dev_trgt_var_per, dev_cvrt_var_per = dev_cvrt_var_per), file = report_name)
}


file_root <- "~/Documents/projects/200HIV/inputs/datasets"
file_vec <- c(
  "20190524_HIVreservoir_GENT_withRNADNARatio.tsv",
  "metaData_pcntgMnct_ssnlt_CD4CD8TC_CD4NADIR_HIVDuration_20190728.csv",
  "totalDataHans200HivWithPercent.csv"
)

gntp_dtfm <- fread("../../../../inputs/datasets/20190524_HIVreservoir_GENT_withRNADNARatio.tsv", data.table = F, verbose = F, header = T)
sample_names <- row.names(gntp_dtfm)

RNAHIV_LOG <- as.vector(gntp_dtfm$RNAHIV_CD4LOG)
RNAHIV_LOG_samples <- gntp_dtfm$id[is.na(RNAHIV_LOG)]

DNAHIV_LOG <- as.vector(gntp_dtfm$DNAHIV_CD4LOG)
DNAHIV_LOG_samples <- gntp_dtfm$id[is.na(DNAHIV_LOG)]

RNA2DNA_LOG <- as.vector(gntp_dtfm$RNAvsDNA_CD4LOG)
RNA2DNA_LOG_samples <- gntp_dtfm$id[is.na(RNA2DNA_LOG)]

integrate_dtfm <- read_merge(file_root, file_vec, merge_on = "id")
integrate_dtfm["CD4_vs_CD8_ratio"] <- integrate_dtfm$CD4P_T_cells_LMI1_parentPerc__Index8 / integrate_dtfm$CD8P_T_cells_LMI1_parentPerc__Index8

# Correlation by Spearman's R
var_vec <- c("DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "age", "HIV_DURATION", "CD4_NADIR")#, "CD4P_T_cells_LMI1_parentPerc__Index8", "CD8P_T_cells_LMI1_parentPerc__Index8", "CD4_vs_CD8_ratio", "Prol_CD8_PBMC_LMI4_parentPerc__Index81", "Prol_CD4P_Tconv_PBMC_LMI4_parentPerc__Index82", "Prol_CD4P_Treg_PBMC_LMI4_parentPerc__Index83", "Prol_DN_CD4NCD8N_PBMC_LMI4_parentPerc__Index80", "Prol_DP_CD4PCD8P_PBMC_LMI4_parentPerc__Index79") 
new_name <- c( "DNAHIV_CD4LOG", "RNAHIV_CD4LOG", "RNAvsDNA_CD4LOG", "age", "HIV_DURATION", "CD4_NADIR")#, "CD4P", "CD8P", "CD4_vs_CD8_ratio", "Prol_CD8", "Prol_CD4P_Tconv", "Prol_CD4P_Treg", "Prol_DN_CD4NCD8N", "Prol_DP_CD4PCD8P" )
cor_matrix(data = integrate_dtfm, candidates = var_vec, new_name = new_name)

# Correlation by linear regression
model_fit(response = "RNAHIV_CD4LOG", target_predictor = "DNAHIV_CD4LOG", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "RNAvsDNA_CD4LOG", target_predictor = "DNAHIV_CD4LOG", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "RNAvsDNA_CD4LOG", target_predictor = "RNAHIV_CD4LOG", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))

model_fit(response = "RNAHIV_CD4LOG", target_predictor = "CD4P_T_cells_LMI1_parentPerc__Index8", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "DNAHIV_CD4LOG", target_predictor = "CD4P_T_cells_LMI1_parentPerc__Index8", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "RNAvsDNA_CD4LOG", target_predictor = "CD4P_T_cells_LMI1_parentPerc__Index8", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))

model_fit(response = "RNAHIV_CD4LOG", target_predictor = "CD8P_T_cells_LMI1_parentPerc__Index8", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "DNAHIV_CD4LOG", target_predictor = "CD8P_T_cells_LMI1_parentPerc__Index8", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "RNAvsDNA_CD4LOG", target_predictor = "CD8P_T_cells_LMI1_parentPerc__Index8", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))
model_fit(response = "RNAvsDNA_CD4LOG", target_predictor = "Prol_CD8_PBMC_LMI4_parentPerc__Index81", data = integrate_dtfm, cvrt_var_vec = c("age", "gender", "HIV_DURATION", "CD4_NADIR"))

target_variates <- c(
  "age", "gender", "HIV_DURATION", "CD4_NADIR",
  "CD4P_T_cells_LMI1_parentPerc__Index8",
  "CD8P_T_cells_LMI1_parentPerc__Index8",
  "CD4_vs_CD8_ratio",
  "Prol_CD8_PBMC_LMI4_parentPerc__Index81",
  "Prol_CD4P_Tconv_PBMC_LMI4_parentPerc__Index82",
  "Prol_CD4P_Treg_PBMC_LMI4_parentPerc__Index83",
  "Prol_DN_CD4NCD8N_PBMC_LMI4_parentPerc__Index80",
  "Prol_DP_CD4PCD8P_PBMC_LMI4_parentPerc__Index79"
)


#
## All types of cell counts
#
all_variates <- colnames(integrate_dtfm)
target_variates <- all_variates[grepl("_parentPerc__", all_variates)]
target_variates <- c("age", "gender", "HIV_DURATION", "CD4_NADIR", target_variates)

fml <- formula(paste("RNAvsDNA_CD4LOG", "~", paste0(target_variates, collapse = " + ")))
m <- glm(fml, data = integrate_dtfm)
summary(m)
plot(m$residuals)

m1 <- step(m)
summary(m1)
plot(m1$residuals)

anova(m, m1, test = "Chi")

target_predictor <- "CD8P_T_cells_LMI1_parentPerc__Index8"
model_fit(
  response = "RNAvsDNA_CD4LOG",
  target_predictor = target_predictor,
  data = integrate_dtfm,
  cvrt_var_vec = target_variates[!target_variates %in% c(target_predictor)]
)

cor.test(integrate_dtfm$HIV_DURATION, integrate_dtfm$DNAHIV_CD4, method = "spearman")
for (pntp in c("RNAvsDNA_CD4LOG", "HIVRNA_CD4LOG", "HIVDNA_CD4LOG")) {
  for (cvrt in c("age", "HIV_DURATION", "CD4_NADIR")) {
    cor_val <- cor.test(integrate_dtfm[pntp], integrate_dtfm[cvrt], method = "spearman")
    g <- ggplot(data = integrate_dtfm) + theme_bw()
    g <- g + geom_point(aes_string(x = cvrt, y = pntp)) + geom_smooth(aes(x = cvrt, y = pntp), method = "lm")
    ggsave(str_glue("{pntp}_vs_{cvrt}.pdf"), width = 10, height = 10)
  }
}
