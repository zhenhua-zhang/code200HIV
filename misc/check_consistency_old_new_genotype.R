#!/usr/bin/env Rscript

library(ggplot2)
library(stringi)
library(optparse)
library(data.table)

parser <- OptionParser()

parser <- add_option(
     parser, c("-o", "--old-qtls"), action = "store", type = "character",
     dest = "old_qtls", default = "old_qtls.txt",
     help = "The file include old QTLs."
)

parser <- add_option(
     parser, c("-n", "--new-qtls"), action = "store", type = "character",
     dest = "new_qtls", default = "new_qtls.txt",
     help = "The file include new QTLs."
)

parser <- add_option(
     parser, c("-p", "--plot-name"), action = "store", type = "character",
     dest = "plot_name", default = NULL,
     help = "The name of plot. Defalt: [old-qtls-name]_vs_[new-qtls-name].pdf"
)


opts_args <- parse_args2(parser)
opts <- opts_args$options

old_qtls <- opts$old_qtls
old_qtls_dtfm <- fread(
    old_qtls, header = FALSE, data.table = FALSE, stringsAsFactors = FALSE,
    verbose = FALSE
)
colnames(old_qtls_dtfm) <- c("rs", "pval")

new_qtls <- opts$new_qtls
new_qtls_dtfm <- fread(
    new_qtls, header = FALSE, data.table = FALSE, stringsAsFactors = FALSE,
    verbose = FALSE
)
colnames(new_qtls_dtfm) <- c("rs", "pval")

merged_dtfm <- merge(
    old_qtls_dtfm, new_qtls_dtfm, by = "rs", all = TRUE,
    suffixes = c("_old", "_new")
)

merged_dtfm[is.na(merged_dtfm)] <- max(merged_dtfm[, "pval_old"], na.rm = TRUE)

plot_name <- opts_args$plot_name

if (is.null(plot_name)) {
    old_qtls_name <- stri_replace(basename(old_qtls), replacement = "", fixed = ".txt")
    new_qtls_name <- stri_replace(basename(new_qtls), replacement = "", fixed = ".txt")
    plot_name <- paste0(old_qtls_name, "_vs_", new_qtls_name, ".png")
    title <- new_qtls_name
}

g <- ggplot(data = merged_dtfm, aes(x = -log10(pval_old), y = -log10(pval_new))) + theme_bw()
g <- g + geom_point()
g <- g + labs(title = title)

ggsave(plot_name, plot = g)
