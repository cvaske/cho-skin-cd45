#!/usr/bin/env Rscript

library(Seurat)
library(harmony)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

input.h5 <- args[1]
basename <- args[2]

data <- Read10X_h5(input.h5)

data.s <- CreateSeuratObject(counts = data, project = basename, min.cells = 3, min.features = 200)

data.s

data.s[["percent.mt"]] <- PercentageFeatureSet(data.s, pattern = "^MT-")

pdf(paste(basename, ".pdf", sep = ""))
ggplot(data.s@meta.data, aes(nFeature_RNA)) +
    geom_histogram(binwidth=200) +
    ggtitle(basename)

ggplot(data.s@meta.data, aes(nCount_RNA)) +
    geom_histogram(binwidth=500) +
    ggtitle(basename)

ggplot(data.s@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_bin2d(bins=200) +
    ggtitle(basename)

ggplot(data.s@meta.data, aes(percent.mt, nFeature_RNA)) +
    geom_bin2d(bins=200) +
    ggtitle(basename)


dev.off()
