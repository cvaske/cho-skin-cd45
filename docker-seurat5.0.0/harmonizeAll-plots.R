#!/usr/bin/env -S Rscript --verbose

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)

argv <- commandArgs(trailingOnly=TRUE)

samples.md.fn <- argv[1]
label <- argv[2]
outdir <- argv[3]
stopifnot(length(argv) == 3)

umap.dims <- 1:20

outpdf <- paste(outdir, "/", label, ".pdf", sep="")
outtab <- paste(outdir, "/", label, ".tsv", sep="")
outRobj <- paste(outdir, "/", label, ".RData", sep="")

samples.md <- read.table(samples.md.fn, header=TRUE)
print(samples.md)
stopifnot(nrow(samples.md) > 1)

load.subset <- function(h5, project, nFeature.min, nFeature.max, percent.mt.max, sample) {
    so <- CreateSeuratObject(counts = Read10X_h5(h5), project = project)
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    orig.cells <- ncol(so$RNA)
    so <- subset(so, subset = nFeature_RNA < nFeature.max & nFeature_RNA > nFeature.min & percent.mt < percent.mt.max)
    return(list(seurat=so, orig.cells = orig.cells, kept.cells = ncol(so$RNA), percent.kept = ncol(so$RNA) / orig.cells))
}

load.mergeset <- function(md) {
    subsets <- mapply(FUN=load.subset,
                      sample = md$Sample,
                      h5 = md$H5, 
                      project = md$Sample, 
                      nFeature.min = md$nFeature.min, 
                      nFeature.max = md$nFeature.max,
                      percent.mt = md$percent.mt.max,
                      SIMPLIFY=FALSE)
    seurats <- sapply(subsets, '[[', 'seurat')
    subset.df <- data.frame(
        Sample = md$Sample, 
        percent.kept = sapply(subsets, '[[', "percent.kept"),
        orig.cells = sapply(subsets, '[[', "orig.cells"),
        kept.cells = sapply(subsets, '[[', "kept.cells")
    )
    
    genes.intersection <- rownames(seurats[[1]])
    for (i in 2:length(seurats)) {
        genes.intersection <- intersect(genes.intersection, rownames(seurats[[i]]))
    }

    seurat <- merge(seurats[[1]], y = seurats[-1], add.cell.ids = md$Sample, project = label)
    seurat <- seurat[genes.intersection,]

    print("Joining layers...")
    # Very important to subset to [["RNA"]] in the right side here, or it consumes all RAM
    # https://github.com/satijalab/seurat/issues/8208
    seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])

    print("Normalizing Data...")
    seurat <- NormalizeData(seurat)

    VariableFeatures(seurat) <- split(row.names(seurat@meta.data), seurat@meta.data$orig.ident) %>%
        lapply(function(cells_use) {
        seurat[,cells_use] %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
            VariableFeatures()
    }) %>% unlist %>% unique

    print("Scaling...")
    seurat <- seurat %>% ScaleData(verbose = FALSE)
    print("Running PCA...")
    seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 20, verbose = FALSE)
    
    return(list(
        seurat=seurat, 
        df=subset.df
    ))
}

mergeset <- load.mergeset(samples.md)
print("Loaded mergeset")

write.table(mergeset$df, file = outtab, row.names=FALSE, sep="\t")

seurat <- mergeset$seurat

pdf(outpdf)

seurat <- seurat %>% RunUMAP(reduction = "pca", dims=umap.dims)

DimPlot(seurat, reduction = "umap") + labs(title=paste(label, "UMAP on PCA"))
DimPlot(seurat, reduction = "pca") + labs(title=paste(label, "PCA"))
DimPlot(seurat, reduction = "pca", dims=3:4)  + labs(title=paste(label, "PCA"))

seurat <- seurat %>% RunHarmony(group.by.vars = "orig.ident", plot_convergence = TRUE)

DimPlot(seurat, reduction = "harmony") + labs(title=paste(label, "Harmony"))
DimPlot(seurat, reduction = "harmony", dims=3:4)  + labs(title=paste(label, "Harmony"))

seurat <- seurat %>% RunUMAP(reduction = "harmony", dims=umap.dims)
DimPlot(seurat, reduction = "umap") + labs(title=paste(label, "UMAP on Harmony"))
dev.off()

save(seurat, file = outRobj)
