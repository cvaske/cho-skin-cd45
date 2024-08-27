#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(ggplot2)

argv <- commandArgs(trailingOnly=TRUE)

stopifnot(
    length(argv) == 5 &&
    is.character("args: seuratobj label clusters.tab cluster.id outdir")
)

seuratobj <- argv[1]
label <- argv[2]
clustersfn <- argv[3]
cluster.id <- argv[4]
outdir <- argv[5]

## constants
min.cells <- 1
platclasses <- c("Skin","PBMC")
pbmc.pattern <- "PBMC"

clusters <- read.table(
    clustersfn,
    row.names=1,
    header=TRUE, 
    sep=","
)
head(clusters)

load(seuratobj)
seurat <- AddMetaData(object = seurat, metadata=clusters)

target.cells <- rownames(filter(clusters, predicted.ID == cluster.id))

ggsave(
    file.path(outdir, paste(label, "-umap-all-dims.png", sep="")),
    DimPlot(seurat)
    + labs(title=paste(label, "on all Harmony dimensions"))
)

pbmc <- as.integer(grepl(x=seurat@meta.data$orig.ident, pattern=pbmc.pattern))
pbmc <- factor(platclasses[1 + pbmc])

h.df <- data.frame(Source=pbmc, seurat@reductions$harmony@cell.embeddings)
ttests <- lapply(1:20, function(i) {
    t.test(as.formula(paste("harmony_", i, " ~ Source", sep="")), data=h.df)
})
p.values <- sapply(ttests, "[[", "p.value")

pdf(file.path(outdir, paste(label, "-dimension-pvalues.pdf", sep="")))
plot(p.values)
dev.off()

harmony.dims <- which(p.values > median(p.values))
seurat <- RunUMAP(seurat, dims = harmony.dims, verbose = FALSE)

DimPlot(seurat, cells.highlight = target.cells) + labs(title=paste(label, "on most platform-neutral dimensions"))

cat("Writing new UMAP plot with selected cells\n")
ggsave(
    file.path(outdir, paste(label, "-umap-plat-neutral-selected.png", sep=""))
)
cat("Writing new UMAP plot\n")
ggsave(
    file.path(outdir, paste(label, "-umap-plat-neutral.png", sep="")),
    DimPlot(seurat)
    + labs(title=paste(label, "on most platform-neutral dimensions"))
)

cat("Reclustering the data\n")
seurat <- seurat |> 
    FindNeighbors(dims = harmony.dims, reduction="harmony") |> 
    FindClusters(resolution = 0.5)
ggsave(
    file.path(outdir, paste(label, "-umap-plat-neutral-reclustered.png", sep="")),
    DimPlot(seurat, group.by="seurat_clusters")
    + labs(title=paste(label, "on most platform-neutral dimensions"))
)
cat("Writing confusion matrix of old vs. new clusters\n")
write.table(
    table(seurat@meta.data$predicted.ID, seurat@meta.data$seurat_clusters),
    file=file.path(outdir, paste(label, "-predicted-vs-newclst.tsv", sep="")),
    sep="\t"
)
cat("Writing new clustering vs. sample source\n")
write.table(
    table(seurat@meta.data$seurat_clusters, seurat@meta.data$orig.ident),
    file=file.path(outdir, paste(label, "-newclust-vs-plat.tsv", sep="")),
    sep="\t"
)

overlapping_clusters <- names(which(table(seurat@meta.data$predicted.ID, seurat@meta.data$seurat_clusters)[cluster.id,] >= min.cells))

best_cluster <- names(which.max(table(seurat@meta.data$predicted.ID, seurat@meta.data$seurat_clusters)[cluster.id,]))

seurat@meta.data$best_target_cluster <- best_cluster == seurat@meta.data$seurat_clusters

seurat@meta.data$overlap_target_cluster <- seurat@meta.data$seurat_clusters %in% overlapping_clusters

write.table(
    seurat@meta.data,
    file.path(outdir, paste(label, "-metadata.tsv", sep="")),
    sep="\t"
)
