library(Seurat)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(future)
library(harmony)

plan("multiprocess", workers = 30)
options(future.globals.maxSize = 100*(1024^3))

setwd('./')
da.smp <- read.table('/NAS/lyc/02.analysis/NSCLC_BM/d202108/01.singlecell/integrated.harmony.20211214/sample_path.txt', header = T, sep = '\t', stringsAsFactors = F)

my_sce = lapply(1:nrow(da.smp), function(n) {
  sce.data <- Read10X(data.dir = da.smp$data_path[n])
  sce <- CreateSeuratObject(counts = sce.data, project = as.character(da.smp$Sample[n]), min.cells = 3, min.features = 200)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[LS]")
  HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes_total,rownames(sce@assays$RNA))
  HB.genes <- rownames(sce@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  sce[["percent.HB"]] <- PercentageFeatureSet(sce,features=HB.genes)
  f1 <- VlnPlot(sce, features = c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.HB", "percent.ribo"), ncol = 5)
  ggsave(f1,filename = paste0(da.smp$Sample[n],'.qc_raw.png'),device='png', width = 12)
  sfe <- subset(sce, subset = nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 20)
  return(sfe)
})

my_sce = lapply(my_sce, function(se) se = se %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000))

scimm <- merge(my_sce[[1]], y = my_sce[2:length(my_sce)],  project = "NSCLC_immunotherapy", merge.data = TRUE)
VariableFeatures(scimm) <- SelectIntegrationFeatures(object.list = my_sce, nfeatures = 2000)
scimm <- ScaleData(scimm, verbose = FALSE)


#integrated
scimm <- RunPCA(object = scimm, npcs = 50)
scimm <- RunHarmony(object = scimm,  reduction = "pca", dims.use = 1:40, group.by.vars = "orig.ident", plot_convergence = TRUE)
scimm <- RunUMAP(object = scimm, reduction = "harmony", dims = 1:40)
scimm <- FindNeighbors(object = scimm, reduction = "harmony", dims = 1:40)

for(i in seq(0.1, 1, 0.1)){
  scimm <- FindClusters(object = scimm, verbose = T, resolution = i)
  f1 <- DimPlot(scimm, reduction = "umap",label = TRUE)
  ggsave(f1,filename=paste0('PC_40.res_',i,'.umap.pdf'),device='pdf')
}

scimm <- FindClusters(object = scimm, resolution = 0.8)
saveRDS(scimm, file='integrated.harmony.rds')

f2 <- DimPlot(scimm, reduction = "umap",label = TRUE)
ggsave(f2,filename='umap.pdf',device='pdf')
f3 <- DimPlot(scimm, reduction = "umap",label = F,group.by='orig.ident')
ggsave(f3,filename = 'sample.pdf',device = 'pdf')


##cluster marker plot
scimm.markers <- FindAllMarkers(scimm, only.pos = FALSE, min.pct = 0.3, logfc.threshold = 0.5, min.diff.pct =0.15)
write.table(scimm.markers,file='Seurat.diffexp.feature.txt',sep = '\t',quote=F,row.names = F)

scimm.markers %>% group_by(cluster) %>% slice_max(n=10, avg_log2FC) -> top_10_marker
DefaultAssay(scimm) <- 'RNA'
lapply(unique(top_10_marker$cluster), function(clu) {
  top_10_marker %>% filter(cluster==m) -> top_10_marker_part
  lapply(top_10_marker_part$gene, function(ge) {
    f5<-FeaturePlot(scimm, features = ge,min.cutoff=0,max.cutoff=20)
    ggsave(f5,filename = paste0('featureplot.cluster_',clu,'.',ge,'.png'),device = 'png')
    f6<-VlnPlot(scimm, features = ge) + NoLegend()
    ggsave(f6,filename = paste0('vlnplot.cluster_',clu,'.',ge,'.png'),device = 'png')
  })
})

