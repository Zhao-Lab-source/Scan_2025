library(ggplot2)
library(Seurat)
library(readr)
library(tidyverse)
library(DoubletFinder)
library(plyr)
library(remotes)
library(devtools)
library(future) 
library(cowplot)
library(RColorBrewer)

setwd("./")
options(future.globals.maxSize = 100*1000 * 1024^2)
plan(multisession, workers = 4)

dir = c("Mut10/", "Mut8/","WT/")
names(dir) = c('MUT10','MUT8','WT')

scRNAlist <- list()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts,min.cells = 3, min.features = 300)
}

#qc
for (i in seq_along(scRNAlist)) {
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "(?i)^MT-")  
  scRNAlist[[i]][["percent.rp"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "(?i)^RP[SL]")  
  
  scRNAlist[[i]] <- subset(scRNAlist[[i]], 
                           subset = nFeature_RNA > 200 & 
                             nFeature_RNA < 7000 & 
                             percent.mt < 10 & percent.rp < 40)
}

for (i in seq_along(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]]) %>% 
    FindVariableFeatures(., nfeatures = 2000)%>%     
    ScaleData(., vars.to.regress = c("percent.mt", "percent.rp"), verbose = FALSE)
}

#run PCA
for (i in seq_along(scRNAlist)) {
  scRNAlist[[i]] <- RunPCA(scRNAlist[[i]], npcs = 50, verbose = FALSE)}

#cluster
for (i in seq_along(scRNAlist)) {
  scRNAlist[[i]] <- RunUMAP(scRNAlist[[i]], reduction = 'pca', dims = 1:10)%>%
    FindNeighbors(., reduction = 'pca', dims = 1:10)%>%
    FindClusters(., resolution = 0.5)
}

save(scRNAlist,file = "scRNAlist.Rdata")

load("scRNAlist.Rdata")

#find doublet
detectDoublet <- function(obj, dims, expected_doublet_rate, pK_value, ncores = 1, SCTransform = FALSE, Homotypic = FALSE, annotation) {
  
  if (SCTransform) {
    obj <- SCTransform(obj, verbose = FALSE)
  }
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
  
  homotypic <- Homotypic
  if (homotypic) {
    message("Consider isotype duplex cells for testing.")
  }
 
  nExp <- round(expected_doublet_rate * ncol(obj))
  
  obj <- doubletFinder(
    obj, 
    pN = 0.25, 
    pK = pK_value,           
    nExp = nExp,             
    PCs = 1:30  
  )
  
  obj$DF.classify <- obj@meta.data$DF.classifications
  
  return(obj)
}

scRNAlist[[1]] <- detectDoublet(
  scRNAlist[[1]], 
  dims = 1:20, 
  expected_doublet_rate = 0.0255, 
  pK_value = 0.09,                                
  ncores = 2, 
  SCTransform = FALSE, 
  Homotypic = FALSE, 
  annotation = "seurat_clusters"
)

scRNAlist[[2]] <- detectDoublet(
  scRNAlist[[2]], 
  dims = 1:20, 
  expected_doublet_rate = 0.0276, 
  pK_value = 0.09,
  ncores = 2, 
  SCTransform = FALSE, 
  Homotypic = FALSE, 
  annotation = "seurat_clusters"
)

scRNAlist[[3]] <- detectDoublet(
  scRNAlist[[3]], 
  dims = 1:20, 
  expected_doublet_rate = 0.0205, 
  pK_value = 0.09,
  ncores = 2, 
  SCTransform = FALSE, 
  Homotypic = FALSE, 
  annotation = "seurat_clusters"
)

duPlot <- list()
for (i in seq_along(scRNAlist)) {
  p = DimPlot(scRNAlist[[i]], group.by = "DF.classify")+
    ggtitle(unique(scRNAlist[[i]]$orig.ident))
  duPlot[[i]] <- p
}

plot_grid(duPlot[[1]],duPlot[[2]],duPlot[[3]], ncol = 3)

###remove Doublet
for (i in seq_along(scRNAlist)) {
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = (DF.classify == "Singlet"))
}


scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = 2000)
scRNA1 <- IntegrateData(anchorset = scRNA.anchors)

DefaultAssay(scRNA1) <- "RNA"
scRNA1 <- ScaleData(scRNA1, features = rownames(scRNA1))
DefaultAssay(scRNA1) <- "integrated"
scRNA1 <- ScaleData(scRNA1)

scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 2000)
scRNA1 <- RunPCA(scRNA1)
DimPlot(scRNA1, reduction = "pca",group.by ="orig.ident")
ElbowPlot(scRNA1,ndims=50,reduction="pca")
scRNA1 <- FindNeighbors(scRNA1,dims= 1:10)
scRNA1 <- FindClusters(scRNA1, resolution = 0.2)
scRNA1 <- RunUMAP(scRNA1,dims= 1:10)

colors <- c("#AEDDFF","#E3B0A3","#FCE8B4")
colors1 <- brewer.pal(n = 7, name = "Set3")

DimPlot(scRNA1,reduction = "umap",group.by = "orig.ident",cols = colors) +theme(aspect.ratio = 1)
DotPlot(scRNA1,features = c("CD3D" ,"CD3G","CD8A","CD8B"),group.by = "integrated_snn_res.0.2")

save(scRNA1,file = "scRNAlist2.Rdata")

##
#load("scRNAlist2.Rdata")
sce_T <- subset(scRNA1,
                subset = (seurat_clusters %in% c(6)))
sce_T

DefaultAssay(sce_T) <- "RNA"
sce_T <- NormalizeData(sce_T, normalization.method = "LogNormalize", scale.factor = 1e4)
sce_T <- FindVariableFeatures(sce_T, selection.method = "vst", nfeatures = 2000)
sce_T <- ScaleData(sce_T)
sce_T <- RunPCA(sce_T, features = VariableFeatures(object = sce_T))
ElbowPlot(sce_T,ndims=50,reduction="pca")
sce_T <- FindNeighbors(sce_T, reduction = "pca", dims = 1:10)
sce_T  <- FindClusters(sce_T, resolution = 0.6)
sce_T  <- RunUMAP(sce_T ,dims= 1:10)
sce_T  <- RunTSNE(sce_T ,dims= 1:10)

markers<-FindAllMarkers(sce_T,only.pos=TRUE,min.pct=0.25,logfc.threshold=0.25)

DimPlot(sce_T ,reduction = "umap",group.by = "orig.ident",cols = colors) +theme(aspect.ratio = 1)
DimPlot(sce_T ,reduction = "umap",group.by = "RNA_snn_res.0.6",cols = colors1) +theme(aspect.ratio = 1)
save(sce_T,file = "sce_T.Rdata")

###Extended Data
#Teff signature

DefaultAssay(sce_T) <- "RNA"
sce_T <- JoinLayers(sce_T, assay = "RNA", mode = "data")
sce_T@meta.data$orig.ident <- factor(sce_T@meta.data$orig.ident, levels = c("WT", "MUT8", "MUT10"))

Teff.sig <- list(c("IL2RA","TNFRSF8","TNFRSF4","CD69","CCL3","CCL4","CCL5","XCL1"))
sce_T <- AddModuleScore(object = sce_T, features = Teff.sig, name = "Teff.sig")
FeaturePlot(object = sce_T, features = "Teff.sig1")
FeaturePlot(object = sce_T, features = "Teff.sig1", split.by = "orig.ident") + 
  theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 7, "RdBu")),
                                                  limits = c(-0.5, 1.3)) & theme(legend.position = c(0.9, 0.9))
colors <- c("#AEDDFF","#E3B0A3","#FCE8B4")
VlnPlot(sce_T, features = c("Teff.sig1"), group.by = "orig.ident", ncol=1, pt.size=0, cols = colors) + theme(aspect.ratio = 1)
pdf("Teff_VlnPlot.pdf", width = 17, height = 6)
VlnPlot(sce_T, features = c("IL2RA","TNFRSF8","TNFRSF4","CD69","CCL3","CCL4","CCL5","XCL1"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()
dev.off()


Teff1.sig <- list(c("IFNG","TNF","PRF1","GZMA","GZMB","LAMP1"))
sce_T <- AddModuleScore(object = sce_T, features = Teff1.sig, name = "Teff1.sig")
FeaturePlot(object = sce_T, features = "Teff1.sig1") #-0.5-1.0
FeaturePlot(object = sce_T, features = "Teff1.sig1", split.by = "orig.ident") + 
  theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 5, "RdBu")),
                                                  limits = c(-0.5, 1.0)) & theme(legend.position = c(0.9, 0.9))

VlnPlot(sce_T, features = c("Teff1.sig1"), group.by = "orig.ident", ncol=1, pt.size=0, cols = colors) + theme(aspect.ratio = 1)
pdf("Teff1_VlnPlot.pdf", width = 17, height = 6)
VlnPlot(sce_T, features = c("IFNG","TNF","PRF1","GZMA","GZMB","LAMP1"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()
dev.off()


#TRM signature
TRM.sig <- list(c( "TCF7", "LEF1", "BACH2", "FOXO1", "ID3", "IL7R", "BCL2", "MCL1", "SELL", "CCR7", "CD27", "CD28", "KLF2", "SATB1", "RBPJ", "NOTCH1"))
sce_T <- AddModuleScore(object = sce_T, features = TRM.sig, name = "TRM.sig")
FeaturePlot(object = sce_T, features = "TRM.sig1")
FeaturePlot(object = sce_T, features = "TRM.sig1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 2, "PRGn")), limits = c(-0.2, 0.4)) & theme(legend.position = c(0.9, 0.9))
VlnPlot(sce_T, features = c("TRM.sig1"), group.by = "orig.ident", ncol=1, pt.size=0, cols = colors) + theme(aspect.ratio = 1)

pdf("TRM_VlnPlot.pdf", width = 17, height = 6)
VlnPlot(sce_T, features = c( "TCF7", "LEF1", "BACH2", "FOXO1", "ID3", "IL7R", "BCL2", "MCL1", 
                             "SELL", "CCR7", "CD27", "CD28", "KLF2", "SATB1", "RBPJ", "NOTCH1"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()
dev.off()


#TEX signature
Tex.sig <- list(c("TOX", "NR4A3", "CD244", "HAVCR2",
                  "LAG3", "TIGIT", "PDCD1","CTLA4"))
sce_T <- AddModuleScore(object = sce_T, features = Tex.sig, name = "Tex.sig")
FeaturePlot(object = sce_T, features = "Tex.sig1")
FeaturePlot(object = sce_T, features = "Tex.sig1", split.by = "orig.ident") + theme(aspect.ratio = 1) & scale_color_gradientn(colours = rev(brewer.pal(n = 4, "PRGn")), limits = c(-0.4, 0.4)) & theme(legend.position = c(0.9, 0.9))
VlnPlot(sce_T, features = c("Tex.sig1"), group.by = "orig.ident", ncol=1, pt.size=0, cols = colors) + theme(aspect.ratio = 1)
pdf("TEX_VlnPlot.pdf", width = 17, height = 6)
VlnPlot(sce_T, features = c("TOX", "NR4A3", "CD244", "HAVCR2","LAG3", "TIGIT", "PDCD1","CTLA4"), group.by = "orig.ident", stack = TRUE) + theme(aspect.ratio = 1) + NoLegend() + coord_flip()
dev.off()
