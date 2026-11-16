
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
my_colors <- c(brewer.pal(12, "Paired"), brewer.pal(5, "Set3"))
set.seed(777)

so.combined <- ScaleData(so.combined)
so.combined <- RunPCA(so.combined)
so.combined <- RunUMAP(so.combined,dims = 1:30)
so.combined <- FindNeighbors(so.combined,reduction = "pca",dims = 1:30)
so.combined <- FindClusters(so.combined,resolution = 0.1)

ncol(so.combined) # number of cells
table(so.combined$seurat_clusters) # number of each cluster


so.combined$genotype_stage <- paste(so.combined$genotype, so.combined$stage, sep ="_")

so.combined_mt5 <- so.combined
DefaultAssay(so.combined_mt5) <-"RNA"
so.combined_mt5[["percent.mt"]] <- PercentageFeatureSet(so.combined_mt5, pattern = "^mt-")
so.combined_mt5 $ orig_Cluster <- so.combined_mt5 $seurat_clusters
so.combined_mt5 <- subset(so.combined_mt5, nFeature_RNA > 200 & percent.mt <5 )
DefaultAssay(so.combined_mt5) <- "integrated"
so.combined_mt5 <- ScaleData(so.combined_mt5)
so.combined_mt5 <- RunPCA(so.combined_mt5)
so.combined_mt5 <- RunUMAP(so.combined_mt5,dims = 1:30)
so.combined_mt5 <- FindNeighbors(so.combined_mt5,reduction = "pca",dims = 1:30)
so.combined_mt5 <- FindClusters(so.combined_mt5,resolution = 0.1)

ncol(so.combined_mt5)
table(so.combined_mt5 $ seurat_clusters)

#sup7.a
DimPlot(so.combined_mt5, raster = F, label = T, cols = my_colors, label.size = 8)
#sup7.b.c.d.e
DimPlot(so.combined_mt5,split.by = "genotype_stage",raster = FALSE,label=TRUE,pt.size = 0.01,ncol = 2, cols = my_colors)&NoAxes()
#sup7.f
DimPlot(so.combined_mt5,reduction = "umap",group.by = "seurat_clusters",split.by = "library",ncol = 4,pt.size = 0.01, cols = my_colors, label.size = 8)
#sup7.g
DimPlot(so.combined_mt5,reduction = "umap",group.by = "library",shuffle = TRUE,pt.size = 0.01,cols = my_colors, label.size = 8)


DefaultAssay(so.combined_mt5) <- "SCT"
so.combined_mt5_1 <- so.combined_mt5
new.cluster.ids1 <- c("GC", "GCP", "GC", "GCP","SOX9-positive cell", "Undetermined","Inhibitory neuron", "Oligodendrocyte lineage cell","Choroid plexus cell", "Microglia","Endothelial cell","Meningeal cell","Blood cell","Unipolar Brush cell","Pericyte")
names(new.cluster.ids1) <- levels(so.combined_mt5_1)
so.combined_mt5_2 <- RenameIdents(so.combined_mt5_1, new.cluster.ids1)
#Sup7.h
VlnPlot(so.combined_mt5_2,features = c("Pcna","Ccnd1","Hey1","Sfrp1","Reln","Cntn2","Neurod1","Aldoc","Aqp4","Fabp7","Pax2","Gad1","Olig1","Olig2","Ttr","Aif1","Cx3cr1","Pecam1","Vwf","Cxcl12","S100a6","Col1a1","Hba-a1","Hba-a2","Eomes","Rgs5","Vtn","Rsph1","Foxj1","Slc17a6","Tbr1","Calb1","Foxp2"),stack = T, split.by = "seurat_clusters",flip = T, cols = my_colors)+
  coord_cartesian(clip = "off")
#Sup7.i
VlnPlot(so.combined_mt5_2,features = c("Aldoc","Aqp4","Sox2","Sox9","S100b","Fgfr1","Fabp7","Hopx","Meig1","Fam183b","Rsph1"),stack = T, split.by = "seurat_clusters",flip = T, cols = my_colors)+
  coord_cartesian(clip = "off")

saveRDS(so.combined_mt5_2, "~/Desktop/Ichijo_filtermt5/so.combined_mt5_annotated.rds")

ASlin <- subset(so.combined_mt5, idents = 4)
DefaultAssay(ASlin) <- "integrated"
ASlin0.5 <- ASlin %>% RunPCA() %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5)
DimPlot(ASlin0.5)
DefaultAssay(ASlin0.5) <- "SCT"
VlnPlot(ASlin0.5, features = c("Sox2","Sox9","Slc6a11","Gria1","Gdf10","Sept4","Pcna","Mki67"),stack = T, split.by = "seurat_clusters",flip = T, cols = my_colors)+
  coord_cartesian(clip = "off")

ASlin0.5_sc5 <- subset(ASlin0.5, idents = c(0,1,2,6,7)) # remove sc3,4,5 for low slc6a11 Glia1, Gdf10, Sept4 expression level

saveRDS(ASlin0.5_sc5, "~/Desktop/Ichijo_filtermt5/ASlin0.5_sc5.rds")

DefaultAssay(ASlin0.5_sc5) <- "integrated"
ASlin0.5_sc5_0.1 <- ASlin0.5_sc5 %>% RunPCA() %>% RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.1)
#Fig2.a
DimPlot(ASlin0.5_sc5_0.1, raster = F, cols = my_colors,label=T,label.size = 8)
#Fig2.b
DimPlot(ASlin0.5_sc5_0.1, raster = F, cols = my_colors, split.by = "genotype",label = T, label.size = 8)
DefaultAssay(ASlin0.5_sc5_0.1) <- "SCT"
#Fig2.c
VlnPlot(ASlin0.5_sc5_0.1, features = c("Sox2","Sox9","Slc6a11","Gria1","Gdf10","Sept4","Mki67"),stack = T, split.by = "seurat_clusters",flip = T, cols = my_colors)+
  coord_cartesian(clip = "off")

Ast <- subset(ASlin0.5_sc5_0.1, idents = 0)
Idents(Ast) <- "Ast"
BG_like_cell <- subset(ASlin0.5_sc5_0.1, idents = c(1,2))
Idents(BG_like_cell) <- "BGlike"
merge_AstBG <- merge(BG_like_cell,y=Ast)
merge_AstBG$genotype <- factor(merge_AstBG$genotype,levels = c("ctr", "cKO"))

DefaultAssay(merge_AstBG) <- "SCT"
#Fig6.b
VlnPlot(merge_AstBG, features = "Vim",split.by = "genotype") 
VlnPlot(merge_AstBG, features = "Zeb2",split.by = "genotype") 
VlnPlot(merge_AstBG, features = "Tnc",split.by = "genotype")
VlnPlot(merge_AstBG, features = "Slc1a3",split.by = "genotype")

#Sup10.a
VlnPlot(merge_AstBG, features = "Gria1",split.by = "genotype")
VlnPlot(merge_AstBG, features = "Adgrg6", split.by = "genotype")
VlnPlot(merge_AstBG, features = "Sept4", split.by = "genotype")


DefaultAssay(BG_like_cell) <-"SCT"
Idents(BG_like_cell) <- "genotype"
obj <- PrepSCTFindMarkers(BG_like_cell)
markers <- FindMarkers(obj, ident.1="cKO", ident.2 ="ctr",only.pos = FALSE)
write.csv(markers, "~/Desktop/Ichijo_revise_1020/BG_like_cell_ctrvscKO_mt5_findmarkers.csv")



DefaultAssay(Ast) <- "SCT"
Idents(Ast) <- "genotype"
obj0 <- PrepSCTFindMarkers(Ast)
markers0 <- FindMarkers(obj0, ident.1="cKO", ident.2 ="ctr",only.pos = FALSE)
write.csv(markers0, "~/Desktop/Ichijo_revise_1020/Ast_ctrvscKO_mt5_findmarkers.csv")

BG <- subset(ASlin0.5_sc5_0.1, idents = 1)
DefaultAssay(BG) <-"SCT"
Idents(BG) <- "genotype"
obj1 <- PrepSCTFindMarkers(BG)
markers1 <- FindMarkers(obj1, ident.1="cKO", ident.2 ="ctr",only.pos = FALSE)
write.csv(markers1, "~/Desktop/BG_ctrvscKO_mt5_findmarkers.csv")

BGLP <- subset(ASlin0.5_sc5_0.1, idents = 2)
DefaultAssay(BGLP) <-"SCT"
Idents(BGLP) <- "genotype"
obj2 <- PrepSCTFindMarkers(BGLP)
markers2 <- FindMarkers(obj2, ident.1="cKO", ident.2 ="ctr",only.pos = FALSE)
write.csv(markers2, "~/Desktop/BGLP_ctrvscKO_mt5_findmarkers.csv")








