# "T3 treated organoids"
# '2022-08-30'
#Load required packages
library(Seurat)
library(tidyverse)
library(BiocParallel)
library(SeuratWrappers)
library(scDblFinder)
library(SingleCellExperiment)
set.seed(1234)

#Load 10x data
files <- list.dirs("D:/R_files/Human organoid/Raw/", recursive = FALSE)
for (i in files){
  name <- basename(i)
  data <- Read10X(data.dir = i)
  obj <- CreateSeuratObject(counts = data, min.cells = 3,
                            min.features = 200, project = name)
  assign(name, obj)
  rm (obj, data, name,i)
}

#Add metadata
Con1@meta.data$treatment <- "Control" 
Con2@meta.data$treatment <- "Control" 
T3_10day1@meta.data$treatment <- "10day" 
T3_10day2@meta.data$treatment <- "10day" 
T3_200day1@meta.data$treatment <- "200day" 
T3_200day2@meta.data$treatment <- "200day" 
T3_200day3@meta.data$treatment <- "200day" 

#Merge datasets 
hgorg <- merge(Con1, y = c(Con2,T3_10day1,T3_10day2,T3_200day1,T3_200day2,T3_200day3), add.cell.ids = c("Con1","Con2","T3_10day1","T3_10day2","T3_200day1","T3_200day2","T3_200day3"), project = "T3_treated")

#Add percent mitocondrial gene expression
hgorg[["percent.mt"]] <- PercentageFeatureSet(hgorg, pattern = "^MT-")
hgorg[["percent.ribo"]] <- PercentageFeatureSet(hgorg, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

#Order meta data
order1<- c("Con1","Con2","T3_10day1","T3_10day2","T3_200day1","T3_200day2","T3_200day3")
hgorg@meta.data$orig.ident <- factor(hgorg@meta.data$orig.ident, levels = order1)
order2<- c("Control","10day","200day")
hgorg@meta.data$treatment <- factor(hgorg@meta.data$treatment, levels = order2)

#Add cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
hgorg <- CellCycleScoring(hgorg, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#Add number of genes per UMI for each cell to metadata
hgorg@meta.data$log10GenesPerUMI <- log10(hgorg@meta.data$nFeature_RNA)/log10(hgorg@meta.data$nCount_RNA)

#Filter out low quality cells
hgorg <- subset(x = hgorg, subset= percent.mt < 25)
hgorg <- subset(x = hgorg, subset= nCount_RNA > 1000)

#Identify doublets
sce <- as.SingleCellExperiment(hgorg)
sce <- scDblFinder(sce, samples="orig.ident",BPPARAM=SerialParam(RNGseed = 1234)) #compatible with windows machines
hgorg@meta.data$doublet <- plyr::mapvalues(x = rownames(hgorg@meta.data), 
                                           from = rownames(colData(sce)), to = as.character(sce$scDblFinder.class))
hgorg@meta.data$doublet <- factor(hgorg@meta.data$doublet, levels = c("singlet","doublet"))
hgorg@meta.data$doublet_score <- plyr::mapvalues(x = rownames(hgorg@meta.data), 
                                                 from = rownames(colData(sce)), to = sce$scDblFinder.score)
hgorg@meta.data$doublet_score<-as.numeric(hgorg@meta.data$doublet_score)
rm (sce)

#Filter out doublets and cells with high UMI
hgorg <- subset(x = hgorg, subset= nCount_RNA < 20000)
hgorg <- subset(x = hgorg, subset= doublet == "singlet")

#Prepare the data for UMAP dimension reduction and find clusters
hgsct <- SCTransform(hgorg, vars.to.regress = c("percent.mt", "percent.ribo")) %>% RunPCA(verbose = FALSE)
hgsct <- RunUMAP(hgsct, dims = 1:14)
hgsct <- FindNeighbors(hgsct, dims = 1:14)
hgsct <- FindClusters(hgsct, resolution = 1.4)

#Use feature genes from Hoang et al., Science (2020), Clark et al., Neuron (2019) and Lu et al., Dev Cell (2020)
#Classify celltypes 
Idents(hgsct) <- "seurat_clusters"
old<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18",
       "19","20","21","22","23","24","25","26","27","28","29","30")
new<-c("Bipolar","Rod","Rod","Glia-like","Rod","Cone","RPCs","RPCs","RPCs","Cone","Cone","Cone",
       "Muller glia","RPCs","Rod","Amacrine","Rod","RPCs","PR Precursors","Rod", "Cone","Rod",
       "Neurogenic PCs","PR Precursors","Glia-like","Amacrine","Cone","Bipolar","RPCs","Proliferating RPCs","Doublet")
hgsct@meta.data$Celltype <- plyr::mapvalues(x = hgsct@meta.data$seurat_clusters, from = old, to = new)
Idents (hgsct) <- "Celltype"

# Use cell selector to label Horizontal cells and doublets not identified by scDblFinder
hz <- WhichCells(hgsct, expression = LHX1 >= 1)
hz <- CellSelector(FeaturePlot(hgsct, features = "LHX1", cells = hz)) #Select cells in amacrine population
doubl <- CellSelector(DimPlot(hgsct)) #Select cells in near cluster 33 belonging to cluster 8
Idents(hgsct, cells = hz) <- "Horizontal"
Idents(hgsct, cells = doubl) <- "Doublet"
hgsct[["Celltype"]] <- Idents(hgsct)

#Remove remaining doublets and order celltypes
hgsct <- subset(hgsct, idents= "Doublet", invert = T)
order<- c("PR Precursors","Rod","Cone","Bipolar","Amacrine","Horizontal",
          "Muller glia","Glia-like","RPCs", "Proliferating RPCs", "Neurogenic PCs")
hgsct@meta.data$Celltype <- factor(hgsct@meta.data$Celltype, levels = order)
Idents(hgsct) <- "Celltype"

#Remove genes that have rowSum = 0 from count matrix and rerun Normalization and dimension reduction
meta <- hgsct@meta.data
counts <- hgsct@assays$RNA@counts [rowSums(hgsct@assays$RNA@counts)!=0,]
hgsct_clean <- CreateSeuratObject(counts = counts, meta.data = meta)
hgsct_clean <- SCTransform(hgsct_clean, vars.to.regress = c("percent.mt", "percent.ribo")) %>% RunPCA(verbose = FALSE)
hgsct_clean <- RunUMAP(hgsct_clean, dims = 1:10)
Idents(hgsct_clean) <- "Celltype"
#For 3D UMAP
hgsct_clean <- RunUMAP(hgsct_clean, dims = 1:10, n.components = 3L, reduction.name = "UMAP3D",reduction.key = "UMAP3D")

#Subset photoreceptors objects
photo <- subset(hgsct_clean, idents= c("PR Precursors","Rod","Cone"))
pr_cell <- CellSelector(DimPlot(photo)) # Use CellSelector to select cells in PR cluster
photo <- subset(photo, cells = pr_cell)
cone <- subset(photo, idents= c("PR Precursors","Cone"))

#Label cells used for photoreceptor objects
hgsct_clean@meta.data<-hgsct_clean@meta.data %>% 
  mutate(used_for_photo = ifelse(rownames(hgsct_clean@meta.data) %in% rownames(photo@meta.data), "Yes", "No"))
hgsct_clean@meta.data<-hgsct_clean@meta.data %>% 
  mutate(used_for_cones = ifelse(rownames(hgsct_clean@meta.data) %in% rownames(cone@meta.data), "Yes", "No"))


