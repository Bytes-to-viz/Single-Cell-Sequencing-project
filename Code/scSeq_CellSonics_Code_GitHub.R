#R version 4.1.2 (2021-11-01) -- "Bird Hippie"
#Copyright (C) 2021 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)



#TITLE: SCRNASEQ CELLSONICS



# PART 0: PREPARE THE ENVIRONMENT --------------------------------

#1. Load the necessary packages. 
library("Biobase")
library("GEOquery")
library("limma")
library("readr")
library("Seurat")
library("Rtsne") 
library("stats") 
library("ggplot2")
library("GSVA")
library("GSEABase")
library("ggpubr")
library("ggExtra")
library("ComplexHeatmap")
library("dplyr")
library("circlize")
library("ggrepel")
library("networkD3")
library("ggalluvial")
library("DropletUtils")
library("gridExtra")
library("stringr")
library("viridis")
library("escape")
library("dplyr")
library("data.table")
library("Seurat")
library("SingleR")
library("dplyr")
library("celldex")
library("RColorBrewer")
library("SingleCellExperiment")
library("writexl")

#Links to guidelines used: 
#General: https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html
#For DGE between groups: https://satijalab.org/seurat/archive/v3.0/immune_alignment.html



# PART 1: DATA IMPORTATION -----------------------------------------

#1. Import the data. 
adj.matrix <- Read10X("Input/filtered_feature_bc_matrix") #Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X. 

#2. Create a Seurat object. 
srat <- CreateSeuratObject(adj.matrix, project="CellSonics") 
srat
adj.matrix <- NULL #Erase adj.matrix from memory to save RAM. 

#3. Look at the Seurat object a bit closer. str command allows us to see all fields of the class.
str(srat) 
meta <- srat@meta.data #Meta.data is the most important field for next steps. 
dim(meta) 
head(meta) #Right now it shows per cell: dataset ID, number of UMI reads detected per cell (nCount_RNA), and the number of expressed (detected) genes per same cell (nFeature_RNA).



# PART 2: QUALITY CONTROL AND CELL FILTERING -----------------------------------------

#1. Do a basic quality control (QC). 
summary(meta$nCount_RNA) #UMI reads detected per cell. 
summary(meta$nFeature_RNA) #Number of expressed (detected) genes per same cell.
##Calculate the proportion of transcripts mapping to mitochondrial genes. 
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "mt-") #Or pattern = "^MT-" if human dataset. 
##Calculate the proportion of transcripts mapping to ribosomal proteins. 
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]") #Or "^RP[SL]".
head(srat[[]])

#2. Create violin plots of the selected metadata features. 
VlnPlot(srat, features=c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol=4,pt.size=0.1) & 
  theme(plot.title = element_text(size=10))

#3. Plot some of the metadata features against each other and see how they correlate (PCC shown).
p1<-FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt");p1
p2<-FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");p2
p3<-FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb");p3
p4<-FeatureScatter(srat, feature1 = "percent.rb", feature2 = "percent.mt");p4
p1 + p2
#High MT percentage strongly correlates with low UMI counts, and usually is interpreted as dead cells. 
#High ribosomal protein content strongly anti-correlates with MT, and seems to contain biological signal. 

#4. Set a QC column in metadata and define it in an informative way.
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500,'Low_nFeature', 'Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']]) #Most cells (26.0007) pass the quality control. 

#5. plot metadata only for cells that pass the tentative QC.
VlnPlot(subset(srat, subset = QC == 'Pass'), 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol=4, pt.size=0.1) & 
  theme(plot.title = element_text(size=10))

#6. Remove unwanted cells from the dataset.
srat <- subset(srat, subset = QC == 'Pass') 



# PART 3: NORMALIZATION AND DIMENSIONALITY REDUCTION -----------------------------------------

#1.Normalize the data to account for sequencing depth. 
##Conventional way is to scale it to 10,000 (as if all cells have 10k UMIs overall), and log-transform the obtained values. 
##Method: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor (10.000 default). This is then natural-log transformed using log1p.
##Normalized data are stored in srat[['RNA']]@data of the ‘RNA’ assay.
srat <- NormalizeData(srat) 

#2. Identify the most variable features (genes). Focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat), 10) #Identify the 10 most highly variable genes.
top10 

#3. Plot variable features with labels.
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

#4. Convert normalized gene expression to Z-score (values centered at 0 and with variance of 1). 
##This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate.
##The results are stored in srat[['RNA']]@scale.data and used in following PCA. 
##Default is to run scaling only on variable genes.
all.genes <- rownames(srat)
srat <- ScaleData(srat) #Add features = all.genes to scale all the genes rather than the 2.000 with highest variance.

#5. Do PCA for linear dimensionality reduction. By default we use the 2000 most variable genes.
srat <- RunPCA(srat, features = VariableFeatures(object = srat))

#6. Examine and visualize PCA results a few different ways.
print(srat[["pca"]], dims=1:5, nfeatures=5) #PC “loadings” should match markers of distinct populations for well behaved datasets.
VizDimLoadings(srat, dims=1:9, reduction="pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
VizDimLoadings(srat, dims = 1:2, reduction = "pca")
DimHeatmap(srat, dims = 1:9, cells = 500, nfeatures = 10, balanced = TRUE) 
##Can be useful when trying to decide which PCs to include for further downstream analyses. 
##In the heatmap both cells and features are ordered according to their PCA scores. 
##Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which speeds plotting for large datasets. 
##Though clearly a supervised analysis, this is a valuable tool for exploring correlated feature sets.

#7. Visualize all reduced representations. 
DimPlot(srat, reduction = "pca")

#8. Identify how many principal components should we include for the clustering.
##Rank principle components based on the percentage of variance explained by each one. 
##Where an ‘elbow’ is observed, suggests that the majority of true signal is captured in the first "n" PCs.
ElbowPlot(srat) #Answer: 10-15. 



# PART 4: CLUSTERING (can start here directly at step 4) -----------------------------------------

#1. Do clustering. Higher resolution leads to more clusters (default is 0.8). 
srat <- FindNeighbors(srat, dims = 1:15) #Construct a KNN graph based on the euclidean distance in PCA space and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). 
srat <- FindClusters(srat, resolution = 0.5) #Apply modularity optimization techniques such as the Louvain algorithm (default) to iteratively group cells together, with the goal of optimizing the standard modularity function. The resolution parameter sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
head(Idents(srat), 5) #Look at cluster IDs of the first 5 cells.
table(srat@meta.data$seurat_clusters) #Look at cluster sizes.

#2. For visualization purposes, generate UMAP reduced dimensionality representation.
srat <- RunUMAP(srat, dims = 1:15, verbose = F)
DimPlot(srat,label.size = 4,repel = T,label = T) #PDF 5x7. #DimPlot uses UMAP by default, with Seurat clusters as identity.

#3. Save the Seurat object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(srat, file = "Output/Seurat_object_v1.rds")

#4. Reimport the object.
srat <- NULL
srat <- readRDS(file = "Output/Seurat_object_v1.rds")

#5. In order to control for clustering resolution and other possible artifacts, visualize some possible cofounders.
FeaturePlot(srat, features = c("percent.mt", "percent.rb", "nFeature_RNA", "nCount_RNA")) #PDF 8x11.
VlnPlot(srat,features = c("percent.mt", "percent.rb", "nFeature_RNA", "nCount_RNA"), ncol=2, pt.size=0.01) & 
  theme(plot.title = element_text(size=10)) #PDF 8x11.

#6. Let’s remove the cells that did not pass QC and compare plots.
DimPlot(srat,label.size = 4,repel = T,label = T)
srat <- subset(srat, subset = QC == 'Pass')
DimPlot(srat,label.size = 4,repel = T,label = T) #Nothing changes coz we had already removed them. 



# PART 5: SCTRANSFORM NORMALIZATION -----------------------------------------

#1. Since we have performed extensive QC with empty cell removal, we can now apply SCTransform normalization. 
##Single SCTransform command replaces NormalizeData, ScaleData, and FindVariableFeatures. 
##This will correct for confounding factors (i.e. % rb genes) that change dramatically between clusters using vars.to.regres.
##However, I won't do it because Jimmie Ye (expert from UCSF) told that the % rb genes didn't seem an artifact, so no correction should be applied. 
#library("glmGamPoi")
#srat <- SCTransform(srat, method = "glmGamPoi", ncells = 26007, vars.to.regress = c("percent.rb"), verbose = F) 
#srat 



# PART 6: COLOURING BY GENES -----------------------------------------

#1. Let's get a very crude idea of what cell clusters are. 

##Plot some canonical markers. 
FeaturePlot(srat, features = c("Ms4a1", #"Ms4a1: B cells"
                               "Lyz1", #"Lyz1: monocytes"
                               "Cd14", #"Cd14: monocytes"
                               "Nkg7", #"Nkg7: natural killers"
                               "Cd8b1", #"Cd8b1: CD8 T cells"
                               "Cd8a", #"Cd8a: CD8 T cells"
                               "Il7r", #"Il7r: CD4 T cells"
                               "Ptprc", #"Ptprc (CD45)"
                               "Adgre1", #"Adgre1: monocytes-macrophages"
                               "Mpo", 
                               "Ccl5", 
                               "Lgals1")) #PDF 10x17

##Plot some other genes that we & Jimmie Ye were curious to know about.
all.genes <- as.data.frame(rownames(srat))
FeaturePlot(srat, features = c("Gzmb", "Gzmk", "Gzmm", "Prf1", "Pdcd1", 
                               "Cd274", "Lag3", "Havcr2", "Lair1", 
                               "Serpina1a", "Cd47")) #PDF 10x17
FeaturePlot(srat, features = c("Myc", "Klf4", "Pou5f1", "Oct4", "Sox2")) #PDF 8x11.



# PART 7: DIFFERENTIAL EXPRESSION AND MARKER SELECTION -----------------------------------------

#1. Find markers for every cluster by comparing it to all remaining cellS. By default, Wilcoxon Rank Sum test is used. 
##For speed, increase the minimal percentage and log2FC cutoffs.
all.markers <- FindAllMarkers(srat, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)

#2. Take a quick glance at the markers.
dim(all.markers) 
table(all.markers$cluster)
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers
top5_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers
##Save as an excel.
write_xlsx(top5_markers,"Output/top5_markers_15_clusters.xlsx")

#3. Visualize markers expression. 
VlnPlot(srat, features = c("Ms4a1", "Cd79a", "Cd3e"))
VlnPlot(srat, features = c("Nkg7", "Pf4"), slot = "counts", log = TRUE) #You can plot raw counts as well
FeaturePlot(srat, features = c("Ms4a1", "Gnly", "Cd3e", "Cd14", "Fcer1a", 
                               "Fcgr3a", "Lyz", "Ppbp", "Cd8a", "Lgals1", 
                               "Igfbp4", "Hmgn2", "Spp1", "Ccl5", "Mpo"))

#4. Generate an expression heatmap for the top 5 markers (or all markers if less than 5) for each cluster.
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(srat, features = top5$gene) + NoLegend() #PDF 9x18.

#5. Visualize across clusters the expression levels of markers and the % of cells expressing them. 
markers.to.plot_usual.ones <- c("Ms4a1", "Gnly", "Cd3e", "Cd14", "Fcer1a", 
                                "Fcgr3a", "Lyz", "Ppbp", "Cd8a", "Lgals1", 
                                "Igfbp4", "Hmgn2", "Spp1", "Ccl5", "Mpo")
markers.to.plot_top3.ones <- unique(top3_markers$gene)
srat <- SetIdent(srat, value = "seurat_clusters")
DotPlot(srat, features = markers.to.plot_usual.ones, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
DotPlot(srat, features = markers.to.plot_top3.ones, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis() #PDF 6x15.



# PART 8: CELL TYPE ANNOTATION USING SINGLER (can start here directly at step 8) -----------------------------------------

#1. Let’s get reference datasets from celldex package. 
mouse.ref <- celldex::ImmGenData() #The one recommended by Jimmie Ye. 

#2. Convert our Seurat object to single cell experiment (SCE) for convenience. 
sce <- as.SingleCellExperiment(DietSeurat(srat))
sce

#3. Assign cell types to the clusters using SingleR. 
##Note that there are two cell type assignments, label.main and label.fine. WARNING: Very slow. 
mouse.main <- SingleR(test = sce,assay.type.test = 1,ref = mouse.ref,labels = mouse.ref$label.main)
mouse.fine <- SingleR(test = sce,assay.type.test = 1,ref = mouse.ref,labels = mouse.ref$label.fine)

#4. See the summary of general cell type annotations. 
table(mouse.main$pruned.labels)
table(mouse.fine$pruned.labels)

#5. Add the annotations to the Seurat object metadata.
srat@meta.data$mouse.main_ImmGen <- mouse.main$pruned.labels
srat@meta.data$mouse.fine_ImmGen <- mouse.fine$pruned.labels

#6. Visualize the cell type annotations.
srat <- SetIdent(srat, value = "mouse.main_ImmGen")
p1<-DimPlot(srat, label=T , repel=T, label.size=4, cols="Set3", pt.size=0.6) + NoLegend(); p1 #PDF 5x7 
srat <- SetIdent(srat, value = "mouse.fine_ImmGen")
p2<-DimPlot(srat,  label=T , repel=T, label.size=3, cols="Set3", pt.size=0.6) + NoLegend(); p2 #PDF 5x7
p1+p2 #PDF 6.5x15

#7. Save the Seurat object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(srat, file = "Output/Seurat_object_v2.rds")

#8. Reimport the object.
srat <- NULL
srat <- readRDS(file = "Output/Seurat_object_v2.rds")



# PART 9: DGE BETWEEN DISSOCIATION METHODS -----------------------------------------

#1. Create a metadata column indicating sample ID.
meta <- srat@meta.data
meta$sampleID <- rownames(meta)
meta$sampleID <- substr(meta$sampleID, nchar(meta$sampleID) - 1 + 1, nchar(meta$sampleID)) #Extract the last character.
table(meta$sampleID)
srat@meta.data$sampleID <- meta$sampleID
srat <- SetIdent(srat, value = "sampleID")
DimPlot(srat, label = F , repel = T, split.by = "sampleID", pt.size=0.7) + NoLegend() #PDF 5x15

#2. Create a metadata column indicating dissociation technique.
meta$dissociation <- ifelse(meta$sampleID=="1" | meta$sampleID=="2", "Enzymatic", "SimpleFlow")
table(meta$dissociation)                           
srat@meta.data$dissociation <- meta$dissociation
srat <- SetIdent(srat, value = "dissociation")
DimPlot(srat, label=F , repel=T, label.size = 3, shuffle=T) + NoLegend() #PDF 5x7

#3. The next steps show how to identify what genes change in different conditions for cells of the same type. 
##Link to the guidelines used: https://satijalab.org/seurat/archive/v3.0/immune_alignment.html

#4. Create a column in the metadata to hold both the cell type and dissociation. Then switch the current ident to that column. 
srat <- SetIdent(srat, value = "mouse.main_ImmGen")
srat$celltype.arm <- paste(Idents(srat), srat$dissociation, sep = "_")
srat$celltype <- Idents(srat)
Idents(srat) <- "celltype.arm"
table(srat$celltype.arm)

#5. Use FindMarkers to find the genes that are different between enzymatic and SimpleFlow main cell types. 
table(srat$celltype)
table(srat$celltype.arm)
##Neutrophils
Neutrophils <- FindMarkers(srat, ident.1 = "Neutrophils_Enzymatic", ident.2 = "Neutrophils_SimpleFlow", verbose = FALSE)
head(Neutrophils, n = 15)
##Tgd
Tgd <- FindMarkers(srat, ident.1 = "Tgd_Enzymatic", ident.2 = "Tgd_SimpleFlow", verbose = FALSE)
head(Tgd, n = 15)
##Monocytes
Monocytes <- FindMarkers(srat, ident.1 = "Monocytes_Enzymatic", ident.2 = "Monocytes_SimpleFlow", verbose = FALSE)
head(Monocytes, n = 15)
##Stem cells
Stem.cells <- FindMarkers(srat, ident.1 = "Stem cells_Enzymatic", ident.2 = "Stem cells_SimpleFlow", verbose = FALSE)
head(Stem.cells, n = 15)

#6. Create violin plots of some of the identified genes.
plots <- VlnPlot(srat, features = c("Hspa1b", "Dnajb1", "Nr4a1", "Hspa8"), 
                 split.by = "dissociation", group.by = "celltype", 
                 pt.size = 0, combine=F, split.plot=T)
CombinePlots(plots = plots, ncol = 2) #PDF 8x15.

#7. Identify the pathways that are different between enzymatic and SimpleFlow cells. 
##Copy the identified top 15 genes from each of the main cell types. 
##Import them into the DAVID webpage to identify their pathways. 



# PART 10: CELL PROPORTIONS DIFFERENCES BETWEEN DISSOCIATION METHODS -----------------------------------------

#1. Compare the proportions of cell types among samples. 
table(meta$dissociation, meta$mouse.main_ImmGen)
table(meta$sampleID, meta$mouse.main_ImmGen)
table(meta$sampleID)

#2. Do stacked and grouped barplots of cell type (main) frequency per sample. 

##Create a dataset with the counts. 
counts <- data.frame(table(meta$sampleID, meta$mouse.main_ImmGen, meta$dissociation))
counts$sample.ID <- counts$Var1
counts$cell.type <- counts$Var2
counts$dissociation <- counts$Var3
counts$Var1 <- NULL
counts$Var2 <- NULL
counts$Var3 <- NULL

##Count the frequency of each cell type to order them. 
a <- data.frame(table(meta$mouse.main_ImmGen))
a$mouse.main_ImmGen <- a$Var1
a$mouse.main_ImmGen <- reorder(a$mouse.main_ImmGen, a$Freq) 
levels(a$mouse.main_ImmGen)
cell.type.order <- as.list(levels(a$mouse.main_ImmGen))
counts$cell.type <- as.factor(counts$cell.type)
counts$cell.type <- factor(counts$cell.type, cell.type.order)
levels(counts$cell.type)

##Plot the barplots. 
p1<-ggplot(counts, aes(fill=sample.ID, y=Freq, x=cell.type)) +
  geom_bar(position="dodge", stat="identity") + #position = "fill" for percentages, "stack" for counts, "dodge" for grouped.
  labs(x = "Cell type", y = "Number of cells", fill = "Sample ID") +
  rremove("legend")+
  theme(axis.text.x=element_text(angle=45,hjust=1));p1
  #scale_fill_manual(values=c("indianred1", "palegreen3", "skyblue2", "orchid2")) 

##Plot %s (out of total cells per sample) instead of cell numbers. 
table(meta$sampleID)
counts$percent <- ifelse(counts$sample.ID=="1", counts$Freq / 4596 * 100, 
                         ifelse(counts$sample.ID=="2", counts$Freq / 6610 * 100, 
                         ifelse(counts$sample.ID=="3", counts$Freq / 6458 * 100, 
                         ifelse(counts$sample.ID=="4", counts$Freq / 8343 * 100, NA))))
p2<-ggplot(counts, aes(fill=sample.ID, y=percent, x=cell.type)) +
  geom_bar(position="dodge", stat="identity") + #position = "fill" for percentages, "stack" for counts, "dodge" for grouped.
  labs(x = "Cell type", y = "% of cells (out of sample)", fill = "Sample ID") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  rremove("legend")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01));p2 #PDF 5x10
grid.arrange(p1, p2, nrow = 2) #PDF 8x8

#3. Do stacked and grouped barplots of cell type (fine) frequency per sample. 

##Create a dataset with the counts. 
counts <- data.frame(table(meta$sampleID, meta$mouse.fine_ImmGen, meta$dissociation))
counts$sample.ID <- counts$Var1
counts$cell.type <- counts$Var2
counts$dissociation <- counts$Var3
counts$Var1 <- NULL
counts$Var2 <- NULL
counts$Var3 <- NULL

##Count the frequency of each cell type to order them. 
a <- data.frame(table(meta$mouse.fine_ImmGen))
a$mouse.main_ImmGen <- a$Var1
a$mouse.main_ImmGen <- reorder(a$mouse.main_ImmGen, a$Freq) 
levels(a$mouse.main_ImmGen)
cell.type.order <- as.list(levels(a$mouse.main_ImmGen))
counts$cell.type <- as.factor(counts$cell.type)
counts$cell.type <- factor(counts$cell.type, cell.type.order)
levels(counts$cell.type)

##Plot the barplots. 
ggplot(counts, aes(fill=sample.ID, y=Freq, x=cell.type)) +
  geom_bar(position="dodge", stat="identity") + #position = "fill" for percentages, "stack" for counts, "dodge" for grouped.
  labs(x = "Cell type", y = "Number of cells", fill = "Sample ID") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) #PDF 5x25
  #scale_fill_manual(values=c("indianred1", "palegreen3", "skyblue2", "orchid2")) 

##Group the least frequent cell types together. 
levels(counts$cell.type)
lowest.cells <- as.list(head(levels(counts$cell.type), n=110))
counts$cell.type2 <- ifelse(counts$cell.type %in% lowest.cells, "Other", as.character(counts$cell.type))
counts$cell.type2 <- factor(counts$cell.type2, cell.type.order)
p3<-ggplot(counts, aes(fill=sample.ID, y=Freq, x=cell.type2)) +
  geom_bar(position="dodge", stat="identity") + #position = "fill" for percentages, "stack" for counts, "dodge" for grouped.
  rremove("legend")+
  labs(x = "Cell type", y = "Number of cells", fill = "Sample ID") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7));p3 #PDF 5x10

##Plot %s (out of total cells per sample) instead of cell numbers. 
table(meta$sampleID)
counts$percent <- ifelse(counts$sample.ID=="1", counts$Freq / 4596 * 100, 
                         ifelse(counts$sample.ID=="2", counts$Freq / 6610 * 100, 
                         ifelse(counts$sample.ID=="3", counts$Freq / 6458 * 100, 
                         ifelse(counts$sample.ID=="4", counts$Freq / 8343 * 100, NA))))
p4<-ggplot(counts, aes(fill=sample.ID, y=percent, x=cell.type2)) +
  geom_bar(position="dodge", stat="identity") + #position = "fill" for percentages, "stack" for counts, "dodge" for grouped.
  labs(x = "Cell type", y = "% of cells (out of sample)", fill = "Sample ID") +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7)) +
  rremove("legend")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01));p4 #PDF 5x10
grid.arrange(p3, p4, nrow = 2) #PDF 8x8
grid.arrange(p1, p3, p2, p4, nrow = 2) #PDF 8x12



# PART 11: SANKEY DIAGRAMS OF CELL TYPES -----------------------------------------

#1. Create a dataset with cells as rows and the columns "main cell type" and "fine cell type".
metadata <- srat[[]]
metadata$cellID <- rownames(metadata)
names(metadata)
SankeyDF <- metadata[,c("mouse.main_ImmGen", "mouse.fine_ImmGen", "seurat_clusters"), ]

#2. Edit the variables. 
table(SankeyDF$mouse.main_ImmGen)
table(SankeyDF$mouse.fine_ImmGen)

##Simplify the main labels. 
SankeyDF$main <- SankeyDF$mouse.main_ImmGen
a <- data.frame(table(metadata$mouse.main_ImmGen))
a$mouse.main_ImmGen <- a$Var1
a$mouse.main_ImmGen <- reorder(a$mouse.main_ImmGen, a$Freq) 
levels(a$mouse.main_ImmGen)
cell.type.order <- as.list(levels(a$mouse.main_ImmGen))
SankeyDF$main <- as.factor(SankeyDF$main)
SankeyDF$main <- factor(SankeyDF$main, cell.type.order)
levels(SankeyDF$main)
table(SankeyDF$main)
lowest.cells <- as.list(head(levels(SankeyDF$main), n=11)) #Remove last 11 levels.
SankeyDF$main <- ifelse(SankeyDF$main %in% lowest.cells, "Other", as.character(SankeyDF$main))
table(SankeyDF$main)

##Simplify the fine labels.
SankeyDF$fine <- gsub("\\s*\\([^\\)]+\\)","",as.character(SankeyDF$mouse.fine_ImmGen)) #Remove part under (). 
table(SankeyDF$fine)  
SankeyDF$mouse.fine_ImmGen <- as.factor(SankeyDF$mouse.fine_ImmGen)
a <- data.frame(table(SankeyDF$fine))
a$mouse.main_ImmGen <- a$Var1
a$mouse.main_ImmGen <- reorder(a$mouse.main_ImmGen, a$Freq) 
levels(a$mouse.main_ImmGen)
cell.type.order <- as.list(levels(a$mouse.main_ImmGen))
SankeyDF$fine <- as.factor(SankeyDF$fine)
SankeyDF$fine <- factor(SankeyDF$fine, cell.type.order)
levels(SankeyDF$fine)
table(SankeyDF$fine)
lowest.cells <- as.list(head(levels(SankeyDF$fine), n=11)) #Remove last 11 levels.
SankeyDF$fine <- ifelse(SankeyDF$fine %in% lowest.cells, "Other", as.character(SankeyDF$fine))
table(SankeyDF$fine)

##Simplify the fine labels keeping details from cell types with >500 cells.
SankeyDF$mouse.fine_ImmGen <- as.character(SankeyDF$mouse.fine_ImmGen)
SankeyDF$fine_Tdetails <- ifelse(SankeyDF$fine=="T cells", SankeyDF$mouse.fine_ImmGen, SankeyDF$fine)
table(SankeyDF$fine_Tdetails)
SankeyDF$fine_Tdetails <- ifelse(SankeyDF$mouse.fine_ImmGen=="T cells (T.DN2)" | 
                                 SankeyDF$mouse.fine_ImmGen=="T cells (T.DN2A)" | 
                                 SankeyDF$mouse.fine_ImmGen=="T cells (T.DN2B)" |
                                 SankeyDF$mouse.fine_ImmGen=="Stem cells (GMP)" |
                                 SankeyDF$mouse.fine_ImmGen=="Neutrophils (GN)" |
                                 SankeyDF$mouse.fine_ImmGen=="Neutrophils (GN.ARTH)" |
                                 SankeyDF$mouse.fine_ImmGen=="Neutrophils (GN.URAC)" |
                                 SankeyDF$mouse.fine_ImmGen=="Monocytes (MO.6C+II-)",
                                 SankeyDF$mouse.fine_ImmGen, SankeyDF$fine)
SankeyDF$fine_Tdetails <- ifelse(SankeyDF$fine_Tdetails=="T cells", "T cells (other)",
                                 ifelse(SankeyDF$fine_Tdetails=="Stem cells", "Stem cells (other)",
                                 ifelse(SankeyDF$fine_Tdetails=="Neutrophils", "Neutrophils (other)",
                                 ifelse(SankeyDF$fine_Tdetails=="Monocytes", "Monocytes (other)",
                                 SankeyDF$fine_Tdetails))))
table(SankeyDF$fine_Tdetails)                                 
                     
#3. Create a dataset with the frequency of each combination of the 2 variables. 
##Complex one.
SankeyDF1 <- SankeyDF[,c("mouse.main_ImmGen", "mouse.fine_ImmGen"), ]
SankeyDF1 <- as.data.frame(table(SankeyDF1))
##Simplified one.
SankeyDF2 <- SankeyDF[,c("main", "fine"), ]
SankeyDF2 <- as.data.frame(table(SankeyDF2))
##Intermediate one.
SankeyDF3 <- SankeyDF[,c("main", "fine_Tdetails"), ]
SankeyDF3 <- as.data.frame(table(SankeyDF3))
##Cluster one. 
SankeyDF4 <- SankeyDF[,c("seurat_clusters", "fine_Tdetails"), ]
SankeyDF4 <- as.data.frame(table(SankeyDF4))

#4. Create Sankey diagrams.
##Complex one.
ggplot(as.data.frame(SankeyDF1), aes(y = Freq, axis1 = mouse.main_ImmGen, axis2 = mouse.fine_ImmGen)) +
  geom_alluvium(aes(fill = mouse.main_ImmGen), width = 1/12) +
  geom_stratum(width = 1/3, color="grey", aes(fill=mouse.main_ImmGen)) +
  geom_label(stat = "stratum",  aes(label = after_stat(stratum))) 
##Simplified one. 
ggplot(data = as.data.frame(SankeyDF2),
       aes(axis1 = main, axis2 = fine, y = Freq)) +
    geom_alluvium(aes(fill = main)) +
    geom_stratum() +
    #geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_label(stat = "stratum",  aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("main", "fine"), expand = c(0.15, 0.05)) +
    theme_void() #PDF 10x15.
##Intermediate one. 
ggplot(data = as.data.frame(SankeyDF3),
       aes(axis1 = main, axis2 = fine_Tdetails, y = Freq)) +
  geom_alluvium(aes(fill = main)) +
  geom_stratum() +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_label(stat = "stratum",  aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("main", "fine"), expand = c(0.15, 0.05)) +
  theme_void() #PDF 10x15.
##Cluster one. 
ggplot(data = as.data.frame(SankeyDF4),
       aes(axis1 = seurat_clusters, axis2 = fine_Tdetails, y = Freq)) +
  geom_alluvium(aes(fill = seurat_clusters)) +
  geom_stratum() +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_label(stat = "stratum",  aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("main", "fine"), expand = c(0.15, 0.05)) +
  theme_void() #PDF 10x12.


