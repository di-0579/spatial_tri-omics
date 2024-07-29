rm(list=ls())
library(Seurat)
library(grid)
library(ggplot2)
library(R.matlab)
library(patchwork)
library(Matrix)
library(data.table)
set.seed(123)

getwd()
setwd("./mergrRNA_11smaple/")



library(SeuratObject)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(Matrix)
library(harmony)

spatial.Q43P0 <- readRDS("/results/20240215/Q43_P0_4/object/00sptialRepairRNA_brain.rds")
spatial.Q43P2 <- readRDS("/results/20240215/Q43_2/object/00sptialRepairRNA_brain.rds")
spatial.Q43P5 <- readRDS("/results/20240215/Q43_P5/object/00sptialRepairRNA_brain.rds")
spatial.Q43P10 <- readRDS("/results/20240215/Q43e_P10/object/00sptialRepairRNA_brain.rds")
spatial.Q43P21 <- readRDS("/results/20240215/R44_P21/object/00sptialRepairRNA_brain.rds")




spatial.Q43P0@meta.data$Age <- "P0"
spatial.Q43P0@meta.data$sampleNames <- "Q43_P0"
spatial.Q43P0@meta.data$SampleID <- "C_series_P0"



spatial.Q43P2@meta.data$Age <- "P2"
spatial.Q43P2@meta.data$sampleNames <- "Q43_P2"
spatial.Q43P2@meta.data$SampleID <- "C_series_P2"




spatial.Q43P5@meta.data$Age <- "P5"
spatial.Q43P5@meta.data$sampleNames <- "Q43_P5"
spatial.Q43P5@meta.data$SampleID <- "C_series_P5"



spatial.Q43P10@meta.data$Age <- "P10"
spatial.Q43P10@meta.data$sampleNames <- "Q43_P10"
spatial.Q43P10@meta.data$SampleID <- "C_series_P10"



spatial.Q43P21@meta.data$Age <- "P21"
spatial.Q43P21@meta.data$sampleNames <- "Q43_P21"
spatial.Q43P21@meta.data$SampleID <- "C_series_P21"


mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.Q43P0), value = TRUE)
mito.genes
spatial.Q43P0[["percent.mito"]] <- PercentageFeatureSet(spatial.Q43P0, features = mito.genes)
DefaultAssay(spatial.Q43P0) <- "Spatial"
spatial.Q43P0

mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.Q43P2), value = TRUE)
mito.genes
spatial.Q43P2[["percent.mito"]] <- PercentageFeatureSet(spatial.Q43P2, features = mito.genes)
DefaultAssay(spatial.Q43P2) <- "Spatial"
spatial.Q43P2

mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.Q43P5), value = TRUE)
mito.genes
spatial.Q43P5[["percent.mito"]] <- PercentageFeatureSet(spatial.Q43P5, features = mito.genes)
DefaultAssay(spatial.Q43P5) <- "Spatial"
spatial.Q43P5

mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.Q43P10), value = TRUE)
mito.genes
spatial.Q43P10[["percent.mito"]] <- PercentageFeatureSet(spatial.Q43P10, features = mito.genes)
DefaultAssay(spatial.Q43P10) <- "Spatial"
spatial.Q43P10

mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.Q43P21), value = TRUE)
mito.genes
spatial.Q43P21[["percent.mito"]] <- PercentageFeatureSet(spatial.Q43P21, features = mito.genes)
DefaultAssay(spatial.Q43P21) <- "Spatial"
spatial.Q43P21



###################################################################################################################

spatial.C8P0 <- readRDS("/P0_C8_4/object/00sptialRepairRNA_brain.rds")
spatial.B8P0 <- readRDS("/P0_B8_4/object/00sptialRepairRNA_brain.rds")
spatial.B8P2 <- readRDS("/P2_B8_2/object/00sptialRepairRNA_brain.rds")
spatial.C8P5 <- readRDS("/P5_C8_3/object/00sptialRepairRNA_brain.rds")
spatial.C8P10 <- readRDS("/P10_C8_5/object/00sptialRepairRNA_brain.rds")
spatial.C8P22 <- readRDS("/P21/object/00sptialRepairRNA_brain.rds")




spatial.C8P0@meta.data$Age <- "P0"
spatial.C8P0@meta.data$sampleNames <- "C8_P0"
spatial.C8P0@meta.data$SampleID <- "C_series_P0"

mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.C8P0), value = TRUE)
mito.genes
spatial.C8P0[["percent.mito"]] <- PercentageFeatureSet(spatial.C8P0, features = mito.genes)
DefaultAssay(spatial.C8P0) <- "Spatial"
spatial.C8P0

spatial.B8P0@meta.data$Age <- "P0"
spatial.B8P0@meta.data$sampleNames <- "B8_P0"
spatial.B8P0@meta.data$SampleID <- "B_series_P0"

mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.B8P0), value = TRUE)
mito.genes
spatial.B8P0[["percent.mito"]] <- PercentageFeatureSet(spatial.B8P0, features = mito.genes)
DefaultAssay(spatial.B8P0) <- "Spatial"
spatial.B8P0




spatial.B8P2@meta.data$Age <- "P2"
spatial.B8P2@meta.data$sampleNames <- "B8_P2"
spatial.B8P2@meta.data$SampleID <- "B_series_P2"



mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.B8P2), value = TRUE)
mito.genes
spatial.B8P2[["percent.mito"]] <- PercentageFeatureSet(spatial.B8P2, features = mito.genes)
DefaultAssay(spatial.B8P2) <- "Spatial"
spatial.B8P2



spatial.C8P5@meta.data$Age <- "P5"
spatial.C8P5@meta.data$sampleNames <- "C8_P5"
spatial.C8P5@meta.data$SampleID <- "C_series_P5"


mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.C8P5), value = TRUE)
mito.genes
spatial.C8P5[["percent.mito"]] <- PercentageFeatureSet(spatial.C8P5, features = mito.genes)
DefaultAssay(spatial.C8P5) <- "Spatial"
spatial.C8P5



spatial.C8P10@meta.data$Age <- "P10"
spatial.C8P10@meta.data$sampleNames <- "C8_P10"
spatial.C8P10@meta.data$SampleID <- "C_series_P10"



mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.C8P10), value = TRUE)
mito.genes
spatial.C8P10[["percent.mito"]] <- PercentageFeatureSet(spatial.C8P10, features = mito.genes)
DefaultAssay(spatial.C8P10) <- "Spatial"
spatial.C8P10



spatial.C8P22@meta.data$Age <- "P21"
spatial.C8P22@meta.data$sampleNames <- "C8_P21"
spatial.C8P22@meta.data$SampleID <- "C_series_P21"





mito.genes <- grep(pattern = "^mt-", x = rownames(spatial.C8P22), value = TRUE)
mito.genes
spatial.C8P22[["percent.mito"]] <- PercentageFeatureSet(spatial.C8P22, features = mito.genes)
DefaultAssay(spatial.C8P22) <- "Spatial"
spatial.C8P22



spatial.merged <- merge(spatial.Q43P0,c(spatial.C8P0,spatial.B8P0,spatial.B8P2,spatial.Q43P2,
                                        spatial.Q43P5,spatial.C8P5,spatial.C8P10,
                                        spatial.Q43P10, spatial.Q43P21,spatial.C8P21))


head(colnames(spatial.merged))
tail(colnames(spatial.merged))
spatial.merged




spatial.merged <- spatial.merged %>%
  NormalizeData(assay = "Spatial") %>%
  FindVariableFeatures(selection.method = "vst", assay = "Spatial",nfeatures = 3000) %>% 
  ScaleData(assay = "Spatial") %>%
  SCTransform(assay = "Spatial",vars.to.regress = c("percent.mito"))

spatial.merged




spatial.merged <- RunPCA(spatial.merged, assay = "SCT", npcs = 50)


table(spatial.merged$sampleNames)

table(spatial.merged$Age)

table(spatial.merged$SampleID)



harmonized_seurat <- RunHarmony(spatial.merged, 
                                group.by.vars = c("sampleNames", "Age","SampleID"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

ElbowPlot(harmonized_seurat, ndims = 50)

names(spatial.merged@reductions)
names(harmonized_seurat@reductions)
harmonized_seurat
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:50)

harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")

names(harmonized_seurat@graphs)
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 0.8, graph.name = "SCT_snn")


table(harmonized_seurat$SCT_snn_res.0.8)


col_method <- c(
     "P0"="#FDE725FF", 
     "P2"="#D3F24BFF",  
     "LPC5_I13_4"="#AADC32FF", 
     "P5"="#6FCF97FF", 
     "LPC21_H13_2"="#27AD81FF", 
     "P10"="#1A9875FF",  
     "P21"="#3B528BFF", 
     "P22"="#440154FF"
   )

col_method <- c(
  "P22"="#FDE725FF", 
#  "P21"="#D3F24BFF",  
  "P21"="#AADC32FF", 
  "P10"="#6FCF97FF", 
#  "LPC21_H13_2"="#27AD81FF", 
  "P5"="#1A9875FF",  
  "P2"="#3B528BFF", 
  "P0"="#440154FF"
)




p2 <- DimPlot(harmonized_seurat, group.by = "Age", pt.size = 0.1, reduction = "umap",cols =col_method,  label = FALSE)
p2



sampleNames <- "harmonyRNA"
ggsave(paste0("plot/", sampleNames, "_AGE_RNAUMP.png"), plot = p2, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_AGE_RNAUMP.pdf"), plot = p2, width = 9, height = 9)






p2 <-DimPlot(harmonized_seurat, group.by = "sampleNames", pt.size = 0.1, reduction = "umap",label = FALSE)
p2 
ggsave(paste0("plot/", sampleNames, "_sampleNames_RNAUMP.png"), plot = p2, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_sampleNames_RNAUMP.pdf"), plot = p2, width = 9, height = 9)



library(grid)

library(paletteer)#
col <- paletteer_d("ggthemes::Classic_20")





table(harmonized_seurat$seurat_clusters)



p1 <- DimPlot(harmonized_seurat, group.by = "SCT_snn_res.1", pt.size = 0.1, reduction = "umap",cols = bright_colors_32,label = FALSE)
p1
p1+p2



ggsave(paste0("plot/", sampleNames, "_11RNAUMP_lab.png"), plot = p1, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_11RNAUMP.png"), plot = p1, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_11RNAUMP.pdf"), plot = p1, width = 9, height = 9)





df_spatial.merged <- harmonized_seurat@meta.data
Assays(harmonized_seurat)
DefaultAssay(harmonized_seurat) <- 'Spatial'
Idents(harmonized_seurat) <- "SCT_snn_res.1.5"

library(grid)
source("/gpfs/gibbs/project/fan/dz286/uesful/code/SpatialPlot_new.R")


table(harmonized_seurat$sampleNames)



meta.data <- as.data.frame(harmonized_seurat@meta.data)
table(meta.data$SCT_snn_res.1)
table(meta.data$sampleNames)



Q43_P21_df <- meta.data[meta.data$sampleNames=='Q43_P0',]
new_row_names <- row.names(Q43_P21_df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("_.*","", x)))
row.names(Q43_P21_df) <- new_row_names


Q43_P21_df <-as.data.frame(Q43_P21_df)


spatial.Q43P0 <- AddMetaData(spatial.Q43P0, metadata = Q43_P21_df)


library(grid)

library(paletteer)#
# col <- paletteer_d("ggthemes::Classic_20")


bright_colors_32 <- c( "0"="#e7298a","5"= "#D5F2F5", "28"= "#FFD130", "9"="#0091FF", "25"="#F868A1","16"="#000000",
                        "12"="#FF0000", "29"="#B3B3B3", "33"= "#8DD3C8","2"= "#B8FFFE", 
                       "10"= "#9CFFB2", "34"="#5CF0FE", "21"="#AEC7E8FF" ,"14"= "#00c724","23"= "#DCF2F5",#FDCDE5",
                       "19"= "#00FF00", "3"= "#3ade59", "18"= "#D9FFFE", "8"= "#FF9500","13"= "#FF836C", 
                       "1"="#313CFF", "15"="#FFF200", "22"="#689B94", "6"="#009E00", "11"="#F0F2CC", 
                       "24"= "#FFF2CC", "26"="#0059CA","120"= "#8EFF05", 
                       "20"="#FF05E1","7"="#CDFFFE", 
                       "30"="#7570b3",  "32"="#92c5de",  "4"="#00FFFF", "33"="#20B2AA",    "17"="#FFD700"
)

p51 <- SpatialPlot_new(spatial.Q43P0,group.by = 'SCT_snn_res.1', pt.size.factor = 1.1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 5)+ ggtitle('Q43P0')
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 



sampleNames <- "Q43P0"

ggsave(paste0("./spatialUMP/", sampleNames, "_spatialRNAUMP.png"), plot = p51, width = 9, height = 9)
ggsave(paste0("./spatialUMP/", sampleNames, "_spatialRNAUMP.pdf"), plot = p51, width = 9, height = 9)



saveRDS(spatial.Q43P0, file = paste0("./object/",sampleNames,"mergedspatialRNAUMP.rds"))
rm(spatial.Q43P0)




p52 <- SpatialPlot_new(spatial.Q43P21,group.by = 'SCT_snn_res.1', pt.size.factor = 1.1, label = TRUE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 5)+ ggtitle('Q43P21')
p52$layers[[1]]$aes_params <- c(p52$layers[[1]]$aes_params, shape=22) 
p52 +p51


p52 <- SpatialPlot_new(spatial.Q43P10,group.by = 'SCT_snn_res.1', pt.size.factor = 1.1, label = TRUE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 5)+ ggtitle('Q43P10')
p52$layers[[1]]$aes_params <- c(p52$layers[[1]]$aes_params, shape=22) 
p52 +p51







sampleNames <- "Q43_P21"

ggsave(paste0("./plotupdataATACCUT/", sampleNames, "_spatialADTUMP.png"), plot = p51, width = 9, height = 9)
#ggsave(paste0("plot/", sampleNames, "_spatial35RNAUMP.pdf"), plot = p51, width = 9, height = 9)


#saveRDS(spatial.LPC21_4, file = paste0("./object/",sampleNames,"ident_merged.rds"))

harmonized_seurat


sampleNames <- "harmonyRNA"

saveRDS(harmonized_seurat, file = paste0("./object/",sampleNames,"harmonyRNAUMP.rds"))
rm(harmonized_seurat)