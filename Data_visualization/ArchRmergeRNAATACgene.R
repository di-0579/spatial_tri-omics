rm(list=ls())

library(Seurat)

library(ArchR)
library(grid)

library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggplot2)
library(patchwork)
library(future)
set.seed(123)


getwd()
setwd("/gpfs/gibbs/pi/fan/dz286/results/20240215/M39N40/")


threads = 28
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- c("/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_2/LPC5_N40_2_ARC/outs/atac_fragments.tsv.gz",  
                "/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_4/LPC5_N40_4/outs/atac_fragments.tsv.gz",
                "/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_2/M39_2_ARC/outs/atac_fragments.tsv.gz",
                "/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_4/M39_4_ARC/outs/atac_fragments.tsv.gz")
sampleNames <- c('LPC5_I13_2', 'LPC5_I13_4','LPC21_H13_2','LPC21_H13_4')

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 0, #4,
  minFrags = 0, #1000,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,  minFragSize = 1,
  maxFragSize = 20000,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, 
#   knnMethod = "UMAP", corCutOff = 0.999,
#   LSIMethod = 1
# )

output_dir <- 'LPCbrain_all'

proj_ALL <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = output_dir,
  copyArrows = TRUE
)
proj_ALL


meta.data <- as.data.frame(getCellColData(ArchRProj = proj_ALL))
table(meta.data$Sample)


source("../archR_fix.R")



#Import RNA
seRNA2 <- import10xFeatureMatrix2(
  input = c("/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_2/Repair_RNA_feature_bc_matrix.h5",
            "/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_4/Repair_RNA_feature_bc_matrix.h5",
            "/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_2/Repair_RNA_feature_bc_matrix.h5",
            "/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_4/Repair_RNA_feature_bc_matrix.h5"),
  names = c('LPC5_I13_2', 'LPC5_I13_4','LPC21_H13_2','LPC21_H13_4')
)
seRNA2


seRNA <- import10xFeatureMatrix(
  input = c("/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_2/LPC5_N40_2_ARC/outs/raw_feature_bc_matrix.h5",
            "/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_4/LPC5_N40_4/outs/raw_feature_bc_matrix.h5",
            "/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_2/M39_2_ARC/outs/raw_feature_bc_matrix.h5",
            "/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_4/M39_4_ARC/outs/raw_feature_bc_matrix.h5"),
  names = c('LPC5_I13_2', 'LPC5_I13_4','LPC21_H13_2','LPC21_H13_4')
)
seRNA


assay(seRNA) <- assay(seRNA2)[match(rownames(seRNA), rownames(seRNA2)), match(colnames(seRNA), colnames(seRNA2))]



proj_ALL <- addGeneExpressionMatrix(input = proj_ALL, seRNA = seRNA, force = TRUE)

proj_ALL

getAvailableMatrices(proj_ALL)






image.N40_2 <- Read10X_Image(image.dir ="/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_2/spatial/", 
                             filter.matrix = TRUE)
bc.spatial <- row.names(image.N40_2@coordinates)
new_bc.N40_2<- paste(bc.spatial, '-1', sep = '')

image.N40_4 <- Read10X_Image(image.dir ="/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_4/spatial/", 
                          filter.matrix = TRUE)
bc.spatial <- row.names(image.N40_4@coordinates)
new_bc.N40_4<- paste(bc.spatial, '-1', sep = '')

image.M39_2 <- Read10X_Image(image.dir ="/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_2/spatial/", 
                          filter.matrix = TRUE)
bc.spatial <- row.names(image.M39_2@coordinates)
new_bc.M39_2 <- paste(bc.spatial, '-1', sep = '')

image.M39_4 <- Read10X_Image(image.dir ="/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_4/spatial/", 
                          filter.matrix = TRUE)
bc.spatial <- row.names(image.M39_4@coordinates)
new_bc.M39_4 <- paste(bc.spatial, '-1', sep = '')


sampleNames <- c('LPC5_I13_2', 'LPC5_I13_4','LPC21_H13_2','LPC21_H13_4')


cells_selected <- c(paste0('LPC5_I13_2#', new_bc.N40_2),
                    paste0('LPC5_I13_4#', new_bc.N40_4),
                    paste0('LPC21_H13_2#', new_bc.M39_2),
                     paste0('LPC21_H13_4#', new_bc.M39_4))


projHeme1 <- proj_ALL[(proj_ALL$cellNames)%in%cells_selected, ]

projHeme1
proj_ALL

getAvailableMatrices(proj_ALL)
table(projHeme1$Sample)



meta.data <- as.data.frame(getCellColData(ArchRProj = projHeme1))
table(meta.data$Sample)


projHeme1 <- projHeme1[(projHeme1$TSSEnrichment) > 0.01, ]
projHeme1
projHeme1 <- projHeme1[(projHeme1$Gex_nUMI) > 10, ]
projHeme1





saveArchRProject(ArchRProj = projHeme1, outputDirectory = paste0("Save-", output_dir), load = FALSE)
#projHeme1 <- loadArchRProject(path = paste0("Save-", output_dir), force = FALSE, showLogo = TRUE)
#projHeme1

#projHeme1 <- readRDS("/vast/palmer/scratch/fan/dz286/testATAC/testP/Save-motifProjHeme1026/Save-ArchR-Project.rds")
#projHeme1
# 
# 
# p1 <- plotGroups(
#   ArchRProj = projHeme1, 
#   groupBy = "Sample", 
#   colorBy = "cellColData", 
#   name = "TSSEnrichment",
#   plotAs = "ridges"
# )
# p1
# 
# p2 <- plotGroups(
#   ArchRProj = projHeme1, 
#   groupBy = "Sample", 
#   colorBy = "cellColData", 
#   name = "TSSEnrichment",
#   plotAs = "violin",
#   alpha = 0.4,
#   addBoxPlot = TRUE
# )
# p2
# 
# p3 <- plotGroups(
#   ArchRProj = projHeme1, 
#   groupBy = "Sample", 
#   colorBy = "cellColData", 
#   name = "log10(nFrags)",
#   plotAs = "ridges"
# )
# p3
# p4 <- plotGroups(
#   ArchRProj = projHeme1, 
#   groupBy = "Sample", 
#   colorBy = "cellColData", 
#   name = "log10(nFrags)",
#   plotAs = "violin",
#   alpha = 0.4,
#   addBoxPlot = TRUE
# )
# 
# p4
# plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 4, height = 4)
# 
# p1 <- plotFragmentSizes(ArchRProj = projHeme1)
# p1
# 
# p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
# p2
# plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
# 
# 
# 



#LSI-ATAC
projHeme1 <- addIterativeLSI(
  ArchRProj = projHeme1, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  dimsToUse = 2:30, 
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)



projHeme1 <- addHarmony(
  ArchRProj = projHeme1,
  reducedDims = "LSI_ATAC",
  name = "Harmony_ATAC",
  groupBy = "Sample",
  dimsToUse = 2:30, 
  scaleDims = TRUE,
  corCutOff = 0.60, 
  verbose = TRUE,
  force = TRUE
)


#LSI-RNA
projHeme1 <- addIterativeLSI(
  ArchRProj = projHeme1, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 3000,
  dimsToUse = 1:30,
  firstSelection = "variable",
  outlierQuantiles = c(0.02, 0.98),
  filterQuantile = 0.995,
  binarize = FALSE,
  name = "LSI_RNA"
)


projHeme1 <- addHarmony(
  ArchRProj = projHeme1,
  reducedDims = "LSI_RNA",
  dimsToUse = 1:30,
  name = "Harmony_RNA",
  groupBy = "Sample"
)




#UMAPs
projHeme1 <- addUMAP(projHeme1, reducedDims = "Harmony_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
#################################
projHeme1 <- addUMAP(projHeme1, reducedDims = "Harmony_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)



p2 <- plotEmbedding(projHeme1, name = "Sample", embedding = "UMAP_RNA", size = 0.1, labelAsFactors=F, labelMeans=F)
p2


p3 <- plotEmbedding(projHeme1, name = "Sample", embedding = "UMAP_ATAC", size = 0.1, labelAsFactors=F, labelMeans=F)
p3

projHeme1 <- addCombinedDims(projHeme1, reducedDims = c("Harmony_ATAC", "Harmony_RNA"), name =  "LSI_Combined")
#all gene
source("/gpfs/gibbs/project/fan/dz286/uesful/code/getGeneScore_ArchR2.R")
source("/gpfs/gibbs/project/fan/dz286/uesful/code/getGeneExpression_ArchR2.R")

getAvailableMatrices(projHeme1)


#UMAPs
#projHeme1 <- addUMAP(projHeme1, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
#projHeme1 <- addUMAP(projHeme1, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)


#Add ATAC_Clusters
projHeme1 <- addClusters(projHeme1, reducedDims = "Harmony_ATAC", name = "ATAC_Clusters", resolution = 1.2, force = TRUE)
table(projHeme1$ATAC_Clusters)

#Add RNA_Clusters
projHeme1 <- addClusters(projHeme1, reducedDims = "Harmony_RNA", name = "RNA_Clusters", resolution = 1.2, force = TRUE)
table(projHeme1$RNA_Clusters)
#Add Combined_Clusters
projHeme1 <- addClusters(projHeme1, reducedDims = "LSI_Combined", name = "Combined_Clusters", resolution = 1, force = TRUE)


#Plot Embedding
p1 <- plotEmbedding(projHeme1, name = "ATAC_Clusters", embedding = "UMAP_ATAC", size = 0.2, labelAsFactors=F, labelMeans=F)
p1
p2 <- plotEmbedding(projHeme1, name = "RNA_Clusters", embedding = "UMAP_RNA", size = 0.2, labelAsFactors=F, labelMeans=F)
p2
p3 <- plotEmbedding(projHeme1, name = "Combined_Clusters", embedding = "UMAP_RNA", size = 0.1, labelAsFactors=F, labelMeans=T,pal = bright_colors_25)
p3


ggsave("plot_Combined_ClustersUMAP.png", plot = p3, width = 10, height = 8, dpi = 300)
ggsave("plot__ClustersUMAP.pdf", plot = p3, width = 10, height = 8)

p3+p2

p4 <- plotEmbedding(projHeme1, name = "Sample", embedding = "UMAP_RNA", size = 0.1, labelAsFactors=F, labelMeans=F)
p4

ggsave("plot_Combined_SampleUMAP.png", plot = p4, width = 10, height = 8, dpi = 300)
ggsave("plot__SampleUMAP.pdf", plot = p4, width = 10, height = 8)



+p3




plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)



########################################

projHeme1 <- addImputeWeights(projHeme1, reducedDims = "LSI_RNA")


## Identify the marker genes for each cluster 
ATAC_markersGS <- getMarkerFeatures(
  ArchRProj = projHeme1, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Combined_Clusters",
  testMethod = "wilcoxon"
)

dd <- getMarkers(ATAC_markersGS, cutOff = "FDR <= 1")
dd
markerGenes_dd <- list()
for (i in seq_len(length(dd))) {
  markerGenes_dd <- c(markerGenes_dd, dd[[i]]$name)
}
# 
markerGenes_dd <- unlist(markerGenes_dd)
markerGenes_dd <- unique(markerGenes_dd)


projHeme1 <- addImputeWeights(projHeme1, reducedDims = "LSI_RNA")
all_ATACgene_score <- getGeneScore_ArchR(ArchRProj = projHeme1, name = markerGenes_dd, imputeWeights = getImputeWeights(projHeme1))
# # 
write.csv(all_ATACgene_score, "./data/all_ATACgene_score.csv", row.names = T, col.names = T, quote=F)










samples <- c("LPC21_H13_2", "LPC21_H13_4", "LPC5_I13_2",   "LPC5_I13_4")


for (sample in samples) {
  sample_cols <- grep(paste0("^", sample, "#"), colnames(all_ATACgene_score))
  sample_df <- all_ATACgene_score[, sample_cols]
  new_col_names <- gsub(".*#", "", colnames(sample_df))  
  new_col_names <- gsub("-.*", "", new_col_names)        
  colnames(sample_df) <- new_col_names
  write.csv(sample_df, paste0("./data/", sample, "_ATACgene_score_Combined_Clusters.csv"), row.names = TRUE, col.names = TRUE, quote = FALSE)

}

#######################################
# RNA
########################################

projHeme1 <- addImputeWeights(projHeme1, reducedDims = "LSI_RNA")


## Identify the marker genes for each cluster 
RNA_markersGS <- getMarkerFeatures(
  ArchRProj = projHeme1, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "Combined_Clusters",
  testMethod = "wilcoxon"
)

dd <- getMarkers(RNA_markersGS, cutOff = "FDR <= 1")
dd
markerGenes_dd <- list()
for (i in seq_len(length(dd))) {
  markerGenes_dd <- c(markerGenes_dd, dd[[i]]$name)
}
# 
markerGenes_dd <- unlist(markerGenes_dd)
markerGenes_dd <- unique(markerGenes_dd)


projHeme1 <- addImputeWeights(projHeme1, reducedDims = "LSI_RNA")
all_RNAgene_score <- getGeneExpression_ArchR(ArchRProj = projHeme1, name = markerGenes_dd, imputeWeights = getImputeWeights(projHeme1))
# # 
write.csv(all_RNAgene_score, "./data/all_RNAgeneExp_score.csv", row.names = T, col.names = T, quote=F)


table(projHeme1$Sample)



samples <- c("LPC21_H13_2", "LPC21_H13_4", "LPC5_I13_2",  "LPC5_I13_4")


for (sample in samples) {
  sample_cols <- grep(paste0("^", sample, "#"), colnames(all_RNAgene_score))
  sample_df <- all_RNAgene_score[, sample_cols]
  new_col_names <- gsub(".*#", "", colnames(sample_df))  
  new_col_names <- gsub("-.*", "", new_col_names)        
  colnames(sample_df) <- new_col_names
  write.csv(sample_df, paste0("./data/", sample, "_RNAexpgene_score_Combined_Clusters.csv"), row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}




saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-Combinedclusters0611", load = FALSE)


projHeme1 <- loadArchRProject("Save-Combinedclusters0611", force = FALSE, showLogo = TRUE)
projHeme1


meta.data <- as.data.frame(getCellColData(ArchRProj = projHeme1))
table(meta.data$Combined_Clusters)
table(meta.data$Sample)


samples <- c("LPC21_H13_2", "LPC21_H13_4", "LPC5_I13_2",  "LPC5_I13_4")

for (sample in samples) {
  sample_df <- meta.data[meta.data$Sample == sample, ]
  
  new_row_names <- rownames(sample_df)
  new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#", "", x)))
  new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*", "", x)))
  rownames(sample_df) <- new_row_names
  
  write.csv(sample_df, paste0("./data/", sample, "_UMPclusters_Combined.csv"), row.names = TRUE, col.names = TRUE, quote = FALSE)
}






############################################################
source("/gpfs/ycga/project/fan/dz286/useful/code/SpatialPlot_new.R")
source("/gpfs/ycga/project/fan/dz286/useful/code/getGeneScore_ArchR.R")



samples <- c("LPC21_H13_2", "LPC21_H13_4", "LPC5_I13_2",   "LPC5_I13_4")
sample <- samples[1]
sample

RNAgene_score <- read.csv(paste0("./data/", sample, "_RNAexpgene_score_Combined_Clusters.csv"),row.names = 1)
ATAC_score <- read.csv(paste0("./data/", sample, "_ATACgene_score_Combined_Clusters.csv"), row.names = 1)

meta.data <- read.csv(paste0("./data/", sample, "_UMPclusters_Combined.csv"),row.names = 1)




data.dir <-"/gpfs/gibbs/pi/fan/dz286/results/20240215/M39H13/M39_2/"
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = RNAgene_score, assay = assay, meta.data = meta.data)


image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object


# 
# p6 <- SpatialPlot_new(spatial.obj, label = FALSE, label.size = 1.8, group.by = 'RNA_Clusters', pt.size.factor = 1, cols = bright_colors_25, image.alpha = 1, stroke = 0)
# p6$layers[[1]]$aes_params <- c(p6$layers[[1]]$aes_params, shape=22)
# p6
# 
# 
# sample
# 
# 
# ggsave("LPC5_I13_4_spatialRNAUMP-ArchRs.png", plot = p6, width = 9, height = 9)





# reduced_bright_colors_13 <- c(
#   "C13" = "#7570b3",
#   "C11" = "#ff836c", 
#   "C3" = "#009e00", 
#   "C4" = "#ffe400",  
#   "C10" = "#0091ff", 
#   "C2" = "#ff0000",  
#   "C5" = "#91eda4",   
#   "C7" = "#5cf0fe",
#   "C6" = "#00ff00", 
#   "C1" = "#000000",  
#   "C12" = "#20b2aa",  
#   "C8" = "#ff9500",  
#   "C9" =  "#6e85a3"
# )
# 
# sample
# 
# p7 <- SpatialPlot_new(spatial.obj, label = TRUE, label.size = 10, group.by = 'ATAC_Clusters', pt.size.factor = 1, cols = reduced_bright_colors_13, image.alpha = 1, stroke = 0)
# p7$layers[[1]]$aes_params <- c(p7$layers[[1]]$aes_params, shape=22)
# p7

#ggsave("LPC5_I13_4_spatialATACUMP-ArchRs.png", plot = p7, width = 9, height = 9)

# table(projHeme1$Combined_Clusters)
# 
# n_clusters <- length(unique(projHeme1$Combined_Clusters))
# cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
# names(cols) <- paste0('C', seq_len(n_clusters))
# cols
# 




bright_colors_25 <- c( "C1" = "#ff0000", "C2" = "#FFF200", "C3" = "#f9e5d4", "C4" = "#20b2aa", "C5" = "#00ff00",
                       "C6" = "#91eda4", "C7" = "#7570b3", "C8" = "#0059ca", "C9" = "#ff8000", "C10" = "#313cff",
                       "C11" = "#8eff05", "C12" = "#5cf0fe", "C13" = "#cff9fa", "C14" = "#ffa6e3", "C15" = "#009e00",
                       "C16" = "#0091ff", "C17" = "#ff05e1", "C18" = "#f868a1", "C19" = "#9500ff",
                       "C21" = "#ffd9f3" ,  "C20" = "#ff836c"
 )


spatial.obj <- readRDS("./object_ArchR/LPC21_H13_4_spatialRNA_geneExp.rds")


p7 <- SpatialPlot_new(spatial.obj, label = F, label.size = 5, group.by = 'Combined_Clusters',
                      pt.size.factor = 1.2, cols = bright_colors_25, image.alpha = 1, stroke = 0)+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p7$layers[[1]]$aes_params <- c(p7$layers[[1]]$aes_params, shape=22)
p7

#sample



df <- spatial.obj@meta.data


df <- df["Combined_Clusters"]
table(df$Combined_Clusters)
sampleNames <- "M39_4_LPC21"

write.csv(df, file =paste0("./data/", sampleNames, "_Combined_Clusters.csv"), row.names = TRUE, col.names = TRUE)







ggsave("./object_ArchR/LPC5_I13_4_spatialCombined_ClustersUMP.png", plot = p7, width = 9, height = 9,dpi = 600)





features_DORC <-c("Itgax",'Siglec1',"Ptprc",#'Cd19',"Gzmb","Dnah12",'Cd1c','Cd209','Cd207','Cxcr1','Ccr2',"Cd8a","Fcer1a",'Clec9a',
                  'Cd80','Cd86','Cd14',
                  'Xcr1','Tmem119','Aif1',"Rgs9","Rarb","Mag","Mal","Slc17a7","Mef2c",
                  "Ano1","Dgkg","Flt1","Mecom","Tox3","Tacr1","Cpne4","Cd9","Itgam","Ms4a1","Dpp4","Cd36","Cd44","Ptprc","Cd48","Itgb2","Itga4","Itga6","Icam1","Cd63",
            "Fcgr1","Cd68","Vcam1","Tnfrsf1b","Il7r","Tnfrsf4","Sdc1","Siglec1","Cxcr5","Cd200","Tnfrsf14","Pdcd1",
            "Nrp1","Cx3cr1","Adgre1","Jam2","Klrg1","Cd86","Xcr1","Cx3cr1","Itgam","Itga6", "Cd9","Cd86","Sdc1","Sox10","Myrf","Itpr2","Cd44",#"H2afb1",                                      #H2AB1
  "Tbr1","Bcl11b","Fezf2",                     #TBR1   CTIP2 
  "Pdgfra",                                    #PDGFR2
  "Pax6","Sox2","Neurog2",#"Tbr2",              #PAX6
  "Olig2","Sox10","Mbp","Mag",                 #OLIGO2
  "Rbfox3",#"NeuN",                                 #NeuN
  "Mog","Plp1","Olig1", "Cnp",                        #MOG/MBP
  "Aif1","Cd68","Cx3cr1",                      #IBA1
  "Gfap","Aqp4","Slc1a3","Slc1a2","Vim","Aldh1l1",        #GFAP
  "Satb2","Sox5","Nr2f1","Nr2f2",              #CTIP2
  "Apc")






features_DORC <- "Itgax"

features_DORC

plot_features <- function(feature){
  p <- SpatialFeaturePlot(spatial.obj, features = feature, pt.size.factor = 1, 
                          image.alpha = 0, stroke = 0, alpha = c(1, 1), min.cutoff = "q2", max.cutoff = "q98") + 
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_DORC, plot_features)
ggList[[1]]



sample


saveRDS(spatial.obj, file = "./object_ArchR/LPC5_I13_4_spatialRNA_geneExp.rds")


rm(spatial.obj)

####################################################
#


#data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = ATAC_score, assay = assay, meta.data = meta.data)


image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object




source("/gpfs/ycga/project/fan/dz286/useful/code/SpatialPlot_DORC.R")

features_DORC
plot_features <- function(feature){
  p <- SpatialPlot_DORC(spatial.obj, features = feature, pt.size.factor = 1, 
                        image.alpha = 0, stroke = 0, alpha = c(1, 1), min.cutoff = "q2", max.cutoff = "q98") + 
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_DORC, plot_features)
ggList[[1]]



sample

saveRDS(spatial.obj, file = "./object_ArchR/LPC5_I13_4_spatialATAC_gene.rds")

rm(spatial.obj)

#######################################################################










# # 
#######################################
library(BSgenome.Mmusculus.UCSC.mm10)



meta.data <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))
table(meta.data$Sample)
getAvailableMatrices(projHeme1)
##################################################################################################  
proj_in_tissue <- projHeme1 
#proj_in_tissue <- readRDS("./M39N40/Save-atacrnaclusters0117/Save-ArchR-Project.rds")
proj_in_tissue
## Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Combined_Clusters")


pathToMacs2  <- "./conda_envs/R413/bin/macs2"



pathToMacs2 <- findMacs2()

table(proj_in_tissue$RNA_Clusters)

library(presto)
proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "Combined_Clusters", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

getPeakSet(proj_in_tissue)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

getAvailableMatrices(proj_in_tissue)

#devtools::install_github("GreenleafLab/chromVARmotifs")
proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)



proj_in_tissue@peakAnnotation

#if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
# proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
#}

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "PeakMatrix", 
  groupBy = "Combined_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.01")

markerList



enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj_in_tissue,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.01"
)
enrichMotifs

####################################
## ChromVAR Deviatons Enrichment
proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj_in_tissue, name = "MotifMatrix", plot = TRUE)
plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_in_tissue, addDOC = FALSE)


markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "MotifMatrix", 
  groupBy = "Combined_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)
markersMotifs 
markersMotifs@metadata


motifPositions <- getPositions(proj_in_tissue)
motifPositions



saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "Save-motifProjHeme", load = FALSE)

##############################################################################






#####################################################
markerMotifsList <- getMarkers(markersMotifs, cutOff = "FDR <= 0.1 & MeanDiff >= 0.1")
#
markerMotifsList <- getMarkers(markersMotifs, cutOff = "FDR <= 1")
# for (i in seq_len(length(markerMotifsList))) {
#   write.table(markerMotifsList[[i]], file=paste0('./data/motifs_list/', sampleNames, '_C', i, '_motifs.txt'),
#               quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
# }




motifs <- list()
for (i in seq_len(length(markerMotifsList))) {
  motifs <- c(motifs, markerMotifsList[[i]]$name)
}
motifs <- unlist(motifs)
motifs <- paste0('z:', motifs)



## Spatial plots
library(ggplot2)
library(patchwork)
library(dplyr)

source("./uesful/code/getDeviation_ArchR2.R")
source("./useful/code/SpatialPlot_new.R")


proj_in_tissue <- addImputeWeights(proj_in_tissue, reducedDims = "LSI_RNA")
dev_score <- getDeviation_ArchR(ArchRProj = proj_in_tissue, name = motifs, imputeWeights = getImputeWeights(proj_in_tissue))
dev_score[is.na(dev_score)] <- 0 #min(dev_score, na.rm = TRUE)

write.csv(dev_score, "data/xiufuallmitf_matrix.csv", row.names = T, col.names = T, quote=F)
# 

table(proj_in_tissue$Sample)


samples



for (sample in samples) {
  sample_cols <- grep(paste0("^", sample, "#"), colnames(dev_score))
  sample_df <- dev_score[, sample_cols]
  new_col_names <- gsub(".*#", "", colnames(sample_df))  
  new_col_names <- gsub("-.*", "", new_col_names)        
  colnames(sample_df) <- new_col_names
  write.csv(sample_df, paste0("./data/", sample, "_dev_scor_4sample.csv"), row.names = TRUE, col.names = TRUE, quote = FALSE)
}



###########################################################

source(./useful/code/SpatialPlot_MT.R")


samples <- c("LPC21_H13_2", "LPC21_H13_4", "LPC5_I13_2",   "LPC5_I13_4")
sample <- samples[4]
sample


Motif_score <- read.csv(paste0("./data/", sample, "_dev_scor_4sample.csv"),row.names = 1)

data.dir <-"/gpfs/gibbs/pi/fan/dz286/results/20240215/N40I13/N40_4/"
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = Motif_score, assay = assay)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image


spatial.obj <- object



############################################################
source("./useful/code/SpatialPlot_new.R")

features_DORC <- "Neurod6-786"

features_DORC

plot_features <- function(feature){
  p <- SpatialFeaturePlot(spatial.obj, features = feature, pt.size.factor = 1, 
                          image.alpha = 0, stroke = 0, alpha = c(1, 1), min.cutoff = "q2", max.cutoff = "q98") + 
    theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))+
    theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_DORC, plot_features)
ggList[[1]]



sample
saveRDS(spatial.obj, file = "./object_ArchR/LPC5_I13_4_spatial_Motif_score.rds")

rm(spatial.obj)



















#############################################################
proj_in_tissue
table(proj_in_tissue$Combined_Clusters)
########################################


proj_in_tissue <- addImputeWeights(proj_in_tissue, reducedDims = "LSI_RNA")


## Identify the marker genes for each cluster 
ATAC_markersGS <- getMarkerFeatures(
  ArchRProj = projHeme1, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Combined_Clusters",
  testMethod = "wilcoxon"
)

dd <- getMarkers(ATAC_markersGS, cutOff = "FDR <= 0.05& Log2FC >= 0.5")
dd
markerGenes_dd <- list()
for (i in seq_len(length(dd))) {
  markerGenes_dd <- c(markerGenes_dd, dd[[i]]$name)
}
# 
markerGenes_dd <- unlist(markerGenes_dd)
markerGenes_dd <- unique(markerGenes_dd)

write.csv(markerGenes_dd, "./data/ATAC_gene_list.csv", row.names = FALSE, quote = FALSE)

#gene_list_from_csv <- read.csv("ATAC_gene_list.csv", stringsAsFactors = FALSE)
#gene_list <-gene_list_from_csv$x

#markerGenes_dd <- gene_list



sampleNames <- "LPC_ATAC"


markerList_pos <- getMarkers(ATAC_markersGS, cutOff = "FDR <= 0.05& Log2FC >= 0.5")

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}

markerGenes <- unlist(markerGenes)


for (i in seq_len(length(markerList_pos))) {
  write.table(markerList_pos[[i]], file=paste0('./data/LPC_ATAC_markers_list/', sampleNames, '_C', i, '_markers.txt'),
              quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
}






markerGenes  <- c("Enpp2","Nfe2l2","Pdgfra","Olig3","Egfr","Ptprz1",
                  "Mobp","Sox10","Olig2","Mog","Olig1","Plp1","Cnp","Nkx2-2","Fa2h","Enpp2","Mbp","Mal","Tcf4",
                  "Cubn",
                  "Id3","Smad6","Notch1" ,"Zic2","Zic4","Zic5","Mertk", "Apoe","Foxl1","Foxf1","Tppp2",
                  "Nfe2l2","Tgfbr2","Cx3cr1" ,"Mertk","Csf1r","Hoxa3","Il10ra","P2ry12","Runx1","Ets1","Tnfrsf1b",
                  "Ccl3","Spi1","Pik3ap1","Prox1","Amotl1","Sox2","Foxp1","Bcl11b",
                  "Slc17a9","Satb2","Icam5","Cux2","Mef2c","Rorb","Npas4",
                  "Gria1","Nrxn3")


heatmapGS <- markerHeatmap(
  seMarker = ATAC_markersGS, 
  cutOff = "FDR <= 0.05& Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

reorder_row = c('C1','C6',  'C5','C4','C3',
                'C2',  'C8','C9','C10','C12',
                'C13','C18','C16','C17','C15','C14',  'C7','C11','C19') 

ComplexHeatmap::draw(heatmapGS,row_order=reorder_row,heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-ATAC", width = 8, height = 6, ArchRProj = projHeme1, addDOC = FALSE)
#######################################


# RNA
########################################

projHeme1 <- addImputeWeights(projHeme1, reducedDims = "LSI_RNA")


## Identify the marker genes for each cluster 
RNA_markersGS <- getMarkerFeatures(
  ArchRProj = projHeme1, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "Combined_Clusters",
  testMethod = "wilcoxon"
)

dd <- getMarkers(RNA_markersGS, cutOff = "FDR <= 0.05& Log2FC >= 0.5")
dd
markerGenes_dd <- list()
for (i in seq_len(length(dd))) {
  markerGenes_dd <- c(markerGenes_dd, dd[[i]]$name)
}
# 
markerGenes_dd <- unlist(markerGenes_dd)
markerGenes_dd <- unique(markerGenes_dd)
"Klk6" %in% markerGenes_dd




sampleNames <- "LPC_ATAC"


markerList_pos <- getMarkers(RNA_markersGS, cutOff = "FDR <= 0.05& Log2FC >= 0.5")

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}

markerGenes <- unlist(markerGenes)


for (i in seq_len(length(markerList_pos))) {
  write.table(markerList_pos[[i]], file=paste0('./data/LPC_RNA_markers_list/', sampleNames, '_C', i, '_markers.txt'),
              quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
}





heatmapGS <- markerHeatmap(
  seMarker = RNA_markersGS, 
  cutOff = "FDR <= 0.05& Log2FC >= 0.5", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

reorder_row = c('C1','C18','C19','C8','C10','C9','C4','C5','C2','C7','C13', 'C12','C11', 'C16','C17','C14', 'C15',
                'C6',  'C3') 

ComplexHeatmap::draw(heatmapGS,row_order=reorder_row,heatmap_legend_side = "bot", annotation_legend_side = "bot")





plotPDF(heatmapGS, name = "GeneExpressionMatrix-Marker-Heatmap-RNA", width = 8, height = 6, ArchRProj = projHeme1, addDOC = FALSE)
#######################################


###############################################################################
proj_in_tissue

getAvailableMatrices(proj_in_tissue)

seGroupMotif <- getGroupSE(ArchRProj = proj_in_tissue, useMatrix = "MotifMatrix", groupBy = "Combined_Clusters")
seGroupMotif
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
  ArchRProj = proj_in_tissue,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_RNA"
)

corGSM_MM

corGIM_MM <- correlateMatrices(
  ArchRProj = proj_in_tissue,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_RNA"
)

corGIM_MM

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])




data_yes <- subset(data.frame(corGSM_MM), TFRegulator == "YES")

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  geom_text(data = data_yes, aes(label = GeneScoreMatrix_name), vjust = -1) + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p


corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*", "", corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
maxDelta_quantile <- quantile(corGIM_MM$maxDelta, 0.75, na.rm = TRUE)

corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & !is.na(corGIM_MM$cor) & 
                              corGIM_MM$padj < 0.01 & !is.na(corGIM_MM$padj) & 
                              corGIM_MM$maxDelta > maxDelta_quantile & !is.na(corGIM_MM$maxDelta))] <- "YES"

sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p


library(ggplot2)


data_yes <- subset(data.frame(corGIM_MM), TFRegulator == "YES")

p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  geom_text(data = data_yes, aes(label = GeneExpressionMatrix_name), vjust = -1) + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p



##############################################################################
markerGenes <- c("Olig2","Sox10","Neurod6")

markerGenes <- c("Mog","Cx3cr1")

markerGenes <- "Csf1r"
#markerGenes <- c("Lef1","Tcf7l2")

############################################################


proj_in_tissue
getAvailableMatrices(proj_in_tissue)

proj_in_tissue <- addPeak2GeneLinks(
  ArchRProj = proj_in_tissue,
  reducedDims = "LSI_RNA",
  dimsToUse = 1:30,
  useMatrix = "GeneExpressionMatrix"
)



p2g <- getPeak2GeneLinks(
  ArchRProj = proj_in_tissue,
  corCutOff = 0.2,
  FDRCutOff = 1e-04,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = T
)
p2g

getPeak2GeneLinks(proj_in_tissue)




p <- plotBrowserTrack(
  ArchRProj = proj_in_tissue, 
  groupBy = "Combined_Clusters", 
  useGroups = c('C1',  'C2', 'C5',"C12","C14","C15","C16","C17"),
  geneSymbol = markerGenes, 
  upstream = 20000,
  downstream = 20000,
  loops = p2g,
  tileSize = 225,
  pal =bright_colors_25 
)
grid::grid.newpage()


plotPDF(plotList = p, 
        name = "Peak2GeneLinksallpiex_Sox10_V2.pdf", 
        ArchRProj = proj_in_tissue, 
        addDOC = FALSE, width = 5, height = 5)


p <- plotPeak2GeneHeatmap(ArchRProj = proj_in_tissue, 
                          groupBy = "Combined_Clusters" ,
                          corCutOff = 0.2
)
p



plotPDF(p, name = "P2G-Heatmap-Combined_Clusters", ArchRProj = proj_in_tissue)



saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "Save-Peak2GeneLinks", load = FALSE)



proj_in_tissue <- loadArchRProject("Save-Peak2GeneLinks", force = FALSE, showLogo = TRUE)
proj_in_tissue

##############################################################
#######################################################################


proj_in_tissue <- addCoAccessibility(
  ArchRProj = proj_in_tissue,
  dimsToUse = 1:30,
  reducedDims = "LSI_RNA"
)



cA <- getCoAccessibility(
  ArchRProj = proj_in_tissue,
  corCutOff = 0.2,
  resolution = 1,
  returnLoops = T
)

cA

table(proj_in_tissue$Combined_Clusters)

p <- plotBrowserTrack(
  ArchRProj = proj_in_tissue, 
  groupBy = "Combined_Clusters", 
  useGroups = c('C1',  'C2', 'C5',"C12","C14","C15","C16","C17"),
  geneSymbol = markerGenes, 
  upstream = 20000,
  downstream = 20000,
  loops = cA,
  tileSize = 225,
  pal =bright_colors_25  ,
)


grid::grid.newpage()

plotPDF(plotList = p, 
        name = "Peak2peaks1_Combined_Sox10_V2.pdf", 
        ArchRProj = proj_in_tissue, 
        addDOC = FALSE, width = 5, height = 5)


saveArchRProject(ArchRProj = proj_in_tissue, outputDirectory = "Save-ATAC_copeka2peak", load = FALSE)


##############################################################################


proj_in_tissue <-  loadArchRProject("Save-ATAC_copeka2peak", force = FALSE, showLogo = TRUE)
proj_in_tissue


tile_matrix <- getMatrixFromProject(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix",
  verbose = TRUE,
  binarize = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)


saveRDS(tile_matrix, 'LPC_ATAC_tile_matrix.rds')

table(proj_in_tissue$ATAC_Clusters)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_in_tissue,
  useMatrix = "PeakMatrix",
  groupBy = "ATAC_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.1")

peaks_list <- markerList$D16
peaks_bed <- data.frame(seqnames=peaks_list$seqnames, start=peaks_list$start, end=peaks_list$end)
write.table(peaks_bed, paste0('./data/maker_peaks/D16_marker_peaks.bed'), row.names = FALSE, col.names = FALSE, quote=FALSE)




promoterRegion = c(2000, 100)
geneAnnotation = getGeneAnnotation(proj_in_tissue)
genomeAnnotation = getGenomeAnnotation(proj_in_tissue)

ArchR:::.requirePackage(genomeAnnotation$genome)
ArchR:::.requirePackage("Biostrings",source="bioc")
BSgenome <- eval(parse(text = genomeAnnotation$genome))
BSgenome <- validBSgenome(BSgenome)


peaks_file = paste0('./data/maker_peaks/D16_marker_peaks.bed')
peaks_df <- as.data.frame(read.table(peaks_file, header = FALSE, stringsAsFactors=FALSE, quote=""))
colnames(peaks_df) <- c('chr', 'start', 'end')

peaks <- makeGRangesFromDataFrame(peaks_df)
peaks

peaks <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion)
peaks_table <- as.data.frame(prop.table(table(peaks@elementMetadata$peakType)))
peaks_table$Group <- 'D16'
colnames(peaks_table) <- c('PeakAnno', 'Freq', 'Group')


peaks_table_all <- peaks_table


peaks_file = paste0('./data/maker_peaks/D16_marker_peaks.bed')
peaks_df <- as.data.frame(read.table(peaks_file, header = FALSE, stringsAsFactors=FALSE, quote=""))
colnames(peaks_df) <- c('chr', 'start', 'end')

peaks <- makeGRangesFromDataFrame(peaks_df)
peaks

peaks <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion)
peaks_table <- as.data.frame(prop.table(table(peaks@elementMetadata$peakType)))
peaks_table$Group <- 'D16'
colnames(peaks_table) <- c('PeakAnno', 'Freq', 'Group')

peaks_table_all <- rbind(peaks_table_all, peaks_table)

pal <- c(Distal = "#60BA64", Exonic = "#73C6FF", Intronic = "#620FA3", Promoter = "#FFC554")
p <- ggplot(peaks_table_all, aes(x = Group, y = Freq*100, fill = PeakAnno)) + 
  geom_bar(stat = "identity") + theme_ArchR(xText90 = FALSE) + 
  ylab("Percentage (%)") + xlab("") + theme(legend.position = "bottom", 
                                            legend.key = element_rect(size = 2), legend.box.background = element_rect(color = NA)) + 
  scale_fill_manual(values = pal)
p




ggsave("./plot/ME13-C7-C8PeakAnno.png", plot = p, width = 9, height = 9)



