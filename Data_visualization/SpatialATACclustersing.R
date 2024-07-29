rm(list=ls())
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
#library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(EnsDb.Mmusculus.v79)



getwd()
setwd("./sginacmerger_11smaple/")




sample_names <- c("P0S1", "P2S1", "P5S1", "P10S1", "P21S1",
           "P0S2", "P0S3", "P2S2", "P5S2", "P10S2", "P21S2")

peaks <- list()

for (sample in sample_names) {
  file_path <- sprintf("./mergeATAC_11smaple/%s/atac_peaks.bed", sample)
  peaks[[sample]] <- read.table(
    file = file_path,
    col.names = c("chr", "start", "end")
  )
}




# Convert to genomic ranges
genomic_ranges <- lapply(peaks[c("P0S1", "P2S1", "P5S1", "P10S1", "P21S1",
           "P0S2", "P0S3", "P2S2", "P5S2", "P10S2", "P21S2")], 
                         makeGRangesFromDataFrame)

# Combine peaks and reduce
combined.peaks <- reduce(do.call(c, genomic_ranges))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

library(GenomicRanges)

combined.peaks <- combined.peaks[grepl("^chr", seqnames(combined.peaks))]
combined.peaks
################################################################################

md_list <- list()
images_list <- list()

for (sample in sample_names) {

  md_path <- sprintf("./mergeATAC_11smaple/%s/per_barcode_metrics.csv", sample)
  md_list[[sample]] <- read.table(
    file = md_path,
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ]  # remove the first row
  
  image_dir <- sprintf("./mergeATAC_11smaple/%s/spatial/", sample)
  images_list[[sample]] <- Read10X_Image(
    image.dir = image_dir,
    filter.matrix = TRUE
  )
  
  bc_spatial <- row.names(images_list[[sample]]@coordinates)
  new_bc <- paste(bc_spatial, '-1', sep = '')
  md_list[[sample]] <- md_list[[sample]][row.names(md_list[[sample]]) %in% new_bc, ]
}

##################################################################################################
fragments_list <- list()

for (sample in sample_names) {
  fragment_path <- sprintf("./mergeATAC_11smaple/%s/atac_fragments.tsv.gz", sample)
  
  fragments_list[[sample]] <- CreateFragmentObject(
    path = fragment_path,
    cells = rownames(md_list[[sample]])
  )
}

#################################################################################################
counts_list <- list()

for (sample in sample_names) {
  counts_list[[sample]] <- FeatureMatrix(
    fragments = fragments_list[[sample]],
    features = combined.peaks,
    cells = rownames(md_list[[sample]])
  )
}


#############################
assay_list <- list()
seurat_objects <- list()

for (sample in sample_names) {
  assay_list[[sample]] <- CreateChromatinAssay(
    counts = counts_list[[sample]],
    fragments = fragments_list[[sample]]
  )
  
  seurat_objects[[sample]] <- CreateSeuratObject(
    assay_list[[sample]],
    assay = "ATAC",
    meta.data = md_list[[sample]]
  )
}


#####################################
for (sample in sample_names) {
  seurat_objects[[sample]]$dataset <- sample
}


for (sample in names(seurat_objects)) {
  age <- sub(".*P", "P", sample)
  seurat_objects[[sample]]$Age <- age
}


############################################
combined <- merge(
  x = seurat_objects[["P0S1"]],
  y = list(
    seurat_objects[["P2S1"]],
    seurat_objects[["P5S1"]],
    seurat_objects[["P10S1"]],
    seurat_objects[["P21S1"]],
    seurat_objects[["P0S2"]],
    seurat_objects[["P0S3"]],
    seurat_objects[["P2S2"]],
    seurat_objects[["P5S2"]],
    seurat_objects[["P10S2"]],
    seurat_objects[["P21S2"]]
  ),
  add.cell.ids = c("P0S1", "P2S1", "P5S1", "P10S1", "P21S1",
           "P0S2", "P0S3", "P2S2", "P5S2", "P10S2", "P21S2")
)


combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined,n=30)
DepthCor(combined)

combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')

col_method <- c("P0"="#FDE725FF","P2"="#AADC32FF","P5"="#27AD81FF","P10"="#3B528BFF","P21"="#440154FF")


col_method <- c(
  "P22"="#FDE725FF", 
  "P21"="#AADC32FF", 
  "P10"="#6FCF97FF", 
  "P5"="#1A9875FF",  
  "P2"="#3B528BFF", 
  "P0"="#440154FF"
)



p1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1,cols = col_method)
p1

brain.list <- SplitObject(combined, split.by = "dataset")
brain.list
brain.list <- lapply(X = brain.list, FUN = function(x) {
  x <- RunTFIDF(x)
  x <- FindTopFeatures(x, min.cutoff = 'q5')
  x <- RunSVD(x)
})

length(rownames(P5))

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = brain.list,
  #  reference = c(1, 2),
  anchor.features = rownames(P10),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p22 <- DimPlot(integrated, group.by = "dataset",pt.size = 0.1,cols = col_method)

(p1 + ggtitle("Merged")) | (p22 + ggtitle("Integrated"))


integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = FALSE,resolution = 1.5, algorithm = 3)




bright_colors_32 <- c( "0"= "#66c2a5", "1"= "#fc8d62", "2"="#8da0cb", "3"="#e78ac3","4"="#a6d854",
                       "5"="#ffd92f", "6"="#e5c494", "7"="#b3b3b3", "8"= "#8dd3c7","9"= "#bebada", 
                       "10"= "#fb8072", "11"="#80b1d3", "13"= "#fccde5","14"= "#bc80bd","12"="#000000" ,#"12"="#fdb462",
                       "15"= "#ffed6f", "16"= "#f768a1", "17"= "#7fc97f", "18"= "#beaed4","19"= "#fdc086", 
                       "20"="#386cb0", "21"="#bf5b17", "22"="#f0027f", "23"="#66a61e", "24"="#f781bf", 
                       "25"= "#e41a1c", "26"="#377eb8","27"= "#4daf4a", "28"="#92c5de", "29"="#d95f02", "30"="#7570b3", "31"="#e7298a"
)



table(integrated$seurat_clusters)

bright_colors_13 <- c(
  "0"= "#e41a1c",  
  "1"= "#377eb8",  
  "2"= "#4daf4a",  
  "3"= "#984ea3",  
  "4"= "#ff7f00", 
  "5"= "#ffff33",  
  "6"= "#a65628",  
  "7"= "#f781bf", 
  "8"= "#999999",  
  "9"= "#17becf",  
  "10"= "#b15928", 
  "11"= "#6a3d9a", 
  "12"= "#33a02c" ,"14"="#f0027f",
  "13"="#000000"
)

p2 <- DimPlot(object = integrated,label = FALSE,cols = bright_colors_32) + NoLegend()
p2
p <- p22|p2
p
ggsave("plot/integratedUMAP.png", plot = p, width = 18, height = 9)

#saveRDS(integrated, file = "./object/integrated1013.rds")

integrated<- readRDS("./object/integrated1013-11sample.rds")
#####################################################

library(paletteer)#
col <- paletteer_d("ggthemes::Classic_20")#[c(6,9,3,7,18,8,5,1,          
#  13,2,17,15,10,11,12,16)]
col
col <- c("#1F77B4FF","#AEC7E8FF","#FFBB78FF","#9EDAE5FF","#2CA02CFF",
         "#98DF8AFF","#FF7F0EFF","#FF9896FF","#D62728FF")


integrated.matadata <- integrated@meta.data


table(integrated$dataset)
table(integrated$seurat_clusters)
table(integrated$ident_merged)
integrated$ident_merged <- integrated$seurat_clusters

bright_colors_32 <- c( "0"= "#66c2a5", "1"= "#fc8d62", "2"="#8da0cb", "3"="#e78ac3","4"="#a6d854",
                       "5"="#ffd92f", "6"="#e5c494", "7"="#b3b3b3", "8"= "#8dd3c7","9"= "#bebada", 
                       "10"= "#fb8072", "11"="#80b1d3", "13"= "#fccde5","14"= "#bc80bd","12"="#000000" ,#"12"="#fdb462",
                       "15"= "#ffed6f", "16"= "#AEC7E8FF", "17"= "#AEC7E8FF", "18"= "#AEC7E8FF","19"= "#AEC7E8FF", 
                       "20"="#AEC7E8FF", "21"="#AEC7E8FF", "22"="#f0027f", "23"="#66a61e", "24"="#f781bf", 
                       "25"= "#e41a1c", "26"="#377eb8","27"= "#4daf4a", "28"="#92c5de", "29"="#d95f02", "30"="#7570b3", "31"="#e7298a"
)


Idents(integrated) <- "ident_merged"
table(integrated$ident_merged)
table(integrated$Age)

p2 <- DimPlot(object = integrated,label = FALSE , group.by = "ident_merged",  cols = bright_colors_32) 
p2



ggsave(paste0("plot/", sampleNames, "_11RNAUMP_lab.png"), plot = p2, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_11RNAUMP.png"), plot = p2, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_11RNAUMP.pdf"), plot = p2, width = 9, height = 9)






p3 <- DimPlot(object = integrated,label = FALSE , group.by = "Age",  cols = col_method) 
p3

sampleNames <- "11sampleATAC"
ggsave(paste0("plot/", sampleNames, "_AGE_ATACUMP.png"), plot = p3, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_AGE_ATACUMP.pdf"), plot = p3, width = 9, height = 9)




p4 <- DimPlot(object = integrated,label = FALSE , group.by = "dataset") 
p4

ggsave(paste0("plot/", sampleNames, "_sampleNames_ATACUMP.png"), plot = p4, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_sampleNames_ATACUMP.pdf"), plot = p4, width = 9, height = 9)





datasets <- sample_names
datasets

for (dataset in datasets) {
  p <- DimPlot(integrated[, integrated$dataset == dataset], label = TRUE, reduction = "umap", cols = bright_colors_32) + ggtitle(dataset)
  file_name <- sprintf("plot/integrated%sUMAP.png", dataset)
  ggsave(file_name, plot = p, width = 9, height = 9)
}


datasets 

for (dataset in datasets) {
  seurat_objects[[dataset]]$ident_merged <- integrated[, integrated$dataset == dataset]@active.ident
  
  saveRDS(seurat_objects[[dataset]], file = paste0("./object/integrated_", dataset, ".rds"))
  
  df <- seurat_objects[[dataset]]@meta.data
  write.csv(df, file = paste0("./data/integrated_", dataset, ".csv"), row.names = TRUE, col.names = TRUE)
}


#######################################################################
seurat_objects[["Q43P0"]]


####################################################################
library(grid)
source("/gpfs/gibbs/pi/fan/dz286/ArchRTutorial/code/SpatialPlot_new.R")


spatial.P0 <- readRDS(file = "./object/Q43_P0_spatialfilterADT.rds")
spatial.P0
df <- seurat_objects[["Q43P0"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P0 <- AddMetaData(spatial.P0, metadata = df)




bright_colors_32 <- c( "120"= "#66c2a5","0"="#a6d854", "14"= "#5cf0fe","1"= "#ff9500", "9"="#00ffff","2"="#0091ff",#"#8da0cb", 
                       "6"="#fff900", "17"= "#5cf0fe", "10"= "#00ff00","7"=  "#ff05e1", "5"="#ffd92f", 
                       "12"= "#fb8b72", "11"="#ff0000", "33"= "#fccde5","3"= "#bc80bd","8"="#000000" ,#"12"="#fdb462",
                       "15"= "#ffed6f",  "17"= "#AEC7E8FF", "18"= "#AEC7E8FF","19"= "#AEC7E8FF", "4" = "#313cff",
                       "13"="#AEC7E8FF", "21"="#AEC7E8FF", "62"="#f0027f", "23"="#66a61e", "24"="#f781bf", 
                       "11111"= "#e41a1c", "26"="#377eb8","422"= "#4daf4a", "28"="#92c5de", "29"="#d95f02", "30"="#7570b3", "31"="#e7298a"
)





bright_colors_25 <- c("C1" = "#7570b3", "C2" = "#ff836c",  "C3" = "#009e00",  "C4" = "#ffe400", "C5" = "#99cada", 
                      "C6" = "#0091ff",  "C7" = "#ff0000",  "C8" = "#ff05e1","C9" = "#cff9fa",  "C10" = "#000000", 
                      "C11" = "#8eff05", "C12" = "#5cf0fe",    "C14" = "#e7298a",  "C15" = "#20b2aa", 
                      "C16" = "#91eda4",  "C25" = "#00ff00",  "C17" = "#00ffff",   "C18" = "#f9e5d4",  "C19" = "#ff9500", 
                      "C20" = "#0059ca", "C21" = "#f868a1",  "C22" = "#b3b3b3",  "C23" = "#ffd130",  "C24" = "#689b94" 
)


Idents(spatial.P0) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P0,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 10)+ ggtitle('Q43P0')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 



ggsave("plot/spatialUMAP/integratedQ43P0spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedQ43P0spatialUMAP.pdf", plot = p51, width = 9, height = 9)





spatial.P2 <- readRDS(file = "./object/Q43_P2_spatialfilterUMPADT.rds")

spatial.P2
df <- seurat_objects[["Q43P2"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P2 <- AddMetaData(spatial.P2, metadata = df)

Idents(spatial.P2) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P2,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 3)+ ggtitle('Q43P2')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 




ggsave("plot/spatialUMAP/integratedQ43P2spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedQ43P2spatialUMAP.pdf", plot = p51, width = 9, height = 9)








spatial.P5 <- readRDS(file = "./object/Q43_P5_spatialfilterADT.rds")
spatial.P5
df <- seurat_objects[["Q43P5"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P5 <- AddMetaData(spatial.P5, metadata = df)

Idents(spatial.P5) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P5,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 3)+ ggtitle('Q43P5')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 


ggsave("plot/spatialUMAP/integratedQ43P5spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedQ43P5spatialUMAP.pdf", plot = p51, width = 9, height = 9)





spatial.P10 <- readRDS(file = "./Q43e_P10/object/Q43_P10_spatialfilterADT.rds")
spatial.P10
df <- seurat_objects[["Q43P10"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P10 <- AddMetaData(spatial.P10, metadata = df)

Idents(spatial.P10) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P10,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 6)+ ggtitle('Q43P10')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 





ggsave("plot/spatialUMAP/integratedQ43P10spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedQ43P10spatialUMAP.pdf", plot = p51, width = 9, height = 9)





spatial.P21 <- readRDS(file = "./object/00sptialRepairRNA_brain.rds")
spatial.P21
df <- seurat_objects[["Q43P21"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P21 <- AddMetaData(spatial.P21, metadata = df)

table(spatial.P21$ident_merged)
spatial.P21$ident_merged[spatial.P21$ident_merged == "14"] <- "4"

Idents(spatial.P21) <- "ident_merged"


p51 <- SpatialPlot_new(spatial.P21,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE, cols = bright_colors_32, image.alpha = 0,
                       label.size = 7)+ ggtitle('Q43P21')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 




ggsave("plot/spatialUMAP/integratedQ43P21spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedQ43P21spatialUMAP.pdf", plot = p51, width = 9, height = 9)





p2 <- DimPlot(object = integrated,label = TRUE , group.by = "ident_merged",  cols = bright_colors_32) 
p2



ggsave(paste0("plot/", sampleNames, "_11ATACUMP_lab.png"), plot = p2, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_11ATACUMP.png"), plot = p2, width = 9, height = 9)
ggsave(paste0("plot/", sampleNames, "_11ATACUMP.pdf"), plot = p2, width = 9, height = 9)


########################################################################################

spatial.P0 <- readRDS(file = "./mergrRNA_11smaple/object/C8P0mergedspatialRNAUMP.rds")
spatial.P0
df <- seurat_objects[["C8P0"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P0 <- AddMetaData(spatial.P0, metadata = df)


Idents(spatial.P0) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P0,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 10)+ ggtitle('C8P0')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 



ggsave("plot/spatialUMAP/integratedC8P0spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedC8P0spatialUMAP.pdf", plot = p51, width = 9, height = 9)




spatial.P0 <- readRDS(file = "./mergrRNA_11smaple/object/B8P0mergedspatialRNAUMP.rds")
spatial.P0
df <- seurat_objects[["B8P0"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P0 <- AddMetaData(spatial.P0, metadata = df)


Idents(spatial.P0) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P0,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 10)+ ggtitle('B8P0')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 



ggsave("plot/spatialUMAP/integratedB8P0spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedB8P0spatialUMAP.pdf", plot = p51, width = 9, height = 9)








spatial.P2 <- readRDS(file ="./mergrRNA_11smaple/object/B8P2mergedspatialRNAUMP.rds")

spatial.P2
df <- seurat_objects[["B8P2"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P2 <- AddMetaData(spatial.P2, metadata = df)

Idents(spatial.P2) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P2,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 3)+ ggtitle('B8P2')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 




ggsave("plot/spatialUMAP/integratedB8P2spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedB8P2spatialUMAP.pdf", plot = p51, width = 9, height = 9)








spatial.P5 <- readRDS(file = "./mergrRNA_11smaple/object/C8P5mergedspatialRNAUMP.rds")
spatial.P5
df <- seurat_objects[["C8P5"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P5 <- AddMetaData(spatial.P5, metadata = df)

Idents(spatial.P5) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P5,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 3)+ ggtitle('C8P5')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 


ggsave("plot/spatialUMAP/integratedC8P5spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedC8P5spatialUMAP.pdf", plot = p51, width = 9, height = 9)





spatial.P10 <- readRDS(file = "./mergrRNA_11smaple/object/C8P10mergedspatialRNAUMP.rds")
spatial.P10
df <- seurat_objects[["C8P10"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P10 <- AddMetaData(spatial.P10, metadata = df)

Idents(spatial.P10) <- "ident_merged"
p51 <- SpatialPlot_new(spatial.P10,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE,cols = bright_colors_32, image.alpha = 0,
                       label.size = 6)+ ggtitle('C8P10')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 





ggsave("plot/spatialUMAP/integratedC8P10spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedC8P10spatialUMAP.pdf", plot = p51, width = 9, height = 9)





spatial.P21 <- readRDS(file = "./mergrRNA_11smaple/object/C8P22mergedspatialRNAUMP.rds")
spatial.P21
df <- seurat_objects[["C8P22"]]@meta.data
table(df$dataset)
new_row_names <- row.names(df)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(df) <- new_row_names
df <- df["ident_merged"]
spatial.P21 <- AddMetaData(spatial.P21, metadata = df)

table(spatial.P21$ident_merged)
spatial.P21$ident_merged[spatial.P21$ident_merged == "13"] <- "12"

Idents(spatial.P21) <- "ident_merged"


p51 <- SpatialPlot_new(spatial.P21,group.by = 'ident_merged', pt.size.factor = 1, label = FALSE, cols = bright_colors_32, image.alpha = 0,
                       label.size = 7)+ ggtitle('C8P22')+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p51$layers[[1]]$aes_params <- c(p51$layers[[1]]$aes_params, shape=22) 
p51 




ggsave("plot/spatialUMAP/integratedC8P22spatialUMAP.png", plot = p51, width = 9, height = 9)
ggsave("plot/spatialUMAP/integratedC8P22spatialUMAP.pdf", plot = p51, width = 9, height = 9)






###########################################################################################
df <- integrated@meta.data

write.csv(df, "./data/ATAC_11sample_QC_ATACcluster.csv", row.names = T, col.names = T, quote=F)
