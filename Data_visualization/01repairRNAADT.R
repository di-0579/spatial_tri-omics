rm(list=ls())
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(grid)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
set.seed(1234)


source("./00repair.R")
getwd()
setwd("/Q43_P5/")


sampleNames <- "Q43_P5"



counts <- Read10X_h5("./Processed_data/P5_Q43_ARC/outs/raw_feature_bc_matrix.h5")

mat_tmp <- as.matrix(counts$`Gene Expression`)
mat_tmp <- as.data.frame(mat_tmp)


spatial_data <- mat_tmp
new_col_names <- colnames(spatial_data)
new_col_names <- unlist(lapply(new_col_names, function(x) gsub("-.*", "", x)))
colnames(spatial_data) <- new_col_names


data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = spatial_data, assay = assay,min.cells = 10)#, meta.data = meta.data)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object



plot2 <- SpatialFeaturePlot(spatial.obj, features = "nCount_Spatial",pt.size.factor = 1,max.cutoff = 800) + theme(legend.position = "right")
plot2

plot2 <- SpatialFeaturePlot(spatial.obj, features = "nCount_Spatial",pt.size.factor = 1,max.cutoff = 5000) + theme(legend.position = "right")
plot2



coord_tmp <- read.csv("./spatial/tissue_positions_list.csv",header = FALSE)
coord_tmp$V1 <- paste(coord_tmp$V1, '-1', sep = '')


names(coord_tmp) <- c("BC", "in_tissue", "y", "x", "width_reversed", "high_loc")
rownames(coord_tmp) <- coord_tmp$BC



mat_tmp <- mat_tmp[, rownames(coord_tmp)]


col_to_fix <- 74
mat_fixed <- imputeCol(mat_tmp = mat_tmp, coord_tmp = coord_tmp, col_to_fix = col_to_fix, direction = "both")




spatial_data <- mat_fixed
new_col_names <- colnames(spatial_data)
new_col_names <- unlist(lapply(new_col_names, function(x) gsub("-.*", "", x)))
colnames(spatial_data) <- new_col_names




data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = spatial_data, assay = assay,min.cells = 10)#, meta.data = meta.data)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object



plot2 <- SpatialFeaturePlot(spatial.obj, features = "nCount_Spatial",pt.size.factor = 1,max.cutoff = 800) + theme(legend.position = "right")
plot2

plot2 <- SpatialFeaturePlot(spatial.obj, features = "nCount_Spatial",pt.size.factor = 1,max.cutoff = 5000) + theme(legend.position = "right")
plot2

spatial.obj <- SCTransform(spatial.obj, assay = "Spatial", verbose = FALSE)
spatial.obj <- RunPCA(spatial.obj, assay = "SCT", verbose = FALSE)
spatial.obj <- FindNeighbors(spatial.obj, reduction = "pca", dims = 1:30)
spatial.obj <- FindClusters(spatial.obj, verbose = FALSE)
spatial.obj <- RunUMAP(spatial.obj, reduction = "pca", dims = 1:30)
p2 <- SpatialDimPlot(spatial.obj , label = TRUE, pt.size.factor = 1,label.size = 3)
p2




spatial_tissue <- read.csv("./spatial/tissue_positions_list.csv",header = FALSE)
spatial_tissue <- spatial_tissue[spatial_tissue$V2 != 0, ]

BC_tissue <- spatial_tissue$V1
spatial.tissue <- subset(spatial.obj, cells = BC_tissue)

if (!dir.exists("object")) {
  dir.create("object")
}



saveRDS(spatial.tissue, file = "./object/00sptialRepairRNA_brain.rds")


################################################################
#################################################################

extracted_data <- spatial_data[, BC_tissue]
write.csv(extracted_data, file =paste0("./data/", sampleNames, "_RepairRNA_matrix.csv"), row.names = TRUE, col.names = TRUE)


library(DropletUtils)
library(Matrix)

if (is.null(rownames(mat_fixed3))) {
  stop("The matrix must have rownames representing gene IDs.")
}
if (is.null(colnames(mat_fixed3))) {
  stop("The matrix must have colnames representing barcodes.")
}

mat_fixed3_sparse <- as(Matrix(as.matrix(mat_fixed3), sparse = TRUE), "dgCMatrix")

write10xCounts(
  path = "Repair_RNA_feature_bc_matrix.h5",
  x = mat_fixed3_sparse,
  barcodes = colnames(mat_fixed),
  gene.id = rownames(mat_fixed),
  gene.symbol = rownames(mat_fixed),  
  overwrite = TRUE,
  type = "HDF5",
  genome = "GRCm38",  
  gene.type = "Gene Expression",
  version = "3",
  chemistry = "Single Cell 3' v3"
)


##################################################################################################################
################################################################################################

data.dir <-"./Processed_data/cite/P5_Q43_citeouts/"
ADT_data <- Seurat::Read10X(data.dir = paste(data.dir, "/umi_count/",sep =""),gene.column = 1)
my_df <-as.data.frame(t(as.matrix(ADT_data)))

data_filtered <-as.data.frame(t(my_df))
new_row_names <- row.names(data_filtered)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-[^-]*$","", x)))
row.names(data_filtered) <- new_row_names


#data_filtered <- data_filtered[rownames(data_filtered) != "unmapped",]


mat_tmp <- data_filtered

coord_tmp <- read.csv("./spatial/tissue_positions_list.csv",header = FALSE)


names(coord_tmp) <- c("BC", "in_tissue", "y", "x", "width_reversed", "high_loc")
rownames(coord_tmp) <- coord_tmp$BC



mat_tmp <- mat_tmp[, rownames(coord_tmp)]



col_to_fix <- 74
mat_fixed <- imputeCol(mat_tmp = mat_tmp, coord_tmp = coord_tmp, col_to_fix = col_to_fix, direction = "both")




mat_fixed3 <- mat_fixed[, rownames(coord_tmp)]
row_to_fix <- 7
mat_fixed3 <- imputeRow(mat_tmp = mat_fixed3, coord_tmp = coord_tmp, row_to_fix = row_to_fix, direction = "both")



spatial_data <-mat_fixed3




data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = spatial_data, assay = assay,min.cells = 10)#, meta.data = meta.data)
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object



plot2 <- SpatialFeaturePlot(spatial.obj, features = "nCount_Spatial",pt.size.factor = 1,max.cutoff = 800) + theme(legend.position = "right")
plot2


spatial_tissue <- read.csv("./spatial/tissue_positions_list.csv",header = FALSE)
spatial_tissue <- spatial_tissue[spatial_tissue$V2 != 0, ]

BC_tissue <- spatial_tissue$V1

#################################################################

extracted_data <- spatial_data[, BC_tissue]
write.csv(extracted_data, file =paste0("./data/", sampleNames, "_RepairADT_intissue_matrix.csv"), row.names = TRUE, col.names = TRUE)
#################################################################

spatial_tissue <- read.csv("./spatial/tissue_positions_list.csv",header = FALSE)
spatial_tissue <- spatial_tissue[spatial_tissue$V2 != 1, ]

BC_tissue <- spatial_tissue$V1

#################################################################

extracted_data <- spatial_data[, BC_tissue]
write.csv(extracted_data, file =paste0("./data/", sampleNames, "_RepairADT_not_intissue_matrix.csv"), row.names = TRUE, col.names = TRUE)


