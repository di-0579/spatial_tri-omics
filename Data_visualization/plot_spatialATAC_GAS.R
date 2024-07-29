rm(list=ls())
library(SeuratObject)
library(Seurat)

library(ArchR)
library(grid)
library(scales)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggplot2)
library(patchwork)
library(future)
set.seed(123)

getwd()
setwd("./Matrix_data/LPC_data/ATAC/object/")
#############
############################################################
source("./ArchRTutorial/code/SpatialPlot_new.R")
source("./ArchRTutorial/code/getGeneScore_ArchR.R")


gradient_colors <- c("#0641e7","#0f93ff","#0bd1ff","#fed30b","#fc0808", "#d00405")





spatial.P0 <- readRDS("5DPL_spatialATAC_gene_V2.rds")

spatial.P2<- readRDS("21DPL_spatialATAC_gene_V2.rds")


spatial.RNAmerge <- merge(spatial.P0,spatial.P2)


spatial.P2 <- ScaleData(spatial.P2, features = rownames(spatial.P2))
spatial.P0 <- ScaleData(spatial.P0, features = rownames(spatial.P0))




feature_names <- c("Itgax","Itgam")


custom_scale <- function(lower_bound, upper_bound, gradient_colors) {
  scale_fill_gradientn(
    colors = gradient_colors,
    limits = c(lower_bound, upper_bound),
    oob = scales::squish
  )
}



# Loop over each feature
for (feature in feature_names) {
  
  all_expression_values <- FetchData(spatial.RNAmerge, vars = feature)[[feature]]
  gene_exprs <- lapply(list(spatial.P0, spatial.P2), function(x) x@assays$Spatial@scale.data[feature, ])
  lapply(gene_exprs, summary)
  lower_bound <- mean(unlist(lapply(gene_exprs, min)))
  upper_bound <- mean(unlist(lapply(gene_exprs, max)))
  
  lower_bound <- -mean(c(abs(lower_bound), abs(upper_bound)))
  upper_bound <- mean(c(abs(lower_bound), abs(upper_bound)))
  
  gradient_colors <- c("#0641e7","#0f93ff","#0bd1ff","#fed30b","#fc0808", "#d00405" )
  
  # Plot for each time point
  time_points <- list(spatial.P0, spatial.P2)
  time_labels <- c("LPC21_M39_2", "LPC5_N40_4")
  
  for (i in 1:length(time_points)) {
    p <- SpatialFeaturePlot(time_points[[i]], features = feature, pt.size.factor = 1.2, image.alpha = 0,
                            stroke = 0, alpha = c(1, 1), slot = "scale.data") +
      theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) +
      theme_void() + xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
    p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) 
    
    p1 <- p + custom_scale(lower_bound, upper_bound, gradient_colors)
    
    # Save plot
    file_path_png <- paste0("plot_Fig5/", feature, "/", time_labels[i], "_ATAC.png")
    ggsave(file_path_png, plot = p1, width = 9, height = 9)
  }
}
