rm(list=ls())
library(SeuratObject)
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
setwd("/Matrix_data/LPC_data/ATAC/object/")
#############
############################################################
source("/ArchRTutorial/code/SpatialPlot_new.R")
source("/ArchRTutorial/code/getGeneScore_ArchR.R")



library(ggplot2)
library(scales)
library(SingleCellExperiment)

# Read data
spatial.P0 <- readRDS("5DPL_spatialRNA_geneExp_V2.rds")
spatial.P2 <- readRDS("21DPL_spatialRNA_geneExp_V2.rds")





bright_colors_25 <- c( "C1" = "#ff0000", "C2" = "#FFF200", "C3" = "#f9e5d4", "C4" = "#20b2aa", "C5" = "#00ff00",
                       "C6" = "#91eda4", "C7" = "#7570b3", "C8" = "#0059ca", "C9" = "#ff8000", "C10" = "#313cff",
                       "C11" = "#8eff05", "C12" = "#5cf0fe", "C13" = "#cff9fa", "C14" = "#ffa6e3", "C15" = "#009e00",
                       "C16" = "#0091ff", "C17" = "#ff05e1", "C18" = "#f868a1", "C19" = "#9500ff"
)

p7 <- SpatialPlot_new(spatial.P0, label = F, label.size = 5, group.by = 'Combined_Clusters',
                      pt.size.factor = 1.2, cols = bright_colors_25, image.alpha = 1, stroke = 0)+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p7$layers[[1]]$aes_params <- c(p7$layers[[1]]$aes_params, shape=22)
p7



p7 <- SpatialPlot_new(spatial.P2, label = F, label.size = 5, group.by = 'Combined_Clusters',
                      pt.size.factor = 1.2, cols = bright_colors_25, image.alpha = 1, stroke = 0)+
  theme_void() +xlab("") + ylab("") +theme(axis.text = element_blank(), axis.ticks = element_blank())
p7$layers[[1]]$aes_params <- c(p7$layers[[1]]$aes_params, shape=22)
p7






# Merge data
spatial.RNAmerge <- merge(spatial.P0, spatial.P2)


# Define features


feature_names <- c("Ank","Anxa5","Aplp2","Hpse","Igf1","Itgax",
                   "Pirb", "Rplp1","Nceh1","Plaur","Pld3","Plin2","Spp1",
                   "Atp1a3","Ephx1","Fabp5","Fam20c",
                   "Colec12","Gm1673",
                   "Csf1","Gpnmb")


feature_names <- c("Itgax","Tmem119","Igf1","Spp1","Itgam")


feature_names <- "Cd44"

# Define custom scale function
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
  lower_bound <- quantile(all_expression_values, 0.0001, na.rm = TRUE)
  upper_bound <- quantile(all_expression_values, 0.9999, na.rm = TRUE)
  
  # Specify gradient colors (this needs to be defined or adjusted as per your preference)
  gradient_colors <- c("#313695", "#0571B0","#91CF60", "#FFFFBF", "#FE9929", "#D73027",  "#8E0152" )
  
  # Plot for each time point
  time_points <- list(spatial.P0, spatial.P2)
  time_labels <- c("5DPL", "21DPL")
  
  for (i in 1:length(time_points)) {
    p <- SpatialFeaturePlot(time_points[[i]], features = feature, pt.size.factor = 1.2, image.alpha = 0, stroke = 0, alpha = c(1, 1)) +
      theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15)) +
      theme_void() + xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
    p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) 
    
    p1 <- p + custom_scale(lower_bound, upper_bound, gradient_colors)
    
    # Save plot
    file_path_png <- paste0("plot_Fig5/", feature, "/", time_labels[i], "_RNA.png")
    ggsave(file_path_png, plot = p1, width = 9, height = 9)
  }
}
