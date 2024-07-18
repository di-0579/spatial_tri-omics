rm(list=ls())
library(Seurat)
library(grid)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(dplyr)


getwd()
setwd("./matlb/")

position <- read.table("./P21_SS.txt", sep = ',', stringsAsFactors = FALSE)
position <- position[-1]
position <- as.data.frame(t(position))
position$in_tisue <- 1
row.names(position) <- position$V1


BC <- read.table("./uesful/barcode_100.txt", sep = '\t', stringsAsFactors = FALSE)
BC$position <- paste0(BC$V2, 'x', BC$V3)
row.names(BC) <- BC$position


df_in_tissue <- merge( position, BC, by=0, all=TRUE)
df_in_tissue[is.na(df_in_tissue)] <- 0
df_in_tissue$V2 <- df_in_tissue$V2 - 1
df_in_tissue$V3 <- df_in_tissue$V3 - 1


w = 3154 # picture width pixel number
l = 3145 # picture length pixel number
w1 = w/199
w2 = w/398
l1 = l/199
l2 = l/398


#For 1, 3 or 5, run this line:
df_in_tissue$width_loc = round((df_in_tissue$V3)*w1*2 + w2)
#For 2 or 4, run this line:
df_in_tissue$width_reversed = round(w - ((df_in_tissue$V3)*w1*2 + w2))

df_in_tissue$high_loc = round ((df_in_tissue$V2)*l1*2 + l2)


df_in_tissue <- df_in_tissue[c("V1.y", "in_tisue", "V3", "V2", "width_loc", "high_loc")]

write.table(df_in_tissue, "tissue_positions_list.csv", sep=",", row.names=FALSE, col.names = FALSE, quote = FALSE)



min_value <- min(w2, l2)
result <- min_value * 4
result_rounded <- round(result, 1)
result_rounded

scale_factors <- list(
  spot_diameter_fullres = result_rounded,
  tissue_hires_scalef = 1.0,
  fiducial_diameter_fullres = result_rounded,
  tissue_lowres_scalef = 1.0
)

json_str <- toJSON(scale_factors, auto_unbox = TRUE, pretty = TRUE)
write(json_str, "scalefactors_json.json")

