rm(list=ls())

setwd("")

#####################################################################
#####################    Settings    ################################
#####################################################################

set.seed(2025)  
library(mgcv)
library(Seurat)
library(chromVAR)
library(doParallel)
library(FigR)
library(FNN)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)
library(ArchR)
library(ggnewscale)
library(RColorBrewer)
library(ggpubr)
addArchRGenome("mm10")
ncore <- 28
addArchRThreads(threads = ncore)


pign <- colorRampPalette(c(rev(brewer.pal(n = 7, name = "RdBu"))[1:3],
                           rev(brewer.pal(n = 7, name = "PRGn"))[c(4:7)]))(100)
brbg <- colorRampPalette(c(rev(brewer.pal(n = 7, name = "PRGn"))[1:3],
                           rev(brewer.pal(n = 7, name = "BrBG"))[c(4:7)]))(100)
pugn <- colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(100)

rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

theme_yx <- function() {
  theme(
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 10, color = "black"),
    axis.text.y  = element_text(size = 10, color = "black")
  )
}

ggplot2::theme_set(theme_bw() + theme_yx() + 
                     theme(axis.text.y = element_text(color = "black"),
                           axis.text.x = element_text(color = "black")) )

input_dir <- "/cortex/"
figures_dir <- "/cortex/figures/"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

#####################################################################
##################### Fit GAM model for RNA-seq #####################
#####################################################################

refined_corpus_region_idx <- read.csv(file.path(input_dir, "cortex_region_idx.csv"),
                                      row.names = 1)




seu_obj_corpus <- readRDS(file.path(input_dir, "cortex_seu_byCluster.rds"))
library_id <- colSums(seu_obj_corpus)

head(library_id)
head(rownames(refined_corpus_region_idx))
library_id <- library_id[rownames(refined_corpus_region_idx)]  #



gene_list <- rownames(seu_obj_corpus)

data_dir <- file.path(figures_dir, "data")
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

write.csv(
  gene_list,
  file = file.path(data_dir, "background_genes.csv"),
  row.names = FALSE
)


table(refined_corpus_region_idx$CellType)
table(seu_obj_corpus$Age)
##############################################################

exprs_mat <- seu_obj_corpus@assays$SCT@counts

exprs_mat <- exprs_mat[, rownames(refined_corpus_region_idx)]


head(colnames(exprs_mat))

sample <- factor(refined_corpus_region_idx$Sample)

spatial <- factor(refined_corpus_region_idx$CellType,
                  levels = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"))
spatial <- as.numeric(spatial)


time <- as.numeric(gsub("P", "", refined_corpus_region_idx$Age))
names(time) <- rownames(refined_corpus_region_idx)


keep_genes <- names(which(rowMeans(exprs_mat != 0) > 0.02))
library(pbmcapply)

gam_fit <- pbmcapply::pbmclapply(keep_genes, function(gene) {
  y_exprs <- exprs_mat[gene, ]
  gam(y_exprs ~ library_id + ti(time, k = 5) + ti(spatial, k = 4) + ti(time, spatial, k = 4), 
      family = nb())
}, mc.cores = 30)


names(gam_fit) <- keep_genes

gam_fit_library_id <- pbmcapply::pbmclapply(keep_genes, function(gene) {
  y_exprs <- exprs_mat[gene, ]
  gam(y_exprs ~ library_id, 
      family = nb())
}, mc.cores = 30)

names(gam_fit_library_id) <- keep_genes


#####################################################################
##################### Fit GAM model for ATAC-seq ####################
#####################################################################

dorcMat.s <- readRDS("refined_cortex_region_DORC_mat.rds")

dorcMat.s <- dorcMat.s[, rownames(refined_corpus_region_idx)]

all(colnames(dorcMat.s) == rownames(refined_corpus_region_idx))

stopifnot(setequal(colnames(dorcMat.s), rownames(refined_corpus_region_idx)))

sample <- factor(refined_corpus_region_idx[colnames(dorcMat.s), ]$Sample)


spatial <- factor(refined_corpus_region_idx[colnames(dorcMat.s), "CellType"],
                  levels = c("ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"))
spatial <- as.numeric(spatial)

time <- as.numeric(gsub("P", "", seu_obj_corpus[, colnames(dorcMat.s)]$Age))


signac_obj_cortex <- readRDS("./cortex/cortex_signac_byCluster.rds")

exprs_mat <- signac_obj_cortex@assays$ATAC@counts

head(colnames(exprs_mat))
head(colnames(dorcMat.s))

library_id <- colSums(exprs_mat)

names(library_id) <- sub("_", "#", names(library_id))

head(names(library_id))
stopifnot(all(colnames(dorcMat.s) %in% names(library_id)))
library_id <- library_id[colnames(dorcMat.s)]  

stopifnot(length(time) == ncol(dorcMat.s))
stopifnot(length(spatial) == ncol(dorcMat.s))
stopifnot(length(library_id) == ncol(dorcMat.s))


library(pbmcapply)
library(mgcv)

all_genes <- rownames(dorcMat.s)
remaining_genes <- all_genes
gam_fit_atac_final <- list()

iteration <- 1
max_iter <- 10  

while (length(remaining_genes) > 0 && iteration <= max_iter) {
  message(sprintf("Fitting iteration %d: %d genes", iteration, length(remaining_genes)))
  
  fit_results <- pbmcapply::pbmclapply(remaining_genes, function(gene) {
    y_exprs <- log2(dorcMat.s[gene, ] + 1)
    data <- data.frame(y_exprs, time, spatial, library_id)
    tryCatch(
      gam(y_exprs ~ ti(time, k = 5) + ti(spatial, k = 4) + ti(time, spatial, k = 4), data = data),
      error = function(e) NULL
    )
  }, mc.cores = 30)
  names(fit_results) <- remaining_genes
  
  success <- fit_results[!sapply(fit_results, is.null)]
  fail <- fit_results[sapply(fit_results, is.null)]
  
  gam_fit_atac_final <- c(gam_fit_atac_final, success)
  
  remaining_genes <- names(fail)
  iteration <- iteration + 1
}

if (length(remaining_genes) > 0) {
  warning(sprintf("Still %d genes failed after %d iterations", length(remaining_genes), iteration - 1))
}

saveRDS(gam_fit_atac, file = file.path(figures_dir, "cortex_region_DORC_GAM_fit_P21.rds"))




gam_fit_atac_library_id <- pbmcapply::pbmclapply(1:nrow(dorcMat.s), function(gene) {
  y_exprs <- log2(dorcMat.s[gene, ] + 1)
  data <- data.frame(y_exprs, time, spatial, library_id)
  gamma_model <- gam(y_exprs ~ library_id,
                     data = data)
}, mc.cores = 30)



names(gam_fit_atac_library_id) <- rownames(dorcMat.s)


####################################
#####################################################################
##################### Analysis of GAM results - RNA #################
#####################################################################

length(gam_fit) == length(gam_fit_library_id)
all(names(gam_fit) == names(gam_fit_library_id))


length(gam_fit_atac) == length(gam_fit_atac_library_id)
all(names(gam_fit_atac) == names(gam_fit_atac_library_id))

head(names(gam_fit))
head(names(gam_fit_atac))


all(colnames(exprs_mat) == colnames(dorcMat.s))  


lrt_test <- lapply(names(gam_fit), function(x) {
  anova(gam_fit[[x]], gam_fit_library_id[[x]], test = "LRT")
})
lrt_test_pvalue <- lapply(lrt_test, function(x) x$`Pr(>Chi)`[2])
lrt_test_pvalue <- unlist(lrt_test_pvalue)
lrt_test_pvalue_adj <- p.adjust(lrt_test_pvalue, method = "BH")
names(lrt_test_pvalue_adj) <- names(gam_fit)


sum_s_table <- lapply(gam_fit, function(x) {
  summary(x)$s.table
})
sum_res <- summary(gam_fit[[2]])$s.table

adjr_squred <- lapply(gam_fit, function(x) {
  summary(x)$r.sq
})
adjr_squred <- unlist(adjr_squred)
summary(adjr_squred)

sum_s_table_time <- do.call(rbind, lapply(sum_s_table, function(x) x[1, ]))
sum_s_table_time <- data.frame(sum_s_table_time)
sum_s_table_time$adj.p.value <- p.adjust(sum_s_table_time$p.value, method = "BH")
sum(sum_s_table_time$adj.p.value < 0.001)
sum_s_table_spatial <- do.call(rbind, lapply(sum_s_table, function(x) x[2, ]))
sum_s_table_spatial <- data.frame(sum_s_table_spatial)
sum_s_table_spatial$adj.p.value <- p.adjust(sum_s_table_spatial$p.value, method = "BH")
sum(sum_s_table_spatial$adj.p.value < 0.001)
sum_s_table_spacetime <- do.call(rbind, lapply(sum_s_table, function(x) x[3, ]))
sum_s_table_spacetime <- data.frame(sum_s_table_spacetime)
sum_s_table_spacetime$adj.p.value <- p.adjust(sum_s_table_spacetime$p.value, method = "BH")
sum(sum_s_table_spacetime$adj.p.value < 0.001)


df_genes <- data.frame(temporal = sum_s_table_time$adj.p.value < 0.01 , 
                       spatial = sum_s_table_spatial$adj.p.value < 0.01 , 
                       spatiotemporal = sum_s_table_spacetime$adj.p.value < 0.01)
df_genes$category <- paste(df_genes$temporal, df_genes$spatial, df_genes$spatiotemporal, sep = "_")
rownames(df_genes) <- rownames(sum_s_table_time)
gene_list <- split(rownames(df_genes), df_genes$category)
table(df_genes$category)

selected_genes <- df_genes$category != "FALSE_FALSE_FALSE" & 
  rownames(df_genes) %in% names(which(adjr_squred > 0.025)) & 
  rownames(df_genes) %in% names(which(lrt_test_pvalue_adj < 0.001))
selected_genes <- rownames(df_genes)[selected_genes]
length(selected_genes)
grep("Rpl|Rps|^mt-", selected_genes, value = TRUE)

selected_genes <- selected_genes[!grepl("^mt-", selected_genes)]
length(selected_genes)

predicted_gam_mat <- lapply(gam_fit, function(x) predict(x, type = "response"))
predicted_gam_mat <- do.call(rbind, predicted_gam_mat)
predicted_gam_mat <- log2(predicted_gam_mat + 1)

all(colnames(predicted_gam_mat) == rownames(refined_corpus_region_idx))
setequal(colnames(predicted_gam_mat), rownames(refined_corpus_region_idx))
predicted_gam_mat <- predicted_gam_mat[, rownames(refined_corpus_region_idx)]
dim(predicted_gam_mat)

####################################################################
##################### Analysis of GAM results - ATAC ################
#####################################################################

lrt_test_atac <- lapply(names(gam_fit_atac), function(x) {
  anova(gam_fit_atac[[x]], gam_fit_atac_library_id[[x]], test = "LRT")
})
lrt_test_atac_pvalue <- lapply(lrt_test_atac, function(x) x$`Pr(>Chi)`[2])
lrt_test_atac_pvalue <- unlist(lrt_test_atac_pvalue)
lrt_test_atac_pvalue_adj <- p.adjust(lrt_test_atac_pvalue, method = "BH")
names(lrt_test_atac_pvalue_adj) <- names(gam_fit_atac)
sum(lrt_test_atac_pvalue_adj <= 0.001, na.rm = TRUE)
length( grep("Rpl|Rps|^mt-", names(which(lrt_test_atac_pvalue_adj <= 0.001)), value = TRUE))


sum_s_table_atac <- lapply(gam_fit_atac, function(x) {
  summary(x)$s.table
})
sum_s_table_time_atac <- do.call(rbind, lapply(sum_s_table_atac, function(x) x[1, ]))
sum_s_table_time_atac <- data.frame(sum_s_table_time_atac)
sum_s_table_time_atac$adj.p.value <- p.adjust(sum_s_table_time_atac$p.value, method = "BH")
sum(sum_s_table_time_atac$adj.p.value < 0.001)
sum_s_table_spatial_atac <- do.call(rbind, lapply(sum_s_table_atac, function(x) x[2, ]))
sum_s_table_spatial_atac <- data.frame(sum_s_table_spatial_atac)
sum_s_table_spatial_atac$adj.p.value <- p.adjust(sum_s_table_spatial_atac$p.value, method = "BH")
sum(sum_s_table_spatial_atac$adj.p.value < 0.001)
sum_s_table_spacetime_atac <- do.call(rbind, lapply(sum_s_table_atac, function(x) x[3, ]))
sum_s_table_spacetime_atac <- data.frame(sum_s_table_spacetime_atac)
sum_s_table_spacetime_atac$adj.p.value <- p.adjust(sum_s_table_spacetime_atac$p.value, method = "BH")
sum(sum_s_table_spacetime_atac$adj.p.value < 0.001)


adjr_squred_atac <- lapply(gam_fit_atac, function(x) {
  summary(x)$r.sq
})
adjr_squred_atac <- unlist(adjr_squred_atac)


df_atac <- data.frame(temporal = sum_s_table_time_atac$adj.p.value < 0.01 , 
                      spatial = sum_s_table_spatial_atac$adj.p.value < 0.01 , 
                      spatiotemporal = sum_s_table_spacetime_atac$adj.p.value < 0.01)
df_atac$category <- paste(df_atac$temporal, df_atac$spatial, df_atac$spatiotemporal, sep = "_")
rownames(df_atac) <- rownames(sum_s_table_time_atac)
gene_list <- split(rownames(df_atac), df_atac$category)
table(df_atac$category)


selected_peaks <- df_atac$category != "FALSE_FALSE_FALSE" & 
  rownames(df_atac) %in% names(which(adjr_squred_atac > 0.025)) &
  lrt_test_atac_pvalue_adj < 0.001 &
  !is.na(lrt_test_atac_pvalue_adj)
selected_peaks <- rownames(df_atac)[selected_peaks]
length(selected_peaks)


predicted_gam_atac_mat <- lapply(gam_fit_atac, function(x) predict(x, type = "response"))
predicted_gam_atac_mat <- do.call(rbind, predicted_gam_atac_mat)

all(colnames(predicted_gam_atac_mat) == rownames(refined_corpus_region_idx))
setequal(colnames(predicted_gam_atac_mat), rownames(refined_corpus_region_idx))

predicted_gam_atac_mat <- predicted_gam_atac_mat[, rownames(refined_corpus_region_idx)]

#####################################################################
##################### Joint clustering  #############################
####################################################################

common_genes <- intersect(intersect(names(gam_fit), names(gam_fit_atac)),
                          union(selected_genes, selected_peaks))

design_grid <- expand.grid(
  Group.1 = c("P0", "P2", "P5", "P7", "P10", "P21"),
  Group.2 = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
)

agg_values <- sapply(1:nrow(design_grid), function(i) {
  age      <- design_grid$Group.1[i]
  celltype <- design_grid$Group.2[i]
  idx <- which(refined_corpus_region_idx$Age == age & refined_corpus_region_idx$CellType == celltype)
  mean(predicted_gam_mat[1, idx], na.rm = TRUE)
})

agg_df <- cbind(design_grid, value = agg_values)

agg_df_combine <- rbind(
  design_grid,  
  design_grid  
)
agg_df_combine$modality <- c(
  rep("RNA", nrow(design_grid)),
  rep("ATAC", nrow(design_grid))
)

head(agg_df_combine)
table(design_grid$Age)

getAggMat <- function(predicted_gam_mat, time_point, spatial_bin) {
  design_grid <- expand.grid(
    Age = c("P0", "P2", "P5","P7", "P10", "P21"),
    CellType = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
  )
  design_levels <- paste(design_grid$Age, design_grid$CellType, sep = "_")
  
  predicted_mat_agg <- apply(predicted_gam_mat, 1, function(expr_vec) {
    agg_result <- aggregate(expr_vec, 
                            by = list(Age = time_point, CellType = spatial_bin), 
                            FUN = mean)
    keys <- paste(agg_result$Age, agg_result$CellType, sep = "_")
    stats::setNames(agg_result$x, keys)
  })
  
  predicted_mat_agg <- t(predicted_mat_agg)
  
  missing_cols <- setdiff(design_levels, colnames(predicted_mat_agg))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      predicted_mat_agg <- cbind(predicted_mat_agg, setNames(rep(NA, nrow(predicted_mat_agg)), col))
    }
  }
  predicted_mat_agg <- predicted_mat_agg[, design_levels]  
  
  gene_means <- rowMeans(predicted_mat_agg, na.rm = TRUE)
  gene_sds   <- apply(predicted_mat_agg, 1, sd, na.rm = TRUE)
  
  predicted_agg_mat_scale <- sweep(predicted_mat_agg, 1, gene_means, "-")
  predicted_agg_mat_scale <- sweep(predicted_agg_mat_scale, 1, gene_sds, "/")
  
  return(list(
    agg   = predicted_mat_agg,
    scale = predicted_agg_mat_scale
  ))
}

length(common_genes)

agg_rna <- getAggMat(predicted_gam_mat[common_genes, ], 
                     refined_corpus_region_idx$Age,     
                     refined_corpus_region_idx$CellType)  

predicted_gam_rna_mat_agg       <- agg_rna$agg
predicted_gam_rna_mat_agg_scale <- agg_rna$scale

agg_atac <- getAggMat(predicted_gam_atac_mat[common_genes, ], 
                      refined_corpus_region_idx$Age, 
                      refined_corpus_region_idx$CellType)
predicted_gam_atac_mat_agg       <- agg_atac$agg
predicted_gam_atac_mat_agg_scale <- agg_atac$scale

predicted_gam_combine_mat_agg_scale <- cbind(
  predicted_gam_rna_mat_agg_scale, 
  predicted_gam_atac_mat_agg_scale
)


#######################################
hclust_res <- hclust(
  1 - proxy::simil(predicted_gam_combine_mat_agg_scale),  
  method = "ward.D2"                                      
)

combine_hclust_cluster <- cutree(hclust_res, k = 30)

table(combine_hclust_cluster)
sample_row <- data.frame(cluster = factor(combine_hclust_cluster))
rownames(sample_row) <- names(combine_hclust_cluster)    

cluster_ids <- unique(combine_hclust_cluster)

rna_cluster_avg <- sapply(cluster_ids, function(cid) {
  colMeans(predicted_gam_rna_mat_agg_scale[names(combine_hclust_cluster)[combine_hclust_cluster == cid], ])
})

atac_cluster_avg <- sapply(cluster_ids, function(cid) {
  colMeans(predicted_gam_atac_mat_agg_scale[names(combine_hclust_cluster)[combine_hclust_cluster == cid], ])
})

rna_cluster_avg <- t(rna_cluster_avg)
atac_cluster_avg <- t(atac_cluster_avg)

h_rna <- hclust(dist(rna_cluster_avg), method = "ward.D2")
RNA_big_cluster <- cutree(h_rna, k = 15)
names(RNA_big_cluster) <- cluster_ids
h_atac <- hclust(dist(atac_cluster_avg), method = "ward.D2")
ATAC_big_cluster <- cutree(h_atac, k = 10)
names(ATAC_big_cluster) <- cluster_ids
combine_hclust_cluster_big <- paste0(
  "R", RNA_big_cluster[as.character(combine_hclust_cluster)],
  "_A", ATAC_big_cluster[as.character(combine_hclust_cluster)]
)

length(unique(combine_hclust_cluster_big))
table(combine_hclust_cluster_big)
names(combine_hclust_cluster) <- rownames(predicted_gam_combine_mat_agg_scale)
names(combine_hclust_cluster_big) <- rownames(predicted_gam_combine_mat_agg_scale)
head(names(combine_hclust_cluster) <- rownames(predicted_gam_combine_mat_agg_scale))
length(combine_hclust_cluster_big)
nrow(predicted_gam_combine_mat_agg_scale)
head(names(combine_hclust_cluster_big))
stopifnot(length(combine_hclust_cluster_big) == nrow(predicted_gam_combine_mat_agg_scale))
gene_order <- match(names(combine_hclust_cluster), rownames(predicted_gam_combine_mat_agg_scale))
combine_hclust_cluster_big <- combine_hclust_cluster_big[gene_order]
names(combine_hclust_cluster_big) <- rownames(predicted_gam_combine_mat_agg_scale)


#####################################################################
##################### Visualization   . #############################

predicted_gam_combine_mat_agg_scale_mean_big <- sapply(
  names(table(combine_hclust_cluster_big)),      
  function(x) {
    colMeans(predicted_gam_combine_mat_agg_scale[combine_hclust_cluster_big == x, ])
  }
)

agg_df_mat <- reshape2::melt(
  data.frame(
    cbind(agg_df_combine, predicted_gam_combine_mat_agg_scale_mean_big)
  ),
  id.vars = c("Group.1", "Group.2", "modality")  
)

agg_df_mat$Group.1 <- factor(
  agg_df_mat$Group.1,
  levels = c("P0", "P2", "P5","P7" ,"P10", "P21")
)

agg_df_mat$Group.2 <- factor(
  agg_df_mat$Group.2,
  levels = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
)

agg_df_mat$modality <- factor(agg_df_mat$modality, levels = c("RNA", "ATAC"))

ggplot(agg_df_mat, 
       aes(
         x = Group.1, 
         y = Group.2, 
         fill = value                      
       )) +
  geom_tile(aes(fill = value),
            data = agg_df_mat[agg_df_mat$modality == "ATAC", ]) +
  scale_fill_gradientn(
    colours = pign,
    limits = c(-2, 2),
    oob = scales::oob_squish,
    name = "ATAC"
  ) +
  new_scale_fill() +
  geom_tile(aes(fill = value),
            data = agg_df_mat[agg_df_mat$modality == "RNA", ]) +
  scale_fill_gradientn(
    colours = rdbu,
    limits = c(-2, 2),
    oob = scales::oob_squish,
    name = "RNA"
  ) +
  theme(
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    legend.position  = "right",
    strip.background = element_blank()
  ) +
  facet_grid(modality ~ variable) +
  theme(aspect.ratio = 1) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")



meta_df <- refined_corpus_region_idx %>%
  mutate(spot = rownames(refined_corpus_region_idx)) %>%
  dplyr::select(spot, x, y, Age, CellType)

colnames(meta_df)[colnames(meta_df) == "Age"] <- "Group.1"
colnames(meta_df)[colnames(meta_df) == "CellType"] <- "Group.2"

meta_rna  <- meta_df %>% mutate(modality = "RNA")
meta_atac <- meta_df %>% mutate(modality = "ATAC")
meta_full <- rbind(meta_rna, meta_atac)
agg_df_mat$cluster <- agg_df_mat$variable

agg_spot <- left_join(meta_full, agg_df_mat, 
                      by = c("Group.1", "Group.2", "modality"))

nrow(meta_full)        
length(unique(agg_df_mat$cluster))  

nrow(agg_spot)

agg_spot$Group.1  <- factor(agg_spot$Group.1, levels = c("P0", "P2", "P5","P7", "P10", "P21"))
agg_spot$Group.2  <- factor(agg_spot$Group.2, levels = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"))
agg_spot$modality <- factor(agg_spot$modality, levels = c("RNA", "ATAC"))
agg_spot$cluster  <- factor(agg_spot$cluster, levels = unique(agg_df_mat$variable))  

ggplot(agg_spot, aes(x = x, y = y)) +
  geom_tile(data = subset(agg_spot, modality == "ATAC"), aes(fill = value)) +
  scale_fill_gradientn(
    colours = pign, limits = c(-2, 2),
    oob = scales::oob_squish, name = "ATAC"
  ) +
  new_scale_fill() +
  geom_tile(data = subset(agg_spot, modality == "RNA"), aes(fill = value)) +
  scale_fill_gradientn(
    colours = rdbu, limits = c(-2, 2),
    oob = scales::oob_squish, name = "RNA"
  ) +
  facet_grid(cluster + modality ~ Group.1) +
  theme_minimal() +
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    strip.background = element_blank(),
    aspect.ratio     = 0.6
  ) +
  labs(x = NULL, y = NULL, fill = "Scaled expression")




##########################
# —— RNA pattern heatmap —— 
predicted_gam_combine_mat_agg_scale_mean_big <- sapply(
  paste0("R", 1:15),  
  function(x) { 
    colMeans(predicted_gam_combine_mat_agg_scale[grepl(x, combine_hclust_cluster_big), ])
  }
)

agg_df_mat <- reshape2::melt(
  data.frame(cbind(agg_df_combine, predicted_gam_combine_mat_agg_scale_mean_big)),
  id.vars = c("Group.1", "Group.2", "modality")
)

agg_df_mat$Group.1 <- factor(
  agg_df_mat$Group.1, 
  levels = c("P21", "P10","P7", "P5", "P2", "P0")
)

ggplot(
  agg_df_mat[agg_df_mat$modality == "RNA", ], 
  aes(
    x = factor(Group.1, levels = rev(levels(agg_df_mat$Group.1))),  
    y = factor(Group.2, levels = c(                               
      "ExciteL6",
      "ExciteL5",
      "ExciteL4",
      "ExciteL23"
    )), 
    fill = value
  )
) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(
    colours = rdbu, 
    limits = c(-2, 2), 
    oob = scales::oob_squish, 
    name = "RNA"
  ) +
  theme(
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    legend.position  = "right",
    strip.background = element_blank()
  ) +
  facet_grid(~variable) +
  theme(aspect.ratio = 1.4) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")

 
predicted_gam_combine_mat_agg_scale_mean_big <- sapply(
  paste0("A", 1:10),
  function(x) { 
    colMeans(predicted_gam_combine_mat_agg_scale[grepl(x, combine_hclust_cluster_big), ])
  }
)

agg_df_mat <- reshape2::melt(
  data.frame(cbind(agg_df_combine, predicted_gam_combine_mat_agg_scale_mean_big)),
  id.vars = c("Group.1", "Group.2", "modality")
)

agg_df_mat$Group.1 <- factor(
  agg_df_mat$Group.1, 
  levels = c("P21", "P10","P7", "P5", "P2", "P0")
)


ggplot(
  agg_df_mat[agg_df_mat$modality == "ATAC", ], 
  aes(
    x = factor(Group.1, levels = rev(levels(agg_df_mat$Group.1))),
    y = factor(Group.2, levels = c(
      "ExciteL6",
      "ExciteL5",
      "ExciteL4",
      "ExciteL23"
    )),
    fill = value
  )
) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(
    colours = pign, 
    limits = c(-2, 2), 
    oob = scales::oob_squish, 
    name = "ATAC"
  ) +
  theme(
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    legend.position  = "right",
    strip.background = element_blank()
  ) +
  facet_grid(~variable) +
  theme(aspect.ratio = 1.4) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")


#################################
#################################
predicted_gam_combine_mat_agg_scale_mean_big_atac <- 
  predicted_gam_combine_mat_agg_scale_mean_big[
    agg_df_combine$modality == "ATAC", 
  ]

agg_df_mat_atac <- reshape2::melt(
  data.frame(
    Age     = agg_df_combine$Group.1[agg_df_combine$modality=="ATAC"], 
    CellType= agg_df_combine$Group.2[agg_df_combine$modality=="ATAC"],
    predicted_gam_combine_mat_agg_scale_mean_big_atac
  ),
  id.vars = c("Age", "CellType")        
)

agg_df_mat_atac$modality <- "ATAC"


typeof(combine_hclust_cluster_big)     
class(combine_hclust_cluster_big)      

combine_hclust_cluster_big_char <- as.character(combine_hclust_cluster_big)
combine_hclust_cluster_big_rna <- sapply(strsplit(combine_hclust_cluster_big_char, "_"), `[`, 1)


combine_hclust_cluster_big_rna <- 
  unlist(lapply(strsplit(combine_hclust_cluster_big, "_"), "[[", 1))

predicted_gam_combine_mat_agg_scale_mean_big_rna <- sapply(
  names(table(combine_hclust_cluster_big_rna)), 
  function(x) {
    colMeans(
      predicted_gam_combine_mat_agg_scale[
        combine_hclust_cluster_big_rna == x, 
      ]
    )
  }
)
predicted_gam_combine_mat_agg_scale_mean_big_rna <- 
  predicted_gam_combine_mat_agg_scale_mean_big_rna[
    agg_df_combine$modality=="RNA", 
  ]

agg_df_mat_rna <- reshape2::melt(
  data.frame(
    Age     = agg_df_combine$Group.1[agg_df_combine$modality=="RNA"], 
    CellType= agg_df_combine$Group.2[agg_df_combine$modality=="RNA"],
    predicted_gam_combine_mat_agg_scale_mean_big_rna
  ),
  id.vars = c("Age", "CellType")
)
agg_df_mat_rna$modality <- "RNA"

agg_df_mat_combine <- rbind(agg_df_mat_rna, agg_df_mat_atac)

colnames(agg_df_mat_combine)[1:4] <- c("Age", "CellType", "cluster", "value", "modality")[1:4]


spatial_key <- rbind(refined_corpus_region_idx, refined_corpus_region_idx)

agg_df_mat_spatial <- merge(
  spatial_key,
  agg_df_mat_combine,
  by = c("Age", "CellType")
)

agg_df_mat_spatial$Sample <- factor(agg_df_mat_spatial$Sample, levels = c("P0",
                                                                          "P2",  "P5","P7",
                                                                          "P10", "P21"))


table(agg_df_mat_spatial$Sample)
agg_df_mat_spatial$Age <- factor(
  agg_df_mat_spatial$Age,
  levels = c("P0","P2","P5","P7",  "P10","P21")
)

agg_df_mat_spatial$CellType <- factor(
  agg_df_mat_spatial$CellType,
  levels = c("ExciteL6","ExciteL5","ExciteL4","ExciteL23")
)

ggplot(agg_df_mat_spatial, aes(x = x, y = y)) +
  geom_tile(
    data = subset(agg_df_mat_spatial, modality == "ATAC"),
    aes(fill = value)
  ) +
  scale_fill_gradientn(
    colours = pign,
    limits  = c(-2, 2),
    oob     = scales::oob_squish,
    name    = "ATAC"
  ) +
  new_scale_fill() +
  geom_tile(
    data = subset(agg_df_mat_spatial, modality == "RNA"),
    aes(fill = value)
  ) +
  scale_fill_gradientn(
    colours = rdbu,
    limits  = c(-2, 2),
    oob     = scales::oob_squish,
    name    = "RNA"
  ) +
  facet_grid(cluster + modality ~ Age) +
  theme_minimal() +
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    strip.background = element_blank(),
    aspect.ratio     = 0.6
  ) +
  labs(x = NULL, y = NULL, fill = "Z‑score")


###########


seu_obj_cortex <- readRDS(file.path(input_dir, "cortex_seu_byCluster.rds"))

dorcMat.s <- readRDS("cortex_region_DORC_mat.rds")


for (gene in gene_list_filtered) {
  message("Plotting gene: ", gene)
  
  dir.create(paste0(figures_dir, "gene_figures_cortex/"), recursive = TRUE, showWarnings = FALSE)
  
  p <- tryCatch(plotEachGene(gene), error = function(e) {
    message("plotEachGene failed for ", gene, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(p)) {
    tryCatch({
      ggsave(paste0(figures_dir, "gene_figures_cortex/", gene, ".pdf"),
             plot = p, width = 16, height = 12)
    }, error = function(e) {
      message("ggsave failed for ", gene, ": ", conditionMessage(e))
    })
  }
  
  rm(p)
  gc()
}


plotHeatmapEachGene <- function(gene) {
  split_info <- strsplit(colnames(predicted_gam_rna_mat_agg_scale), "_")
  meta_df <- do.call(rbind, lapply(split_info, function(x) data.frame(Age = x[1], CellType = x[2])))
  meta_df$Age <- as.character(meta_df$Age)
  meta_df$CellType <- as.character(meta_df$CellType)
  
  Age_levels <- c("P0", "P2", "P5","P7", "P10", "P21")
  CellType_levels <- c("ExciteL23", "ExciteL4", "ExciteL5","ExciteL6" )
  
  df_rna <- meta_df
  df_rna$value <- predicted_gam_rna_mat_agg_scale[gene, ]
  df_rna$Age <- factor(df_rna$Age, levels = Age_levels)
  df_rna$CellType <- factor(df_rna$CellType, levels = rev(CellType_levels))
  
  g_rna <- ggplot(df_rna, aes(x = Age, y = CellType, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = rdbu, limits = c(-2, 2), oob = scales::oob_squish) +
    labs(title = paste0(gene, " — RNA (scaled prediction)"), fill = "Z-score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1)
  
  df_atac <- meta_df
  df_atac$value <- predicted_gam_atac_mat_agg_scale[gene, ]
  df_atac$Age <- factor(df_atac$Age, levels = Age_levels)
  df_atac$CellType <- factor(df_atac$CellType, levels = rev(CellType_levels))
  
  g_atac <- ggplot(df_atac, aes(x = Age, y = CellType, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = pign, limits = c(-2, 2), oob = scales::oob_squish) +
    labs(title = paste0(gene, " — ATAC (scaled prediction)"), fill = "Z-score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          aspect.ratio = 1)
  
  cowplot::plot_grid(g_rna, g_atac, ncol = 2)
}


plotHeatmapEachGene("Mbp")

plotHeatmapGeneList <- function(gene_list, save_dir = "./gene_heatmaps/") {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (gene in gene_list) {
    message("Plotting: ", gene)
    p <- plotHeatmapEachGene(gene)  
    
    ggsave(
      filename = file.path(save_dir, paste0("heatmap_", gene, ".pdf")),
      plot     = p,
      width    = 10,
      height   = 5
    )
  }
}

gene_list <- c("Bcl11b", "Tbr1", "Satb2", "Cux1", "Foxp1", "Lhx2", "Mef2c", "Pou3f2",
               "Bcl11a", "Etv1", "Meis2", "Myt1l", "Neurod2", "Neurod6", "Pbx1", "Satb2", "Sox4", "Sox11",
               "Cux2", "Neurod6", "Fezf2", "Tle4", "Ldb2", "Sox5", "Tbr1", "Ndn", "Tcf4", "Hpcal1",
               "Zfp428", "Chgb", "Plxna4", "Satb2", "Nfe2l3", "Nfib", "Rbfox3", "Tle4",
               "Vcan", "Spon1", "Pcdh15", "2610035D17Rik", "Ptprz1", "Serpine2", "Grin3a",
               "Sox6", "Fyn", "Fam107b", "Frmd4a", "Bcas1", "Sept8", "Nrxn3", "Nrgn", "S100b",
               "Grin2d", "Grin3a", "Grin2a",
               "Egr1", "Mobp", "Mag", "Mbp")

dir.create(file.path(figures_dir, "cortex_heatmaps"), recursive = TRUE, showWarnings = FALSE)

plotHeatmapGeneList(gene_list, save_dir = file.path(figures_dir, "cortex_heatmaps"))

###########################

plotHeatmapGeneSet <- function(
    gene_set,
    rna_mat, atac_mat,
    celltype_levels = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"),
    age_levels = c("P0", "P2", "P5","P7", "P10", "P21"),
    zlim = c(-2, 2)
) {
  stopifnot(all(gene_set %in% rownames(rna_mat)))
  stopifnot(identical(colnames(rna_mat), colnames(atac_mat)))
  
  col_ids <- colnames(rna_mat)
  split_col <- strsplit(col_ids, "_")
  col_meta <- do.call(rbind, lapply(split_col, function(x) {
    data.frame(Age = x[1], CellType = x[2])
  }))
  rownames(col_meta) <- col_ids
  col_meta$Age      <- factor(col_meta$Age, levels = age_levels)
  col_meta$CellType <- factor(col_meta$CellType, levels = rev(celltype_levels))  
  
  df_rna <- reshape2::melt(as.matrix(rna_mat[gene_set, ]))
  colnames(df_rna) <- c("gene", "col_id", "value")
  df_rna <- cbind(df_rna, col_meta[df_rna$col_id, ])
  df_rna$modality <- "RNA"
  
  df_atac <- reshape2::melt(as.matrix(atac_mat[gene_set, ]))
  colnames(df_atac) <- c("gene", "col_id", "value")
  df_atac <- cbind(df_atac, col_meta[df_atac$col_id, ])
  df_atac$modality <- "ATAC"
  
  df_plot <- rbind(df_rna, df_atac)
  df_plot$gene     <- factor(df_plot$gene, levels = gene_set)
  df_plot$Age      <- factor(df_plot$Age, levels = age_levels)
  df_plot$CellType <- factor(df_plot$CellType, levels = rev(celltype_levels))
  
  p <- ggplot(df_plot, aes(x = gene, y = CellType, fill = value)) +
    geom_tile(data = subset(df_plot, modality == "ATAC"), aes(fill = value)) +
    scale_fill_gradientn(
      colours = pign,
      limits  = zlim,
      oob     = scales::oob_squish,
      name    = "ATAC"
    ) +
    new_scale_fill() +
    geom_tile(data = subset(df_plot, modality == "RNA"), aes(fill = value)) +
    scale_fill_gradientn(
      colours = rdbu,
      limits  = zlim,
      oob     = scales::oob_squish,
      name    = "RNA"
    ) +
    facet_grid(modality ~ Age) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid  = element_blank(),
      strip.background = element_blank(),
      axis.title = element_blank(),
      aspect.ratio = 0.6
    ) +
    labs(fill = "Z-score")
  
  return(p)
}






plotHeatmapGeneSet <- function(
    gene_set,
    rna_mat, atac_mat,
    celltype_levels = c("ExciteL23", "ExciteL4", "ExciteL5","ExciteL6"  ),
    age_levels = c("P0", "P2", "P5","P7", "P10", "P21"),
    zlim = c(-2, 2)
) {
  stopifnot(all(gene_set %in% rownames(rna_mat)))
  stopifnot(identical(colnames(rna_mat), colnames(atac_mat)))
  
  col_ids <- colnames(rna_mat)
  split_col <- strsplit(col_ids, "_")
  col_meta <- do.call(rbind, lapply(split_col, function(x) {
    data.frame(Age = x[1], CellType = x[2])
  }))
  rownames(col_meta) <- col_ids
  col_meta$Age      <- factor(col_meta$Age, levels = age_levels)
  col_meta$CellType <- factor(col_meta$CellType, levels = rev(celltype_levels))  
  
  df_rna <- reshape2::melt(as.matrix(rna_mat[gene_set, ]))
  colnames(df_rna) <- c("gene", "col_id", "value")
  df_rna <- cbind(df_rna, col_meta[df_rna$col_id, ])
  df_rna$modality <- "RNA"
  
  df_atac <- reshape2::melt(as.matrix(atac_mat[gene_set, ]))
  colnames(df_atac) <- c("gene", "col_id", "value")
  df_atac <- cbind(df_atac, col_meta[df_atac$col_id, ])
  df_atac$modality <- "ATAC"
  
  df_plot <- rbind(df_rna, df_atac)
  df_plot$gene     <- factor(df_plot$gene, levels = gene_set)
  df_plot$Age      <- factor(df_plot$Age, levels = age_levels)
  df_plot$CellType <- factor(df_plot$CellType, levels = rev(celltype_levels))
  
  p <- ggplot(df_plot, aes(x = gene, y = CellType, fill = value)) +
    geom_tile(data = subset(df_plot, modality == "ATAC"), aes(fill = value)) +
    scale_fill_gradientn(
      colours = pign,
      limits  = zlim,
      oob     = scales::oob_squish,
      name    = "ATAC"
    ) +
    new_scale_fill() +
    geom_tile(data = subset(df_plot, modality == "RNA"), aes(fill = value)) +
    scale_fill_gradientn(
      colours = rdbu,
      limits  = zlim,
      oob     = scales::oob_squish,
      name    = "RNA"
    ) +
    facet_grid(modality ~ Age) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid  = element_blank(),
      strip.background = element_blank(),
      axis.title = element_blank(),
      aspect.ratio = 0.6
    ) +
    labs(fill = "Z-score")
  
  return(p)
}




gene_set <- c("Rbfox3","Sox11","Sox4","Pbx1","Neurod6","Zeb2",
              "Tbr1","Ndn","Tle4","Zfp428","Chgb","Tcf4",
              "Bcl11b","Crym","Ldb2","Fezf2","Etv1","Thy1",
              "Cux1","Cux2","Lhx2","Dok5","Plxna4","Ptprk")


plotHeatmapGeneSet(
  gene_set,
  rna_mat  = predicted_gam_rna_mat_agg_scale,
  atac_mat = predicted_gam_atac_mat_agg_scale
)




Age_levels <- c("P0", "P2", "P5", "P7",  "P10", "P21")
CellType_levels <- c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")

get_expression_trend_df <- function(gene, predicted_rna, predicted_atac, meta) {
  meta$Age <- factor(meta$Age, levels = Age_levels)
  meta$CellType <- factor(meta$CellType, levels = CellType_levels)
  
  df_rna <- meta %>%
    mutate(value = predicted_rna[gene, rownames(meta)],
           modality = "RNA") %>%
    group_by(Age, CellType, modality) %>%
    summarise(expr = mean(value), .groups = "drop")
  
  df_atac <- meta %>%
    mutate(value = predicted_atac[gene, rownames(meta)],
           modality = "ATAC") %>%
    group_by(Age, CellType, modality) %>%
    summarise(expr = mean(value), .groups = "drop")
  
  df_plot <- bind_rows(df_rna, df_atac)
  df_plot$modality <- factor(df_plot$modality, levels = c("RNA", "ATAC"))
  return(df_plot)
}


head(rownames(predicted_gam_mat))
gene <- "Mbp"

gene %in% rownames(predicted_gam_mat)
gene %in% rownames(predicted_gam_atac_mat)


getAggMat_single_gene <- function(expr_vec, time_point, spatial_bin) {
  agg_result <- aggregate(expr_vec, by = list(Age = time_point, CellType = spatial_bin), FUN = mean)
  
  design_grid <- expand.grid(
    Age = c("P0", "P2", "P5","P7", "P10", "P21"),
    CellType = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
  )
  
  full_df <- merge(design_grid, agg_result, by = c("Age", "CellType"), all.x = TRUE)
  
  zscore <- scale(full_df$x)
  full_df$z <- as.numeric(zscore)
  
  full_df$Age <- factor(full_df$Age, levels = c("P0", "P2", "P5", "P7", "P10", "P21"))
  full_df$CellType <- factor(full_df$CellType, levels = c( "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"))
  
  return(full_df)
}


gene <- gene1[1]
gene<-"Mbp"
agg_rna_df  <- getAggMat_single_gene(predicted_gam_mat[gene, ], refined_corpus_region_idx$Age, refined_corpus_region_idx$CellType)
agg_atac_df <- getAggMat_single_gene(predicted_gam_atac_mat[gene, ], refined_corpus_region_idx$Age, refined_corpus_region_idx$CellType)

agg_rna_df$modality  <- "RNA"
agg_atac_df$modality <- "ATAC"
combined_df <- rbind(agg_rna_df, agg_atac_df)


celltype_colors <- c("ExciteL23" = "#ff05e1", "ExciteL5" = "#ffa6e3", 
                     "ExciteL4" = "#0091ff","ExciteL6" = "#009e00")
library(patchwork)

p_rna <- ggplot(agg_rna_df, aes(x = Age, y = z, group = CellType, color = CellType)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(values = celltype_colors) +
  labs(title = paste0(gene, " — RNA dynamics"),
       y = "Z-score", x = "Time", color = "Cell type") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid       = element_blank(),                           
    axis.line        = element_line(color = "black"),             
    axis.ticks       = element_line(color = "black"),             
    axis.text        = element_text(color = "black"),             
    plot.background  = element_rect(fill = "white", color = NA)  
  )


p_rna <- ggplot(agg_rna_df, aes(x = as.numeric(Age), y = z, group = CellType, color = CellType, fill = CellType)) +
  geom_smooth(se = TRUE, method = "lm", size = 1.5, alpha = 0.2) +  
  geom_point(size = 1.5) +
  scale_color_manual(values = celltype_colors) +
  scale_fill_manual(values = celltype_colors) + 
  labs(title = paste0(gene, " — RNA dynamics (linear fit)"),
       y = "Z-score", x = "Time", color = "Cell type", fill = "Cell type") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.text        = element_text(color = "black"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

p_atac <- ggplot(agg_atac_df, aes(x = as.numeric(Age), y = z, group = CellType, color = CellType, fill = CellType)) +
  geom_smooth(se = TRUE, method = "lm", size = 1.5, alpha = 0.2) +  
  geom_point(size = 1.5) +
  scale_color_manual(values = celltype_colors) +
  scale_fill_manual(values = celltype_colors) +  
  labs(title = paste0(gene, " — RNA dynamics (linear fit)"),
       y = "Z-score", x = "Time", color = "Cell type", fill = "Cell type") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.text        = element_text(color = "black"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

fit_stats <- agg_rna_df %>%
  mutate(Age_num = as.numeric(Age)) %>%
  group_by(CellType) %>%
  group_modify(~ {
    fit <- lm(z ~ Age_num, data = .x)
    tibble(
      slope = coef(fit)[2],
      r_squared = summary(fit)$r.squared
    )
  })



label_df <- fit_stats %>%
  mutate(
    label = paste0(CellType, ": slope = ", round(slope, 2), ", R² = ", round(r_squared, 2)),
    x = 1.1,  
    y = seq(2.2, 0.8, length.out = n())  
  )

p_rna <- p_rna +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label, color = CellType),
    inherit.aes = FALSE,
    hjust = 0, size = 3.5
  )

p_rna


fit_stats <- agg_atac_df %>%
  mutate(Age_num = as.numeric(Age)) %>%
  group_by(CellType) %>%
  group_modify(~ {
    fit <- lm(z ~ Age_num, data = .x)
    tibble(
      slope = coef(fit)[2],
      r_squared = summary(fit)$r.squared
    )
  })



label_df <- fit_stats %>%
  mutate(
    label = paste0(CellType, ": slope = ", round(slope, 2), ", R² = ", round(r_squared, 2)),
    x = 1.1,  
    y = seq(2.2, 0.8, length.out = n()) 
  )

p_atac <- p_atac +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label, color = CellType),
    inherit.aes = FALSE,
    hjust = 0, size = 3.5
  )

p_atac



p_rna <- ggplot(agg_rna_df, aes(x = as.numeric(Age), y = z, group = CellType, color = CellType, fill = CellType)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 1.5, linetype = "solid") +  
  geom_smooth(se = TRUE, method = "lm", size = 1.2, alpha = 0.2, linetype = "dashed") +    
  geom_point(size = 2) +
  scale_color_manual(values = celltype_colors) +
  scale_fill_manual(values = celltype_colors) +
  labs(title = paste0(gene, " — RNA dynamics (loess + linear fit)"),
       y = "Z-score", x = "Time", color = "Cell type", fill = "Cell type") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.text        = element_text(color = "black"),
    plot.background  = element_rect(fill = "white", color = NA)
  )




p_atac <- ggplot(agg_atac_df, aes(x = Age, y = z, group = CellType, color = CellType)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(values = celltype_colors) +
  labs(title = paste0(gene, " —  ATAC dynamics"),
       y = "Z-score", x = "Time", color = "Cell type") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid       = element_blank(),                           
    axis.line        = element_line(color = "black"),             
    axis.ticks       = element_line(color = "black"),             
    axis.text        = element_text(color = "black"),             
    plot.background  = element_rect(fill = "white", color = NA)  
  )








