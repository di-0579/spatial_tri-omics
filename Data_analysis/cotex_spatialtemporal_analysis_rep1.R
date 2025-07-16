#####################################################################
#####################    Settings    ################################
#####################################################################

set.seed(2023)  
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
library(reshape2)

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

input_dir <- "/input_dir/"
figures_dir <- "/cortex_final/"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}


data_dir <- file.path(figures_dir, "data")
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

obj_dir <- file.path(figures_dir, "object")
if (!dir.exists(obj_dir)) {
  dir.create(obj_dir, recursive = TRUE)
}



#####################################################################
##################### Fit GAM model for RNA-seq #####################
#####################################################################

refined_corpus_region_idx <- read.csv(file.path(input_dir, "refined_cortex_region_idx.csv"),
                                      row.names = 1)

。
dorcMat.s <- readRDS(file.path(input_dir, "refined_cortex_region_DORC_mat.rds"))

refined_corpus_region_idx <- refined_corpus_region_idx[colnames(dorcMat.s), , drop = FALSE]


seu_obj_corpus <- readRDS(file.path(input_dir, "cortex_seu_byCluster.rds"))
library_id <- colSums(seu_obj_corpus)

names(library_id) <- gsub("_", "#", names(library_id))
library_id <- library_id[colnames(dorcMat.s)]  #

##############################################################
exprs_mat <- seu_obj_corpus@assays$SCT@counts
colnames(exprs_mat) <- gsub("_", "#", colnames(exprs_mat))
exprs_mat <- exprs_mat[, colnames(dorcMat.s)]


sample <- factor(refined_corpus_region_idx$Sample)

spatial <- factor(refined_corpus_region_idx$CellType,
                  levels = c("ExciteL6b", "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"))
spatial <- as.numeric(spatial)

time <- as.numeric(gsub("P", "", refined_corpus_region_idx$Age))
names(time) <- colnames(dorcMat.s)

keep_genes <- names(which(rowMeans(exprs_mat != 0) > 0.02))

library(pbmcapply)
gam_fit <- pbmcapply::pbmclapply(keep_genes, function(gene) {
  y_exprs <- exprs_mat[gene, ]
  gam(y_exprs ~ library_id + ti(time) + ti(spatial) + ti(time, spatial), 
      family = nb())
}, mc.cores = 8)


names(gam_fit) <- keep_genes


gam_fit_library_id <- pbmcapply::pbmclapply(keep_genes, function(gene) {
  y_exprs <- exprs_mat[gene, ]
  gam_fit_null <- gam(y_exprs ~ library_id, 
                      family = nb())
}, mc.cores = 8)

names(gam_fit_library_id) <- keep_genes


#####################################################################
##################### Fit GAM model for ATAC-seq ####################
#####################################################################

signac_obj_cortex <- readRDS(file.path(input_dir,"cortex_signac_byCluster.rds"))

exprs_mat <- signac_obj_cortex@assays$ATAC@counts
library_id <- colSums(exprs_mat)
names(library_id) <- gsub("_", "#", names(library_id))
library_id <- library_id[colnames(dorcMat.s)]  

head(library_id)
head(colnames(dorcMat.s))
head(names(time))

gam_fit_atac <- pbmcapply::pbmclapply(1:nrow(dorcMat.s), function(gene) {
  y_exprs <- log2(dorcMat.s[gene, ] + 1)
  data <- data.frame(y_exprs, time, spatial, library_id)
  gamma_model <- gam(y_exprs ~ ti(time) + ti(spatial) + ti(time, spatial), 
                     data = data)
}, mc.cores = 30)

names(gam_fit_atac) <- rownames(dorcMat.s)


sum(sapply(gam_fit_atac, is.null))


gam_fit_atac_success <- gam_fit_atac[!sapply(gam_fit_atac, is.null)]

failed_genes <- names(gam_fit_atac)[sapply(gam_fit_atac, is.null)]

retry_round <- 1
while (length(failed_genes) > 0) {
  cat("Retry round:", retry_round, "- Remaining failed genes:", length(failed_genes), "\n")
  
  gam_fit_retry <- pbmcapply::pbmclapply(failed_genes, function(gene) {
    y_exprs <- log2(dorcMat.s[gene, ] + 1)
    data <- data.frame(y_exprs, time, spatial, library_id)
    gamma_model <- tryCatch(
      gam(y_exprs ~ ti(time) + ti(spatial) + ti(time, spatial), data = data),
      error = function(e) NULL
    )
    gamma_model
  }, mc.cores = 30)
  
  names(gam_fit_retry) <- failed_genes
  
  gam_fit_retry_success <- gam_fit_retry[!sapply(gam_fit_retry, is.null)]
  gam_fit_atac_success <- c(gam_fit_atac_success, gam_fit_retry_success)
  
  failed_genes <- names(gam_fit_retry)[sapply(gam_fit_retry, is.null)]
  
  retry_round <- retry_round + 1
}

cat(":", length(gam_fit_atac_success), "\n")




gam_fit_atac_library_id <- pbmcapply::pbmclapply(1:nrow(dorcMat.s), function(gene) {
  y_exprs <- log2(dorcMat.s[gene, ] + 1)
  data <- data.frame(y_exprs, time, spatial, library_id)
  gamma_model <- gam(y_exprs ~ library_id,
                     data = data)
}, mc.cores = 30)


names(gam_fit_atac_library_id) <- rownames(dorcMat.s)

#####################################################################
##################### Analysis of GAM results - RNA #################
#####################################################################

length(gam_fit) == length(gam_fit_library_id)
all(names(gam_fit) == names(gam_fit_library_id))


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
。
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
  rownames(df_genes) %in% names(which(adjr_squred > 0.02)) & 
  rownames(df_genes) %in% names(which(lrt_test_pvalue_adj < 0.01))
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

####################################################################
##################### Analysis of GAM results - ATAC ###############
####################################################################

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
  rownames(df_atac) %in% names(which(adjr_squred_atac > 0.02)) &
  lrt_test_atac_pvalue_adj < 0.01 &
  !is.na(lrt_test_atac_pvalue_adj)
selected_peaks <- rownames(df_atac)[selected_peaks]
length(selected_peaks)


predicted_gam_atac_mat <- lapply(gam_fit_atac, function(x) predict(x, type = "response"))
predicted_gam_atac_mat <- do.call(rbind, predicted_gam_atac_mat)
colnames(predicted_gam_atac_mat) <- gsub("_", "#", colnames(predicted_gam_atac_mat))

all(colnames(predicted_gam_atac_mat) == rownames(refined_corpus_region_idx))
setequal(colnames(predicted_gam_atac_mat), rownames(refined_corpus_region_idx))

predicted_gam_atac_mat <- predicted_gam_atac_mat[, rownames(refined_corpus_region_idx)]

#####################################################################
##################### Joint clustering  #############################
####################################################################

common_genes <- intersect(intersect(names(gam_fit), names(gam_fit_atac)),
                          union(selected_genes, selected_peaks))
length(common_genes)

design_grid <- expand.grid(
  Group.1 = c("P0", "P2", "P5", "P10", "P21"),
  Group.2 = c("ExciteL6b", "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
)

agg_values <- sapply(1:nrow(design_grid), function(i) {
  age      <- design_grid$Group.1[i]
  celltype <- design_grid$Group.2[i]
  idx <- which(refined_corpus_region_idx$Age == age & refined_corpus_region_idx$CellType == celltype)
  mean(predicted_gam_mat[1, idx], na.rm = TRUE)
})

agg_df <- cbind(design_grid, value = agg_values)

agg_df_combine <- rbind(
  design_grid,  # RNA
  design_grid   # ATAC
)
agg_df_combine$modality <- c(
  rep("RNA", nrow(design_grid)),
  rep("ATAC", nrow(design_grid))
)

head(agg_df_combine)

getAggMat <- function(predicted_gam_mat, time_point, spatial_bin) {
  design_grid <- expand.grid(
    Age = c("P0", "P2", "P5", "P10", "P21"),
    CellType = c("ExciteL6b", "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
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
RNA_big_cluster <- cutree(h_rna, k = 13)
names(RNA_big_cluster) <- cluster_ids
h_atac <- hclust(dist(atac_cluster_avg), method = "ward.D2")
ATAC_big_cluster <- cutree(h_atac, k = 7)
names(ATAC_big_cluster) <- cluster_ids
combine_hclust_cluster_big <- paste0(
  "R", RNA_big_cluster[as.character(combine_hclust_cluster)],
  "_A", ATAC_big_cluster[as.character(combine_hclust_cluster)]
)


names(combine_hclust_cluster) <- rownames(predicted_gam_combine_mat_agg_scale)
names(combine_hclust_cluster_big) <- rownames(predicted_gam_combine_mat_agg_scale)
head(names(combine_hclust_cluster) <- rownames(predicted_gam_combine_mat_agg_scale))


 table(combine_hclust_cluster_big)
#####################################################################
##################### Visualization.    #############################
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
  levels = c("P0", "P2", "P5", "P10", "P21")
)

agg_df_mat$Group.2 <- factor(
  agg_df_mat$Group.2,
  levels = c("ExciteL6b", "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23")
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


agg_spot$Group.1  <- factor(agg_spot$Group.1, levels = c("P0", "P2", "P5", "P10", "P21"))
agg_spot$Group.2  <- factor(agg_spot$Group.2, levels = c("ExciteL6b", "ExciteL6", "ExciteL5", "ExciteL4", "ExciteL23"))
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



predicted_gam_combine_mat_agg_scale_mean_big <- sapply(
  paste0("R", 1:13),  # 大簇标签
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
  levels = c("P21", "P10", "P5", "P2", "P0")
)

ggplot(
  agg_df_mat[agg_df_mat$modality == "RNA", ], 
  aes(
    x = factor(Group.1, levels = rev(levels(agg_df_mat$Group.1))),  
    y = factor(Group.2, levels = c(                               
      "ExciteL6b",
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
  paste0("A", 1:7),
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
  levels = c("P21", "P10", "P5", "P2", "P0")
)



ggplot(
  agg_df_mat[agg_df_mat$modality == "ATAC", ], 
  aes(
    x = factor(Group.1, levels = rev(levels(agg_df_mat$Group.1))),
    y = factor(Group.2, levels = c(
      "ExciteL6b",
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
                                                                          "P2",  "P5",
                                                                          "P10", "P21"))


agg_df_mat_spatial$Age <- factor(
  agg_df_mat_spatial$Age,
  levels = c("P0","P2","P5","P10","P21")
)

agg_df_mat_spatial$CellType <- factor(
  agg_df_mat_spatial$CellType,
  levels = c("ExciteL6b","ExciteL6","ExciteL5","ExciteL4","ExciteL23")
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



for (i in names(st_gene_lsit)) {
  print(i)
  dir.create(paste0(figures_dir, "gene_figures_cortex/", i, "/"), recursive = TRUE)
  for (genes in st_gene_lsit[[i]]) {
    message("Plotting gene: ", genes)
    
    p <- tryCatch(plotEachGene(genes), error = function(e) {
      message("plotEachGene failed for ", genes, ": ", conditionMessage(e))
      return(NULL)
    })
    
    if (!is.null(p)) {
      tryCatch({
        ggsave(paste0(figures_dir, "gene_figures_cortex/", i, "/", genes, ".pdf"),
               plot = p, width = 16, height = 12)
      }, error = function(e) {
        message("ggsave failed for ", genes, ": ", conditionMessage(e))
      })
    }
    
    rm(p)
    gc()
  }
}





plotHeatmapEachGene <- function(gene) {
  split_info <- strsplit(colnames(predicted_gam_rna_mat_agg_scale), "_")
  meta_df <- do.call(rbind, lapply(split_info, function(x) data.frame(Age = x[1], CellType = x[2])))
  meta_df$Age <- as.character(meta_df$Age)
  meta_df$CellType <- as.character(meta_df$CellType)
  
  Age_levels <- c("P0", "P2", "P5", "P10", "P21")
  CellType_levels <- c("ExciteL23", "ExciteL4", "ExciteL5","ExciteL6", "ExciteL6b" )
  
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

dir.create(file.path(figures_dir, "cortex_heatmaps"), recursive = TRUE, showWarnings = FALSE)

plotHeatmapGeneList(gene_list, save_dir = file.path(figures_dir, "cortex_heatmaps"))