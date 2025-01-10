#####################################################################
#####################    Settings    ################################
#####################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
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

addArchRGenome("mm10")
ncore <- 8
addArchRThreads(threads = ncore)


pign <- colorRampPalette(c(rev(brewer.pal(n = 7, name = "RdBu"))[1:3],
                           rev(brewer.pal(n = 7, name = "PRGn"))[c(4:7)]))(100)
brbg <- colorRampPalette(c(rev(brewer.pal(n = 7, name = "PRGn"))[1:3],
                           rev(brewer.pal(n = 7, name = "BrBG"))[c(4:7)]))(100)
pugn <- colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(100)


ggplot2::theme_set(theme_bw() + theme_yx() + 
                     theme(axis.text.y = element_text(color = "black"),
                           axis.text.x = element_text(color = "black")) )

input_dir <- "output/P21/"
figures_dir <- "figures/P21_final/"
dir.create(figures_dir, recursive = TRUE)





#####################################################################
##################### Fit GAM model for RNA-seq #####################
#####################################################################



# Seurat object
seu_obj_corpus <- readRDS(file.path(input_dir, "MouseDev_11sample_corpus_seu_byCluster_P21.rds"))
library_id <- colSums(seu_obj_corpus)

# Spatial bin 
refined_corpus_region_idx <- read.csv(file.path(input_dir, "MouseDev_11sample_refined_corpus_region_idx_P21.csv"),
                                      row.names = 1)


exprs_mat <- seu_obj_corpus@assays$SCT@counts
sample <- factor(refined_corpus_region_idx$Sample)
spatial <- (refined_corpus_region_idx$final_bin)
time <- as.numeric(gsub("P", "", seu_obj_corpus$Timepoint))
keep_genes <- names(which(rowMeans(exprs_mat != 0) > 0.02))
gam_fit <- pbmcapply::pbmclapply(keep_genes, function(gene) {
  y_exprs <- exprs_mat[gene, ]
  
  gam(y_exprs ~ library_id + ti(time) + ti(spatial) + ti(time, spatial), 
      family = nb())
}, mc.cores = 40)
names(gam_fit) <- keep_genes
saveRDS(gam_fit, file = file.path(input_dir, "MouseDev_11sample_refined_corpus_region_GAM_fit_P21.rds"))



gam_fit_library_id <- pbmcapply::pbmclapply(keep_genes, function(gene) {
  y_exprs <- exprs_mat[gene, ]
  gam_fit_null <- gam(y_exprs ~ library_id, 
                      family = nb())
}, mc.cores = 40)
names(gam_fit_library_id) <- keep_genes
saveRDS(gam_fit_library_id, file = file.path(input_dir, "MouseDev_11sample_refined_corpus_region_GAM_fit_P21_library_id.rds"))



#####################################################################
##################### Fit GAM model for ATAC-seq ####################
#####################################################################


# DORC matrix
dorcMat.s <- readRDS(file.path(input_dir, "MouseDev_11sample_refined_corpus_region_DORC_mat_P21.rds"))

# Spatial bin 
refined_corpus_region_idx <- read.csv(file.path(input_dir, "MouseDev_11sample_refined_corpus_region_idx_P21.csv"),
                                      row.names = 1)

sample <- factor(refined_corpus_region_idx[colnames(dorcMat.s), ]$Sample)
spatial <- (refined_corpus_region_idx[colnames(dorcMat.s), ]$final_bin)
time <- as.numeric(gsub("P", "", seu_obj_corpus[, colnames(dorcMat.s)]$Timepoint))
library_id <- readRDS(file.path(input_dir, "MouseDev_11sample_refined_corpus_region_DORC_mat_P21_library_id.rds"))
names(library_id) <- gsub("_", "#", names(library_id))
library_id <- library_id[colnames(dorcMat.s)]

gam_fit_atac <- pbmcapply::pbmclapply(1:nrow(dorcMat.s), function(gene) {
  y_exprs <- log2(dorcMat.s[gene, ] + 1)
  data <- data.frame(y_exprs, time, spatial, library_id)
  gamma_model <- gam(y_exprs ~ ti(time) + ti(spatial) + ti(time, spatial), 
                     data = data)
}, mc.cores = 40)

names(gam_fit_atac) <- rownames(dorcMat.s)
saveRDS(gam_fit_atac, file = file.path(input_dir, "MouseDev_11sample_refined_corpus_region_DORC_GAM_fit_P21.rds"))


gam_fit_atac_library_id <- pbmcapply::pbmclapply(1:nrow(dorcMat.s), function(gene) {
  y_exprs <- log2(dorcMat.s[gene, ] + 1)
  data <- data.frame(y_exprs, time, spatial, library_id)
  gamma_model <- gam(y_exprs ~ library_id,
                     data = data)
  
  
  
}, mc.cores = 40)
names(gam_fit_atac_library_id) <- rownames(dorcMat.s)
saveRDS(gam_fit_atac_library_id, file = file.path(input_dir, "MouseDev_11sample_refined_corpus_region_DORC_GAM_fit_P21_library_id.rds"))




#####################################################################
##################### Analysis of GAM results - RNA #################
#####################################################################


# Perform LRT
lrt_test <- lapply(names(gam_fit), function(x) {
  anova(gam_fit[[x]], gam_fit_library_id[[x]], test = "LRT")
})
lrt_test_pvalue <- lapply(lrt_test, function(x) x$`Pr(>Chi)`[2])
lrt_test_pvalue <- unlist(lrt_test_pvalue)
lrt_test_pvalue_adj <- p.adjust(lrt_test_pvalue, method = "BH")
names(lrt_test_pvalue_adj) <- names(gam_fit)




# Select genes
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





df_genes <- data.frame(temporal = sum_s_table_time$adj.p.value < 0.01 , 
                       spatial = sum_s_table_spatial$adj.p.value < 0.01 , 
                       spatiotemporal = sum_s_table_spacetime$adj.p.value < 0.01)
df_genes$category <- paste(df_genes$temporal, df_genes$spatial, df_genes$spatiotemporal, sep = "_")
rownames(df_genes) <- rownames(sum_s_table_time)
gene_list <- split(rownames(df_genes), df_genes$category)
table(df_genes$category)

selected_genes <- df_genes$category != "FALSE_FALSE_FALSE" & 
  rownames(df_genes) %in% names(which(adjr_squred > 0.01)) & 
  rownames(df_genes) %in% names(which(lrt_test_pvalue_adj < 0.01))
selected_genes <- rownames(df_genes)[selected_genes]
length(selected_genes)
grep("Rpl|Rps|^mt-", selected_genes, value = TRUE)

selected_genes <- selected_genes[!grepl("^mt-", selected_genes)]



# Predict based on the GAM model
predicted_gam_mat <- lapply(gam_fit, function(x) predict(x, type = "response"))
predicted_gam_mat <- do.call(rbind, predicted_gam_mat)
predicted_gam_mat <- log2(predicted_gam_mat + 1)
predicted_gam_mat <- predicted_gam_mat[, rownames(refined_corpus_region_idx)]




#####################################################################
##################### Analysis of GAM results - ATAC ################
#####################################################################


# LRT
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
sum(sum_s_table_time_atac$adj.p.value < 0.01)
sum_s_table_spatial_atac <- do.call(rbind, lapply(sum_s_table_atac, function(x) x[2, ]))
sum_s_table_spatial_atac <- data.frame(sum_s_table_spatial_atac)
sum_s_table_spatial_atac$adj.p.value <- p.adjust(sum_s_table_spatial_atac$p.value, method = "BH")
sum(sum_s_table_spatial_atac$adj.p.value < 0.01)
sum_s_table_spacetime_atac <- do.call(rbind, lapply(sum_s_table_atac, function(x) x[3, ]))
sum_s_table_spacetime_atac <- data.frame(sum_s_table_spacetime_atac)
sum_s_table_spacetime_atac$adj.p.value <- p.adjust(sum_s_table_spacetime_atac$p.value, method = "BH")
sum(sum_s_table_spacetime_atac$adj.p.value < 0.01)

adjr_squred_atac <- lapply(gam_fit_atac, function(x) {
  summary(x)$r.sq
})
adjr_squred_atac <- unlist(adjr_squred_atac)


# Select peaks

df_atac <- data.frame(temporal = sum_s_table_time_atac$adj.p.value < 0.01 , 
                      spatial = sum_s_table_spatial_atac$adj.p.value < 0.01 , 
                      spatiotemporal = sum_s_table_spacetime_atac$adj.p.value < 0.01)
df_atac$category <- paste(df_atac$temporal, df_atac$spatial, df_atac$spatiotemporal, sep = "_")
rownames(df_atac) <- rownames(sum_s_table_time_atac)
gene_list <- split(rownames(df_atac), df_atac$category)
table(df_atac$category)

selected_peaks <- df_atac$category != "FALSE_FALSE_FALSE" & 
  rownames(df_atac) %in% names(which(adjr_squred_atac > 0.025)) &
  lrt_test_atac_pvalue_adj < 0.01 &
  !is.na(lrt_test_atac_pvalue_adj)
selected_peaks <- rownames(df_atac)[selected_peaks]
length(selected_peaks)



# GAM prediction
predicted_gam_atac_mat <- lapply(gam_fit_atac, function(x) predict(x, type = "response"))
predicted_gam_atac_mat <- do.call(rbind, predicted_gam_atac_mat)
colnames(predicted_gam_atac_mat) <- gsub("_", "#", colnames(predicted_gam_atac_mat))
predicted_gam_atac_mat <- predicted_gam_atac_mat[, rownames(refined_corpus_region_idx)]






#####################################################################
##################### Joint clustering  #############################
#####################################################################



common_genes <- intersect(intersect(names(gam_fit), names(gam_fit_atac)),
                          union(selected_genes, selected_peaks))
length(common_genes)



# Scale data

getAggMat <- function(predicted_gam_mat, time_point, spatial_bin) {
  predicted_mat_agg <- apply(predicted_gam_mat, 1, function(x) {
    aggregate(x, list(time_point, spatial_bin), mean)$x
  })
  
  predicted_mat_agg <- t(predicted_mat_agg)
  gene_means <- rowMeans(predicted_mat_agg)
  gene_sds <- apply(predicted_mat_agg, 1, sd)
  predicted_agg_mat_scale <- sweep(predicted_mat_agg, 1, gene_means, "-")
  predicted_agg_mat_scale <- sweep(predicted_agg_mat_scale, 1, gene_sds, "/")
  return(predicted_agg_mat_scale)
}

predicted_gam_rna_mat_agg_scale <- getAggMat(predicted_gam_mat[common_genes, ], 
                                             refined_corpus_region_idx$Timepoint, 
                                             refined_corpus_region_idx$final_bin)

predicted_gam_atac_mat_agg_scale <- getAggMat(predicted_gam_atac_mat[common_genes, ], 
                                              refined_corpus_region_idx$Timepoint, 
                                              refined_corpus_region_idx$final_bin)

predicted_gam_combine_mat_agg_scale <- cbind(predicted_gam_rna_mat_agg_scale,
                                             predicted_gam_atac_mat_agg_scale)


# Joint clustering

hclust_res <- hclust(1 - proxy::simil((predicted_gam_combine_mat_agg_scale)), method = "ward.D2")
combine_hclust_cluster <- cutree(hclust_res, k = 20)
table(combine_hclust_cluster)
sample_row <- data.frame(cluster = factor(combine_hclust_cluster))
rownames(sample_row) <- names(combine_hclust_cluster)

predicted_gam_combine_mat_agg_scale_mean <- sapply(names(table(combine_hclust_cluster)), function(x) {
  colMeans(predicted_gam_combine_mat_agg_scale[combine_hclust_cluster == x, ])
})



hclust_res <- hclust(dist(t(predicted_gam_combine_mat_agg_scale_mean[agg_df_combine$modality == "RNA", ])),
                     method = "ward.D2")
RNA_big_cluster <- cutree(hclust_res, k = 6)


atac_hclust_res <- hclust(dist(t(predicted_gam_combine_mat_agg_scale_mean[agg_df_combine$modality == "ATAC", ])),
                          method = "ward.D2")
ATAC_big_cluster <- cutree(atac_hclust_res, k = 4)





agg_df <- expand.grid(levels(refined_corpus_region_idx$Timepoint),
                      names(table(refined_corpus_region_idx$final_bin)))
colnames(agg_df) <- c("Group.1", "Group.2")
agg_df_combine <- rbind(agg_df[, c(1:2)],
                        agg_df[, c(1:2)])
agg_df_combine$modality <- c(rep("RNA", nrow(agg_df)),
                             rep("ATAC", nrow(agg_df)))


combine_hclust_cluster_big <- paste(plyr::mapvalues(combine_hclust_cluster, 
                                                    from = names(RNA_big_cluster), 
                                                    to = paste0("R", RNA_big_cluster)),
                                    plyr::mapvalues(combine_hclust_cluster, 
                                                    from = names(ATAC_big_cluster), 
                                                    to = paste0("A", ATAC_big_cluster)),
                                    sep = "_")
length(unique(combine_hclust_cluster_big))




#####################################################################
##################### Visualization.    #############################
#####################################################################

# Combine heatmap


predicted_gam_combine_mat_agg_scale_mean_big <- sapply(names(table(combine_hclust_cluster_big)), function(x) {
  colMeans(predicted_gam_combine_mat_agg_scale[combine_hclust_cluster_big == x, ])
})

agg_df_mat <- melt(data.frame(cbind(agg_df_combine, predicted_gam_combine_mat_agg_scale_mean_big)), 
                   id.vars = c("Group.1", "Group.2", "modality"))


ggplot(agg_df_mat, 
       aes(x = Group.1, y = factor(Group.2, levels = c(10:1)), 
           fill = value)) +
  geom_tile() +
  geom_tile(aes(fill = value),
            data = agg_df_mat[agg_df_mat$modality == "ATAC", ]) +
  scale_fill_gradientn(colours = pign, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "ATAC") +
  new_scale_fill() +
  geom_tile(aes(fill = value),
            data = agg_df_mat[agg_df_mat$modality == "RNA", ]) +
  scale_fill_gradientn(colours = rdbu, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "RNA") +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    strip.background = element_blank() ) +
  facet_grid(modality~variable) +
  theme(aspect.ratio = 1) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")


ggsave(file.path(figures_dir, "stSeq_combine_heatmap_P21_combined.pdf"), width = 25, height = 8)






# RNA pattern heatmap

predicted_gam_combine_mat_agg_scale_mean_big <- sapply(paste0("R", c(1:6)), function(x) { 
  colMeans(predicted_gam_combine_mat_agg_scale[grepl(x, combine_hclust_cluster_big), ])
})

agg_df_mat <- melt(data.frame(cbind(agg_df_combine, predicted_gam_combine_mat_agg_scale_mean_big)), 
                   id.vars = c("Group.1", "Group.2", "modality"))



ggplot(agg_df_mat[agg_df_mat$modality == "RNA", ], 
       aes(y = factor(Group.1, levels = rev(levels(agg_df_mat$Group.1))), 
           x = factor(Group.2, levels = c(1:10)), 
           fill = value)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours = rdbu, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "RNA") +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    strip.background = element_blank() ) +
  facet_grid(~variable) +
  theme(aspect.ratio = 1.4) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")


ggsave(file.path(figures_dir, "stSeq_combine_heatmap_P21_combined_RNA_pattern.pdf"), 
       width = 10, height = 4)


# ATAC pattern heatmap

predicted_gam_combine_mat_agg_scale_mean_big <- sapply(paste0("A", c(1:4)), function(x) { 
  colMeans(predicted_gam_combine_mat_agg_scale[grepl(x, combine_hclust_cluster_big), ])
})


agg_df_mat <- melt(data.frame(cbind(agg_df_combine, predicted_gam_combine_mat_agg_scale_mean_big)), 
                   id.vars = c("Group.1", "Group.2", "modality"))


ggplot(agg_df_mat[agg_df_mat$modality == "ATAC", ], 
       aes(y = factor(Group.1, levels = rev(levels(agg_df_mat$Group.1))), 
           x = factor(Group.2, levels = c(1:10)), 
           fill = value)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours = pign, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "ATAC") +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    strip.background = element_blank() ) +
  facet_grid(~variable) +
  theme(aspect.ratio = 1.4) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")


ggsave(file.path(figures_dir, "stSeq_combine_heatmap_P21_combined_ATAC_pattern.pdf"), 
       width = 10, height = 4)



# Combine spatial pattern
predicted_gam_combine_mat_agg_scale_mean_big_atac <- predicted_gam_combine_mat_agg_scale_mean_big[agg_df_combine$modality == "ATAC", ]
agg_df_mat_atac <- melt(data.frame(cbind(agg_df[, c(1:2)], predicted_gam_combine_mat_agg_scale_mean_big_atac)), id.vars = c("Group.1", "Group.2"))
agg_df_mat_atac$modality <- "ATAC"
combine_hclust_cluster_big_rna <- unlist(lapply(strsplit(combine_hclust_cluster_big, "_"), "[[", 1))
predicted_gam_combine_mat_agg_scale_mean_big_rna <- sapply(names(table(combine_hclust_cluster_big_rna)), function(x) {
  colMeans(predicted_gam_combine_mat_agg_scale[combine_hclust_cluster_big_rna == x, ])
})
predicted_gam_combine_mat_agg_scale_mean_big_rna <- predicted_gam_combine_mat_agg_scale_mean_big_rna[agg_df_combine$modality == "RNA", ]
agg_df_mat_rna <- melt(data.frame(cbind(agg_df[, c(1:2)], predicted_gam_combine_mat_agg_scale_mean_big_rna)), id.vars = c("Group.1", "Group.2"))
agg_df_mat_rna$modality <- "RNA"
agg_df_mat_combine <- rbind(agg_df_mat_rna, agg_df_mat_atac)
colnames(agg_df_mat_combine)[c(1:5)] <- c("Timepoint", "final_bin", "cluster", "value", "modality")


agg_df_mat_spatial <- merge(rbind(refined_corpus_region_idx, refined_corpus_region_idx), agg_df_mat_combine)
agg_df_mat_spatial$Sample <- factor(agg_df_mat_spatial$Sample, 
                                    levels = c("B8P0", "C8P0", "Q43P0",
                                               "B8P2", "Q43P2", "C8P5", "Q43P5",
                                               "C8P10", "Q43P10", "Q43P21", "C8P22"))




ggplot(agg_df_mat_spatial, 
       aes(x = x, y = y, fill = value)) +
  geom_tile() +
  geom_tile(aes(fill = value),
            data = agg_df_mat_spatial[agg_df_mat_spatial$modality == "ATAC", ]) +
  scale_fill_gradientn(colours = pign, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "ATAC") +
  new_scale_fill() +
  geom_tile(aes(fill = value),
            data = agg_df_mat_spatial[agg_df_mat_spatial$modality == "RNA", ]) +
  scale_fill_gradientn(colours = rdbu, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "RNA") +
  facet_grid(cluster +  modality~Timepoint) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12), 
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 0.6,
        strip.background = element_blank()) +
  ylab("Spatial regions") +
  xlab("Timepoints") +
  labs(fill = "Scaled expression")
ggsave(file.path(figures_dir, "stSeq_combine_spatial_P21_combined.pdf"), width = 10, height = 40)










## Ploting selected genes




gene_list <- c("Pcdh15", "Sema3c", "Ptprz1",
               "Vcan", "Sox6", "Ncan", "Pdgfra",
               "Ednrb", "Cspg5", "Kcnip3", "Bcas1",
               "Serpine2", "Cd9", "Gpr17", "Tnr",
               "Mpzl1", "Itpr2", "9630013A20Rik",
               "Rras2", "Cnksr3", "Ctps", "Glul", "Fyn", 
               "Onecut2", "Myrf", "S100b", "Car2", "Prr5l",
               "Cnp", "Egr1", "Mog", "Mobp", "Mag", "Mbp", "Opalin",
               "Kcna2", "Nfasc", 
               "Kcna1", "Cntn2")



agg_df_selected_rna <- cbind(agg_df[, c(1:2)], t(predicted_gam_rna_mat_agg_scale[gene_list, ]))
agg_df_selected_rna <- melt(agg_df_selected_rna, id.vars = c("Group.1", "Group.2"))
agg_df_selected_rna$variable <- factor(agg_df_selected_rna$variable, levels = rev(gene_list))
agg_df_selected_rna$Group.2 <- factor(agg_df_selected_rna$Group.2, levels = c(1:10))
g1 <- ggplot(agg_df_selected_rna, aes(x = Group.2, y = variable, 
                                      fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = rdbu, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "RNA") +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(aspect.ratio = 2) +
  ylab("") +
  xlab("") +
  facet_grid(~Group.1) +
  labs(fill = "Scaled expression") 

agg_df_selected_atac <- cbind(agg_df[, c(1:2)], t(predicted_gam_atac_mat_agg_scale[gene_list, ]))
agg_df_selected_atac <- melt(agg_df_selected_atac, id.vars = c("Group.1", "Group.2"))
agg_df_selected_atac$variable <- factor(agg_df_selected_atac$variable, levels = rev(gene_list))
agg_df_selected_atac$Group.2 <- factor(agg_df_selected_atac$Group.2, levels = c(1:10))
g2 <- ggplot(agg_df_selected_atac, aes(x = Group.2, y = variable, 
                                       fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = pign, limits = c(-2, 2), oob = scales::oob_squish,
                       name = "ATAC") +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(aspect.ratio = 2) +
  ylab("") +
  xlab("") +
  facet_grid(~Group.1) +
  labs(fill = "Scaled expression") 
ggarrange(g1, g2, ncol = 1, nrow = 2, align = "hv")
ggsave(file.path(figures_dir, "stSeq_corpus_selected_gene_st_heatmap_combined.pdf"), 
       width = 14, height = 12, limitsize = FALSE)




