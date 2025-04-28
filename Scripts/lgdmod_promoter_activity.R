
setwd("./G4AP/")
set.seed(1)
options(scipen = 999)

library(data.table)
library(dplyr)
library(stringr)
library(proActiv)
library(edgeR)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(ggrepel)

GrG4 <- function(G4H_path = NULL) {
  G4 <- data.table::fread(G4H_path, sep = "\t", header = FALSE) %>% data.frame
  G4 <- G4[, c(1:3, 5, 7)]
  colnames(G4) <- c("chr", "start", "end", "strand", "score")
  G4 <- GenomicRanges::makeGRangesFromDataFrame(G4,
                                                keep.extra.columns = TRUE, 
                                                ignore.strand = FALSE)
  return(G4)
}

PltVolcano <- function(tmp_ds_tss, save_path) {
  tmp_ds_tss$geneId <- sub("\\..*", "", tmp_ds_tss$geneId)
  
  tmp_ds_tss$significance <- ifelse(
    tmp_ds_tss$FDR < 0.05 & tmp_ds_tss$P.Value < 0.05 & abs(tmp_ds_tss$logFC) > 1, 
    "Both P and FDR Significant", 
    ifelse(tmp_ds_tss$P.Value < 0.05 & abs(tmp_ds_tss$logFC) > 1, 
           "P Significant", "Not Significant"))

  tmp_ds_tss$significance <- factor(tmp_ds_tss$significance, 
                                    levels = c("Both P and FDR Significant", 
                                               "P Significant", 
                                               "Not Significant"))

  p <- ggplot(tmp_ds_tss, aes(x = logFC, y = -log10(P.Value), color = significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Both P and FDR Significant" = "#c01d4f", 
                                  "P Significant" = "#004972", 
                                  "Not Significant" = "#939192")) +
    theme_minimal() +
    labs(x = "log2 Fold Change", y = "-log10(P-value)") +
    theme(legend.title = element_blank())

  ggsave(save_path, p, width = 8, height = 6)
}

PltDens <- function(data, x_var, group_var, 
                    x_label = "Distance to nearest G4 (bp)", 
                    y_label = "Density") {
  
  colors <- c("treat" = "#ee2026", "ctrl" = "#367abe")
  
  density_plot <- ggplot(data, aes(x = .data[[x_var]], fill = .data[[group_var]], color = .data[[group_var]])) +
    geom_density(alpha = 0.85, linewidth = 0.6, color = "black") +
    scale_fill_manual(values = colors) +
    labs(x = x_label, y = y_label) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      legend.position = "top",
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  box_plot <- ggplot(data, aes(y = .data[[x_var]], x = .data[[group_var]], fill = .data[[group_var]])) +
    geom_boxplot(width = 0.4, alpha = 0.85, outlier.shape = NA) +
    geom_jitter(aes(color = .data[[group_var]]), width = 0.2, size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  combined_plot <- density_plot + box_plot + plot_layout(widths = c(3, 1))
  
  return(combined_plot)
}

CalNrstG4Dist <- function(ligand_df, pG4) {
  
  pds_gr <- GRanges(
    seqnames = ligand_df$seqnames,
    ranges = IRanges(start = ligand_df$start, end = ligand_df$end),
    strand = ligand_df$strand
  )
  nearest_g4 <- distanceToNearest(pds_gr, pG4, ignore.strand = TRUE)
  ligand_df$nearest_g4_distance <- mcols(nearest_g4)$distance
  nearest_g4_idx <- subjectHits(nearest_g4)
  ligand_df$nearest_g4_start <- start(pG4)[nearest_g4_idx]
  ligand_df$nearest_g4_end <- end(pG4)[nearest_g4_idx]
  ligand_df$nearest_g4_chr <- as.character(seqnames(pG4)[nearest_g4_idx])
  ligand_df$nearest_g4_score <- score(pG4)[nearest_g4_idx]
  
  return(ligand_df)
}

sample_info_path <- "./data/GSE133419_sample_info.csv"
gtf_path <- "./data/ref/gencode.v19.annotation.gtf"

prmtr_anno <- preparePromoterAnnotation(file = gtf_path, 
                                        species = 'Homo_sapiens')
sample_info <- fread(sample_info_path, sep = ",", header = TRUE) %>% data.frame

pG4 <- GrG4(G4H_path = "./data/G4/G4hunter_w25_s1.5.txt")

pv <- NULL
prmtr_acti <- list()
prmtr_info <- list()
de_list <- list()

sample_info$path <- paste0("./data/RNA-seq/", 
                           unique(sample_info$AccessionID1), 
                           "/align/", 
                           unique(sample_info$AccessionID2, 
                                  "_SJ.out.tab"))
treat_lgd <- setdiff(unique(sample_info$Ligand), "Control")

for (tmp_treat in treat_lgd) {
  
  message(tmp_treat, " start...")
  tmp_ct_path <- sample_info %>% filter(Ligand == "Control") %>% dplyr::select(path) %>% unlist %>% as.vector %>% paste0(., "_SJ.out.tab")
  tmp_tt_path <- sample_info %>% filter(Ligand == tmp_treat) %>% dplyr::select(path) %>% unlist %>% as.vector %>% paste0(., "_SJ.out.tab")
  tmp_files <- c(tmp_ct_path, tmp_tt_path)
    
  tmp_condition <- c(rep("Control", length(tmp_ct_path)), rep(tmp_treat, length(tmp_tt_path))) 
  
  tmp_activ <- proActiv(files = tmp_files, 
                        promoterAnnotation = prmtr_anno,
                        condition = tmp_condition)
  tmp_activ <- tmp_activ[complete.cases(assays(tmp_activ)$promoterCounts),]
  tmp_activ <- tmp_activ[rowData(tmp_activ)$internalPromoter == FALSE, ]
  
  tmp_activ <- tmp_activ[rowData(tmp_activ)$seqnames %in% paste0("chr", c(1:22, "X")), ]
  prmtr_acti[[tmp_treat]] <- assays(tmp_activ)$absolutePromoterActivity %>% data.frame()
  prmtr_info[[tmp_treat]] <- rowData(tmp_activ) %>% data.frame()
 
  tmp_dge <- DGEList(counts=assays(tmp_activ)$promoterCounts, group=colData(tmp_activ)$condition, 
                       genes=rowData(tmp_activ)[, c("geneId", "promoterId")])
    
  tmp_keep <- rowSums(cpm(tmp_dge) > 1 ) >= 1
  tmp_dge <- tmp_dge[tmp_keep, , keep.lib.sizes = FALSE]
  tmp_dge <- calcNormFactors(tmp_dge)
  tmp_design = model.matrix(~colData(tmp_activ)$condition)
  tmp_disp <- estimateDisp(tmp_dge, design = tmp_design, tagwise = FALSE)
  tmp_QLfit <- glmQLFit(tmp_disp, design = tmp_design, robust = TRUE)
  tmp_qlf <- glmQLFTest(tmp_QLfit)
  tmp_qlf_total <- topTags(tmp_qlf, n = 1000000) %>% data.frame
  tmp_ds <- diffSpliceDGE(tmp_QLfit, geneid = "geneId")
  tmp_ds_tss <- topSpliceDGE(tmp_ds, test = "exon", n = 1000000)
  
  PltVolcano(tmp_ds_tss, paste0("./figure/ligand/ligand_volcano_", tmp_treat,".pdf"))
  
  x1 <- tmp_ds_tss %>% filter(logFC > log2(2), FDR < 0.05)
  x2 <- tmp_ds_tss %>% filter(logFC < -log2(2), FDR < 0.05)
  x3 <- length(intersect(unique(x1$geneId), unique(x2$geneId)))
  x4 <- length(unique(x1$geneId)) - x3
  x5 <- length(unique(x2$geneId)) - x3
  pv <- bind_rows(pv, data.frame(ligand = tmp_treat, 
                                 up = x4, down = x5, intersection = x3))
  
  de_list[[tmp_treat]] <- tmp_ds_tss %>% filter(abs(logFC) > log2(2), FDR < 0.05)
  message(tmp_treat, " end...")
}

# promoter location
all_prmtr <- data.frame(prmtr_anno@promoterCoordinates)
all_prmtr <- all_prmtr %>% filter(seqnames %in% paste0("chr", c(1:22, "X")))
de_list$PDS <- left_join(de_list$PDS, all_prmtr, by = "promoterId")
de_list$PhenDC3 <- left_join(de_list$PhenDC3, all_prmtr, by = "promoterId")
de_list$PhenDC6 <- left_join(de_list$PhenDC6, all_prmtr, by = "promoterId")

de_list_PDS <- CalNrstG4Dist(de_list$PDS, pG4)
de_list_PhenDC3 <- CalNrstG4Dist(de_list$PhenDC3, pG4)
de_list_PhenDC6 <- CalNrstG4Dist(de_list$PhenDC6, pG4)

ct_PDS_dist <- CalNrstG4Dist(all_prmtr[sample(1:nrow(all_prmtr), nrow(de_list_PDS), replace = FALSE), ], pG4)
ct_PhenDC3_dist <- CalNrstG4Dist(all_prmtr[sample(1:nrow(all_prmtr), nrow(de_list_PhenDC3), replace = FALSE), ], pG4)
ct_PhenDC6_dist <- CalNrstG4Dist(all_prmtr[sample(1:nrow(all_prmtr), nrow(de_list_PhenDC6), replace = FALSE), ], pG4)

PDS_dist <- bind_rows(data.frame(distance = de_list_PDS$nearest_g4_distance, ligand = "treat"),
                      data.frame(distance = ct_PDS_dist$nearest_g4_distance, ligand = "ctrl"))
PDS_dist$ligand <- factor(PDS_dist$ligand, level = c("treat", "ctrl"))
p1 <- PltDens(PDS_dist, "distance", "ligand")
wilcox.test(de_list_PDS$nearest_g4_distance, ct_PDS_dist$nearest_g4_distance, alternative = "less")

PhenDC3_dist <- bind_rows(data.frame(distance = de_list_PhenDC3$nearest_g4_distance, ligand = "treat"),
                      data.frame(distance = ct_PhenDC3_dist$nearest_g4_distance, ligand = "ctrl"))
PhenDC3_dist$ligand <- factor(PhenDC3_dist$ligand, level = c("treat", "ctrl"))
p2 <- PltDens(PhenDC3_dist, "distance", "ligand")
wilcox.test(de_list_PhenDC3$nearest_g4_distance, ct_PhenDC3_dist$nearest_g4_distance, alternative = "less")

PhenDC6_dist <- bind_rows(data.frame(distance = de_list_PhenDC6$nearest_g4_distance, ligand = "treat"),
                          data.frame(distance = ct_PhenDC6_dist$nearest_g4_distance, ligand = "ctrl"))
PhenDC6_dist$ligand <- factor(PhenDC6_dist$ligand, level = c("treat", "ctrl"))
p3 <- PltDens(PhenDC6_dist, "distance", "ligand")
wilcox.test(de_list_PhenDC6$nearest_g4_distance, ct_PhenDC6_dist$nearest_g4_distance, alternative = "less")

ggsave("./figure/ligand/G4_distance.pdf", p1 / p2 / p3, width = 8.5, height = 8)
