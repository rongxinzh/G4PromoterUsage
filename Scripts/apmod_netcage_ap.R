### 

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(CAGEfightR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
#suppressPackageStartupMessages(library(ggthemes))
#suppressPackageStartupMessages(library(ggpie))
suppressPackageStartupMessages(library(Gviz))
#suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(UpSetR))

setwd("./G4AP/")

PROMO_DISTANCE <- 100

bsg <- BSgenome.Hsapiens.UCSC.hg19
gencode_path <- "./data/ref/gencode.v19.annotation.gtf"
txdb <- makeTxDbFromGFF(file = gencode_path, 
                        format = "gtf", 
                        organism = "Homo sapiens", 
                        taxonomyId = 9606)
odb <- org.Hs.eg.db

k562_G4 <- fread("./data/G4/cellline/K562_hg19.bed", 
                 sep = "\t", header = FALSE) %>% data.frame %>% 
  filter(V1 %in% paste0("chr", c(1:22, "X")))
hepg2_G4 <- fread("./data/G4/cellline/HepG2_hg19.bed", 
                  sep = "\t", header = FALSE) %>% data.frame %>% 
  filter(V1 %in% paste0("chr", c(1:22, "X")))
k562_G4_gr <- k562_G4
k562_G4_gr$score <- 1
colnames(k562_G4_gr)[1:3] <- c("chr", "start", "end")
k562_G4_gr <- makeGRangesFromDataFrame(k562_G4_gr,
                                       keep.extra.columns = TRUE,
                                       ignore.strand = TRUE)
hepg2_G4_gr <- hepg2_G4
hepg2_G4_gr$score <- 1
colnames(hepg2_G4_gr)[1:3] <- c("chr", "start", "end")
hepg2_G4_gr <- makeGRangesFromDataFrame(hepg2_G4_gr,
                                        keep.extra.columns = TRUE,
                                        ignore.strand = TRUE)

# create cage information object
cage_data <- bind_rows(
data.frame(Class = "K562", Name = paste0("kr", 1:6), 
           BigWigPlus = c("GSM3318241_CNhi10920_biologicalRep1-K562-NETCAGE-0_5M_CAC.ctss_positive.bw", 
                          "GSM3318242_CNhi10920_biologicalRep1-K562-NETCAGE-1M_AGT.ctss_positive.bw",
                          "GSM3318243_CNhi10920_biologicalRep1-K562-NETCAGE-2M_GCG.ctss_positive.bw",
                          "GSM3318244_CNhi10920_biologicalRep2-K562-NETCAGE-0_5M_TAC.ctss_positive.bw",
                          "GSM3318245_CNhi10920_biologicalRep2-K562-NETCAGE-1M_ACG.ctss_positive.bw",
                          "GSM3318246_CNhi10920_biologicalRep2-K562-NETCAGE-2M_GCT.ctss_positive.bw"), 
           BigWigMinus = c("GSM3318241_CNhi10920_biologicalRep1-K562-NETCAGE-0_5M_CAC.ctss_negative.bw",
                           "GSM3318242_CNhi10920_biologicalRep1-K562-NETCAGE-1M_AGT.ctss_negative.bw",
                           "GSM3318243_CNhi10920_biologicalRep1-K562-NETCAGE-2M_GCG.ctss_negative.bw",
                           "GSM3318244_CNhi10920_biologicalRep2-K562-NETCAGE-0_5M_TAC.ctss_negative.bw",
                           "GSM3318245_CNhi10920_biologicalRep2-K562-NETCAGE-1M_ACG.ctss_negative.bw",
                           "GSM3318246_CNhi10920_biologicalRep2-K562-NETCAGE-2M_GCT.ctss_negative.bw")),
data.frame(Class = "HepG2", Name = paste0("hr", 1:6), 
           BigWigPlus = c("GSM3318233_CNhi10919_biologicalRep1-HepG2-NETCAGE-0_5_CAC.ctss_positive.bw",
                          "GSM3318234_CNhi10919_biologicalRep1-HepG2-NETCAGE-1M_AGT.ctss_positive.bw",
                          "GSM3318235_CNhi10919_biologicalRep1-HepG2-NETCAGE-2M_GCG.ctss_positive.bw",
                          "GSM3318236_CNhi10919_biologicalRep2-HepG2-NETCAGE-0_5M_TAC.ctss_positive.bw",
                          "GSM3318237_CNhi10919_biologicalRep2-HepG2-NETCAGE-1M_ACG.ctss_positive.bw",
                          "GSM3318238_CNhi10919_biologicalRep2-HepG2-NETCAGE-2M_GCT.ctss_positive.bw"), 
           BigWigMinus = c("GSM3318233_CNhi10919_biologicalRep1-HepG2-NETCAGE-0_5_CAC.ctss_negative.bw",
                           "GSM3318234_CNhi10919_biologicalRep1-HepG2-NETCAGE-1M_AGT.ctss_negative.bw",
                           "GSM3318235_CNhi10919_biologicalRep1-HepG2-NETCAGE-2M_GCG.ctss_negative.bw",
                           "GSM3318236_CNhi10919_biologicalRep2-HepG2-NETCAGE-0_5M_TAC.ctss_negative.bw",
                           "GSM3318237_CNhi10919_biologicalRep2-HepG2-NETCAGE-1M_ACG.ctss_negative.bw",
                           "GSM3318238_CNhi10919_biologicalRep2-HepG2-NETCAGE-2M_GCT.ctss_negative.bw")))
cage_data$Class = factor(cage_data$Class, level = c("K562", "HepG2"))
rownames(cage_data) <- cage_data$Name

# create bw object
bw_plus <- c(list.files("./data/cage/K562/", "*positive.bw$", full.names = TRUE), 
             list.files("./data/cage/HepG2/", "*positive.bw$", full.names = TRUE))
bw_minus <- c(list.files("./data/cage/K562/", "*negative.bw$", full.names = TRUE), 
              list.files("./data/cage/HepG2/", "*negative.bw$", full.names = TRUE))
bw_plus <- BigWigFileList(bw_plus)
bw_minus <- BigWigFileList(bw_minus)
names(bw_plus) <- names(bw_minus) <- cage_data$Name

CTSSs <- quantifyCTSSs(plusStrand = bw_plus,
                       minusStrand = bw_minus,
                       genome = seqinfo(bsg),
                       design = cage_data)
CTSSs <- CTSSs %>% calcTPM() %>% calcPooled()
TCs <- quickTSSs(CTSSs)
Prom <- TCs %>% calcTPM() %>%
        subsetBySupport(inputAssay = "TPM", 
                        unexpressed = 1, 
                        minSamples = 5) 
Prom <- assignTxID(Prom, txModels = txdb, swap = "thick")
Prom <- assignTxType(Prom, txModels = txdb, swap = "thick")

# detect potential enhancers
BCs <- quickEnhancers(CTSSs)
BCs <- subsetBySupport(BCs, 
                       inputAssay = "counts", 
                       unexpressed = 0, 
                       minSamples = 5) 

# remove Prom that overlapped with enhancers
overlaps <- findOverlaps(rowRanges(Prom), rowRanges(BCs))
Prom <- Prom[setdiff(seq_along(rowRanges(Prom)), queryHits(overlaps)), ]
Prom$totalTags <- NULL
Prom <- calcTPM(Prom) %>% calcPooled()

Prom <- assignGeneID(Prom, geneModels = txdb, upstream = 0, downstream = 1)

anno_prom <- promoters(transcripts(txdb, columns = c("TXSTART")), 
                      upstream = 0, downstream = 1)
overlapping <- findOverlaps(anno_prom, Prom, type = "within")
Prom <- Prom[unique(subjectHits(overlapping)), ]
                      
Prom <- Prom %>% subset(!is.na(geneID)) %>% 
  subsetByComposition(inputAssay = "counts", 
                      genes = "geneID", 
                      unexpressed = 0.1, 
                      minSamples = 5)

rowData(Prom)$symbol <- mapIds(odb, 
                               keys = str_split( 
                                 rowData(Prom)$geneID, "\\.", 
                                 simplify = TRUE)[, 1],
                               column = "SYMBOL", 
                               keytype = "ENSEMBL")
Prom <- Prom %>% subset(!is.na(symbol))

coding_genes <- select(org.Hs.eg.db,
                       keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"), 
                       columns = c("GENETYPE", "CHR"), 
                       keytype = "ENSEMBL") %>% 
  subset(GENETYPE == "protein-coding" & CHR %in% c(1:22, "X")) %>%
  .$ENSEMBL %>% unique()

## protein coding genes only!!!
Prom <- Prom %>% 
  subset(str_split(geneID, "\\.", simplify = TRUE)[, 1] %in% coding_genes)

# plot the number of genes with different numbers of Prom
pdf("./figure/ap/Multi_Prom_freq.pdf", width = 6, height = 5)
Prom %>% rowData() %>% as.data.frame() %>% as_tibble() %>% 
  dplyr::count(geneID) %>% 
  ggplot(aes(x = n, fill = n >= 2)) + 
  geom_bar(alpha = 1, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("#1C5B9B", "#A81C32"), 
                    labels = c("Single-Promoter", "Multi-Promoter")) +
  labs(x = "Number of Prom Used per Gene", y = "Frequency", 
       fill = "Promoter Category") +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size = 14, face = "bold"), 
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.background = element_blank(),
    legend.key = element_blank()
  )
dev.off()

########################
# differential used Prom
########################

selected_Prom <- Prom %>% rowData() %>% 
  subset(select = c(score, txType, geneID, symbol)) %>% as.data.frame()

dge <- DGEList(counts = assay(Prom, "counts"), genes = selected_Prom)

dge <- calcNormFactors(dge)
mod <- model.matrix(~ Class, data = colData(Prom))

disp <- estimateDisp(dge, design = mod, tagwise = FALSE)
QLfit <- glmQLFit(disp, design = mod, robust = TRUE)

du_prom <- diffSpliceDGE(QLfit, geneid = "geneID")
du_prom <- topSpliceDGE(du_prom, test = "exon", n = 1000000)

du_prom_gr <- GRanges(
  seqnames = str_split(rownames(du_prom), ":", simplify = TRUE)[, 1],
  ranges = IRanges(
    start = as.numeric(str_split(rownames(du_prom), "(:|-)", simplify = TRUE)[, 2]),
    end = as.numeric(str_split(rownames(du_prom), "(;|-)", simplify = TRUE)[, 2])
  ),
  strand = str_split(rownames(du_prom), ";", simplify = TRUE)[, 2],
  promInfo = rownames(du_prom),
  geneID =  du_prom$geneID,
  symbol = du_prom$symbol,
  logFC = du_prom$logFC,
  P.Value = du_prom$P.Value,
  FDR = du_prom$FDR
)

nearest_hits <- distanceToNearest(du_prom_gr, k562_G4_gr)
distances <- mcols(nearest_hits)$distance
du_prom_gr$K562_G4_proximity <- ifelse(distances < PROMO_DISTANCE, 1, 0)
du_prom_gr$K562_G4_distance <- distances

nearest_hits <- distanceToNearest(du_prom_gr, hepg2_G4_gr)
distances <- mcols(nearest_hits)$distance
du_prom_gr$HepG2_G4_proximity <- ifelse(distances < PROMO_DISTANCE, 1, 0)
du_prom_gr$HepG2_G4_distance <- distances

du_prom_sig <- 
  du_prom_gr[abs(mcols(du_prom_gr)$logFC) > 1 & mcols(du_prom_gr)$FDR < 0.05]

fwrite(data.frame(du_prom_sig)[, 1:3], "./data/du_prom_sig.bed", sep = "\t", 
       col.names = FALSE, row.names = FALSE, quote = FALSE)

distance_data <- as.data.frame(du_prom_sig) %>%
  dplyr::select(K562_G4_distance, HepG2_G4_distance) %>%
  tidyr::pivot_longer(cols = c(K562_G4_distance, HepG2_G4_distance), 
                      names_to = "CellLine", 
                      values_to = "Distance")
distance_data$Distance_log <- log10(distance_data$Distance + 1)

pdf("./figure/ap/Distance_prom_G4.pdf", width = 8, height = 5)
ggplot(distance_data, aes(x = Distance_log, fill = CellLine)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) + 
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  labs(
    title = "Distance from Alternative Promoter to Nearest G4 (Log-scaled)", 
    x = "Distance (log10)", 
    y = "Density",
    fill = "Cell Line"
  )
dev.off()

sigprom_k562_g4 <- normalizeToMatrix(k562_G4_gr, du_prom_sig, 
                                    value_column = "score", extend = 1000, 
                                    mean_mode = "coverage", w = 10)

sigprom_hepg2_g4 <- normalizeToMatrix(hepg2_G4_gr, du_prom_sig, 
                                    value_column = "score", extend = 1000, 
                                    mean_mode = "coverage", w = 10)

pdf("./figure/ap/Distribution_K562_G4.pdf", width = 2, height = 10)
EnrichedHeatmap(sigprom_k562_g4, name = "K562 G4")
dev.off()
pdf("./figure/ap/Distribution_HepG2_G4.pdf", width = 2, height = 10)
EnrichedHeatmap(sigprom_hepg2_g4, name = "HepG2 G4")
dev.off()

# fold change in 
state <- du_prom_sig$K562_G4_proximity + du_prom_sig$HepG2_G4_proximity

G4_state_fc <- bind_rows(
  #G4-Shared Promoters G4-Unique / Unshared Promoters
  data.frame(logfc = abs(du_prom_sig[state == 1]$logFC), type = "G4-Differential AP"),
  data.frame(logfc = abs(du_prom_sig[state == 2]$logFC), type = "G4-Shared AP"))

pdf("./figure/ap/AP_G4_state_logfc.pdf", width = 4, height = 4)
ggplot(G4_state_fc, aes(x = type, y = abs(logfc), fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, color = "black", size = 1) +  
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.4, size = 3) +  
  scale_fill_manual(values = c("G4-Differential AP" = "#199994", "G4-Shared AP" = "#c9aa82")) +  
  scale_x_discrete(limits = c("G4-Differential AP", "G4-Shared AP")) +  
  labs(x = "", y = "Absolute log2 Fold Change") +  
  theme_classic(base_size = 14) +  
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )
dev.off()

logfc_G4_single <- bind_rows(
  data.frame(logFC = mcols(du_prom_sig)[du_prom_sig$K562_G4_proximity == 1 & du_prom_sig$HepG2_G4_proximity == 0, "logFC", drop = FALSE],
             group = "K562-specific"),
  data.frame(logFC = mcols(du_prom_sig)[du_prom_sig$K562_G4_proximity == 0 & du_prom_sig$HepG2_G4_proximity == 1, "logFC", drop = FALSE],
             group = "HepG2-specific")
)


logfc_range <- range(logfc_G4_single$logFC, na.rm = TRUE)

p1 <- ggplot(logfc_G4_single, aes(x = logFC, fill = group, color = group)) +
  geom_histogram(aes(y = after_stat(count)), bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("K562-specific" = "#7570b3", "HepG2-specific" = "#e5a024")) +
  scale_color_manual(values = c("K562-specific" = "#7570b3", "HepG2-specific" = "#e5a024")) +
  coord_cartesian(xlim = logfc_range) + 
  theme_classic(base_size = 18) +
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.position = c(0.8, 0.8),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_blank()
  ) +
  labs(y = "Count", fill = "G4-associated AP Group", 
       color = "G4-associated AP Group")

p2 <- ggplot(logfc_G4_single, aes(x = logFC, y = -1, color = group)) +
  geom_quasirandom(size = 2.5, alpha = 0.8, width = 0.3, method = "quasirandom") +  
  scale_color_manual(values = c("K562-specific" = "#7570b3", 
                                "HepG2-specific" = "#e5a024")) +
  coord_cartesian(xlim = logfc_range) +
  theme_classic(base_size = 18) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Log2 Fold Change")

pdf("./figure/ap/G4_associated_specific_log2fc.pdf", width = 8, height = 5)
p1 / p2 + plot_layout(heights = c(3, 1))  
dev.off()

# hepg2 up-regulated == k562 down-regulated
hepg2_up_genes <- data.frame(du_prom_sig) %>% filter(logFC > 0) %>% 
  pull(geneID) %>% unique()
# hepg2 down-regulated == k562 up-regulated
hepg2_down_genes <- data.frame(du_prom_sig) %>% filter(logFC < 0) %>% 
  pull(geneID) %>% unique()

up_down_df <- data.frame(
  Gene = unique(c(hepg2_up_genes, hepg2_down_genes))
)

up_down_df$`Upregulated in HepG2 \n(Downregulated in K562)` <- 
  up_down_df$Gene %in% hepg2_up_genes
up_down_df$`Downregulated in HepG2 \n(Upregulated in K562)` <- 
  up_down_df$Gene %in% hepg2_down_genes

up_down_df[, 2:3] <- lapply(up_down_df[, 2:3], as.numeric)

pdf("./figure/ap/Sig_Gene_upset.pdf", width = 6, height = 6, onefile = F)
upset(
  up_down_df, 
  sets = c("Upregulated in HepG2 \n(Downregulated in K562)",
           "Downregulated in HepG2 \n(Upregulated in K562)"), 
  keep.order = TRUE,
  sets.bar.color = c("#E41A1C", "#377EB8"),
  main.bar.color = "#8b2251",
  matrix.color = "#914900",
  mainbar.y.label = "Number of Genes",
  text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 2),
  point.size = 3.5,
  line.size = 1.2,
  shade.color = c("#8B7E66")
  )
dev.off()

seqs <- getSeq(bsg, du_prom_sig)
names(seqs) <- mcols(du_prom_sig)$promInfo
writeXStringSet(seqs, "./data/input_sequences.fasta")

labels <- c("G4-Associated AP", "Other AP")
values <- c(length(du_prom_sig) - table(du_prom_sig$K562_G4_proximity+du_prom_sig$HepG2_G4_proximity)[1],
            table(du_prom_sig$K562_G4_proximity+du_prom_sig$HepG2_G4_proximity)[1])
percentages <- round(values / sum(values) * 100, 1)
labels <- paste(labels, "\n", percentages, "%") 
pdf("./figure/ap/G4_assocaited_AP_pie.pdf", width = 5, height = 5)
pie(values, labels = labels, col =  c("#576fa0", "#a7b9d7"), cex = 1)
dev.off()

pdf("./figure/ap/AP_count.pdf", width = 3.5, height = 4.5)
barplot(c(table(du_prom_sig$logFC > 0)[1],
          table(du_prom_sig$logFC > 0)[2]),
        names.arg = c("Up-regulated\nin K562", 
                      "Up-regulated\nin HepG2"), 
        col = c("#a8786e", "#a2a2a2"), 
        ylab = "Count", 
        border = "black", 
        space = 0.5)
dev.off()

plotAP <- function(f_path = NULL, transcripts = NULL, plot_region = plot_region) {
  
  #stop("Please run this function in Rstudio")
  
  axis_track <- GenomeAxisTrack(col = "black", fontcolor = "black")
  tx_track <- GeneRegionTrack(txdb, 
                              name = "x", 
                              col = NA,
                              #fill = "#75c0c1", 
                              shape = c("smallArrow", "box"),
                              showId = FALSE,
                              geneSymbol = FALSE,
                              stacking = "full",
                              fill = "#f99d1e")
  if (length(transcripts) > 0) {
    tx_track <- tx_track[transcript(tx_track)%in%transcripts]
  }
  
  cluster_track <- Prom %>%
    subsetByOverlaps(plot_region) %>%
    trackClusters(name = "x", col = NA, showId = FALSE,
                  plusColor = "#02a2d9", minusColor = "#bd1a36")
  
  track1 <- CTSSs %>%
    subset(select = Class == "K562") %>%
    calcPooled() %>%
    subsetByOverlaps(plot_region) %>%
    trackCTSS(name = "x", plusColor = "#02a2d9", minusColor = "#bd1a36")
  
  track2 <- CTSSs %>%
    subset(select = Class == "HepG2") %>%
    calcPooled() %>%
    subsetByOverlaps(plot_region) %>%
    trackCTSS(name = "x", plusColor = "#02a2d9", minusColor = "#bd1a36")
  
  k562_G4_gr_k <- AnnotationTrack(k562_G4_gr, name = "x", col = NA, fill = "#9567bd")
  hepg2_G4_gr_k <- AnnotationTrack(hepg2_G4_gr, name = "x", col = NA, fill = "#6ba3d6")
  
  png(f_path, width = 5, height = 4, units = "in", res = 300)
  print(plotTracks(list(axis_track,
                        tx_track,
                        cluster_track,
                        track1,
                        track2,
                        k562_G4_gr_k,
                        hepg2_G4_gr_k),
                   from = start(plot_region),
                   to = end(plot_region),
                   chromosome = as.character(seqnames(plot_region)),
                   scale = 0.2,
                   #ylim = ylim, 
                   background.title = "steelblue", fontcolor = "black"))
  dev.off()
}

gn = "PDE4DIP"
plot_region <- Prom %>%
  rowRanges() %>%
  subset(symbol == gn) %>%
  GenomicRanges::reduce(min.gapwidth = 1e6L) %>%
  unstrand() %>%
  add(12e3L)
plotAP(f_path = "./figure/ap/PDE4DIP.png", transcripts = c("ENST00000313431.9", 
                                                           "ENST00000369354.3"), 
       plot_region = plot_region)


gn = "WDR90"
plot_region <- Prom %>%
  rowRanges() %>%
  subset(symbol == gn) %>%
  GenomicRanges::reduce(min.gapwidth = 1e6L) %>%
  unstrand() %>%
  add(3e3L)
plotAP(f_path = "./figure/ap/WDR90.png", transcripts = c("ENST00000549024.1", 
                                               "ENST00000552943.1"), 
       plot_region = plot_region)
