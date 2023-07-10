### 

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(introdataviz))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(ggpie))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(ggforce))


PlotTSSWidth <- function(f.path = NULL, data = NULL, col = NULL) {

  fig <- ggdensity(data, x = "Length", y = "density",
                   color = col, fill = col, alpha = 0.3, add = "mean", rug = TRUE)
  ggsave(f.path, fig, width = 5, height = 3.5)

}

PlotG4Prop <- function(f.path = NULL, RSE.k.tssg4 = NULL, RSE.h.tssg4 = NULL) {

  fig.k <- ggdonut(data = RSE.k.tssg4, group_key = "label", count_type = "full",
                 label_info = "all", label_type = "horizon", donut.label.size = 6,
                 label_size = 6, label_pos = "in", label_threshold = 10, fill_color = c("#a2a2a2", "#d82626"))

  fig.h <- ggdonut(data = RSE.h.tssg4, group_key = "label", count_type = "full",
                 label_info = "all", label_type = "horizon", donut.label.size = 6,
                 label_size = 6, label_pos = "in", label_threshold = 10, fill_color = c("#a2a2a2", "#d82626"))
  ggsave(f.path, ggarrange(fig.k, fig.h), width = 16, height = 5)

} 

PlotSpineG4Prop <- function(f.path = NULL, RSE.tssbar = NULL, col = NULL) {
  
  pdf(f.path)
  par(las = 3)
  spineplot(RSE.tssbar, col = col, xlab = '', ylab = 'Group')
  dev.off()

}

PlotGroupedG4TssExp <- function(f.path = NULL, RSE.tss.exp = NULL) {

  fig <- ggviolin(RSE.tss.exp, "cell", "exp", fill = "group",
           palette = c("#ff9e49", "#ababab"),
           add = "boxplot")
  fig <- facet(fig + theme_bw(), facet.by = "cell", scales = "free")


  ggsave(f.path, fig, width = 10, height = 4)

}


PlotTSSFreq <- function(f.path = NULL, RSE.k.tssp, RSE.h.tssp, fillcol = NULL) {

  k.tssp <-  ggplot(RSE.k.tssp, aes(x = Var1, y = Freq, fill = fillcol)) +
             geom_col() +
             geom_text(
             aes(label = perc), vjust = 0) +
             scale_fill_identity(guide = "none") +
             theme_minimal() + theme(axis.title.x = element_blank())
  h.tssp <-  ggplot(RSE.h.tssp, aes(x = Var1, y = Freq, fill = fillcol)) +
             geom_col() +
             geom_text(
             aes(label = perc), vjust = 0) +
             scale_fill_identity(guide = "none") +
             theme_minimal() + theme(axis.title.x = element_blank())

  ggsave(f.path, ggarrange(k.tssp, h.tssp), width = 12, height = 3)

}

PlotLenIQR <- function(f.path = NULL, high.tss = NULL, col = NULL) {
  
  ggsave(f.path,
  high.tss %>% rowData() %>% data.frame() %>% 
    ggplot(aes(x = IQR)) +
    geom_histogram(binwidth = 1, 
                   fill = col, 
                   alpha = 0.9) +
    geom_vline(xintercept = 10, 
               linetype = "dashed", 
               alpha = 0.75, 
               color = "black") +
    facet_zoom(xlim = c(0,100)) +
    labs(x = "10-90% IQR", 
         y = "Frequency"),
    width = 6, height = 3)
  
}

PlotShapeProp <- function(f.path = NULL, data1 = NULL, data2 = NULL) {

  fig.k <- ggdonut(data = data1, group_key = "shape", count_type = "full",
                   label_info = "all", label_type = "horizon", donut.label.size = 6,
                   label_size = 6, label_pos = "in", label_threshold = 10, c("#d82626", "#a2a2a2"))

  fig.h <- ggdonut(data = data2, group_key = "shape", count_type = "full",
                   label_info = "all", label_type = "horizon", donut.label.size = 6,
                   label_size = 6, label_pos = "in", label_threshold = 10, c("#d82626", "#a2a2a2"))
  ggsave(f.path, ggarrange(fig.k, fig.h), width = 16, height = 5)

}

PlotTSSMotif <- function(f.path = NULL, prmtr.seqs = NULL, high.tss = NULL) {
  
  ggsave(f.path,
  prmtr.seqs %>% as.character %>% split(rowData(high.tss)$shape) %>% 
           ggseqlogo(data = ., ncol = 1, seq_type = "dna") +
           theme_logo() +
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()))

}

PlotTSSExp <- function(f.path = NULL, txdb = NULL, RSE = NULL, CTSSs = NULL, G4.gr = NULL) {
  
  #stop("Please use rstudio for figure visualization")
  axis.track <- GenomeAxisTrack(col = "black", fontcolor = "black")
  tx.track <- GeneRegionTrack(txdb, 
                              name = "Gene Annotation", 
                              col = NA,
                              fill = "#f99d1e",
                              shape = c("smallArrow", "box"),
                              showId = TRUE,
                              geneSymbol = TRUE)

  plot.region <- RSE %>% rowRanges() %>% subset(clusterType == "TSS") %>% .[35] %>% add(500) %>% unstrand()
  ctss.track <- CTSSs %>% rowRanges() %>% subsetByOverlaps(plot.region) %>% trackCTSS(name = "CTSSs", plusColor = "#02a2d9", minusColor = "#bd1a36")
  cluster.track <- RSE %>% subsetByOverlaps(plot.region) %>% trackClusters(name = "Clusters", col = NA, plusColor = "#02a2d9", 
                                                                           minusColor = "#bd1a36", showId = FALSE)
  G4.gr.track <- AnnotationTrack(G4.gr, name = "G4 peak", col = NA, fill = "#9567bd")
  png(f.path, width = 5, height = 4, units = "in", res = 300)
  plotTracks(list(axis.track, 
                  ctss.track,
                  cluster.track,
                  G4.gr.track,
                  tx.track),
            from = start(plot.region), 
            to = end(plot.region), 
            chromosome = as.character(seqnames(plot.region)),
            scale = 0.5,
            background.title = "steelblue", fontcolor = "black")
  dev.off()
}

