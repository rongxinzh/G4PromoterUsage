
set.seed(1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(bedtoolsr))

GetG4 <- function(cl = NULL) {
  path <- switch(
    tolower(cl),
    "k562" = "../data/G4/cellline/K562_hg19.bed",
    "hepg2" = "../data/G4/cellline/HepG2_hg19.bed",
    stop("Please provide a valid cell line name.")
    )
  G4 <- fread(path, sep = "\t", header = FALSE) %>% data.frame
  G4 <- G4 %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
  G4[, 2] <- ceiling((G4[, 2] + G4[, 3]) / 2)
  G4[, 3] <- G4[, 2] + 1

  return(G4)
}

Getgr <- function(data = NULL) {
  data <- data[, 1:3]
  data$score <- 1
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- makeGRangesFromDataFrame(data,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE)

  return(data)
}

k562.G4 <- GetG4("k562")
hepg2.G4 <- GetG4("hepg2")
k562.G4.gr <- Getgr(k562.G4)
hepg2.G4.gr <- Getgr(hepg2.G4)

k.property <- list.files("../data/roadmap_gpeak/", pattern = "E123-*", full.names = TRUE)
h.property <- list.files("../data/roadmap_gpeak/", pattern = "E118-*", full.names = TRUE)

mat.k <- NULL
mat.h <- NULL

for (f.path in k.property) {
  message(f.path)

  tmp.prop <- str_split(f.path, "\\.imputed|-", simplify = TRUE)[, 2]
  tmp.gr <- rtracklayer::import(f.path, format = "bed")
  score(tmp.gr) <- 1

  tmp.mat.k <- normalizeToMatrix(tmp.gr, k562.G4.gr, 
    value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

  tmp.mat.h <- normalizeToMatrix(tmp.gr, hepg2.G4.gr, 
    value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

  mat.k <- bind_rows(mat.k, data.frame(value = colMeans(data.frame(tmp.mat.k)), 
                position = 1:ncol(data.frame(tmp.mat.k)),
                epi = tmp.prop, group = "G4"))

  mat.h <- bind_rows(mat.h, data.frame(value = colMeans(data.frame(tmp.mat.h)), 
                position = 1:ncol(data.frame(tmp.mat.h)),
                epi = tmp.prop, group = "G4"))

} 

p.k <- ggline(mat.k, "position", "value", color = "group", facet.by = "epi", ncol = 8, scales = "free_y", palette = c("#60b8d3"), plot_type = "l", size = 1) + 
rremove("xlab") + ylab("Relative density") + scale_x_continuous(breaks = c(1, 101, 201), labels=c("-1kb", "G4 center", "1kb")) + rremove("legend")
p.h <- ggline(mat.h, "position", "value", color = "group", facet.by = "epi", ncol = 8, scales = "free_y", palette = c("#60b8d3"), plot_type = "l", size = 1) + 
rremove("xlab") + ylab("Relative density") + scale_x_continuous(breaks = c(1, 101, 201), labels=c("-1kb", "G4 center", "1kb")) + rremove("legend")

ggsave("../figure/cage/imp/K562_G4_genomicproperty.pdf", p.k, width = 18, height = 6) 
ggsave("../figure/cage/imp/HepG2_G4_genomicproperty.pdf", p.h, width = 18, height = 6) 


