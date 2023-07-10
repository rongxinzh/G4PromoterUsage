
set.seed(1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(GenomicRanges))

GetTFBS <- function(cl = NULL) {
  path <- switch(
    tolower(cl),
    "k562" = "../data/tfbs/remap2022_nr_macs2_hg19_v1_0_K-562_processed.bed",
    "hepg2" = "../data/tfbs/remap2022_nr_macs2_hg19_v1_0_Hep-G2_processed.bed",
    stop("Please provide a valid cell line name.")
    )
  tfbs <- fread(path, sep = "\t", header = FALSE) %>% data.frame
  tfbs <- tfbs %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
  
  return(tfbs)
}

GetG4 <- function(cl = NULL) {
  path <- switch(
    tolower(cl),
    "k562" = "../data/G4/cellline/K562_hg19.bed",
    "hepg2" = "../data/G4/cellline/HepG2_hg19.bed",
    stop("Please provide a valid cell line name.")
    )
  G4 <- fread(path, sep = "\t", header = FALSE) %>% data.frame
  G4 <- G4 %>% filter(V1 %in% paste0("chr", c(1:22, "X")))
  
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

TSSgr <- function(data = NULL) {

  data.s <- str_split(data$V10, "[:\\-;]", simplify = TRUE)[, 1:3] %>% data.frame
  data.s[, 2] <- as.numeric(data.s[, 2])
  data.s[, 3] <- as.numeric(data.s[, 3])
  data.s[, 2] <- ceiling((data.s[, 2] + data.s[, 3]) / 2)
  data.s[, 3] <- data.s[, 2]
  data.s <- bind_cols(data.s, data$V5, data$label)
  colnames(data.s) <- c("chr", "start", "end", "strand", "label")
  data.s$score <- 1
  data.s.gr <- makeGRangesFromDataFrame(data.s,
                                            keep.extra.columns = TRUE,
                                            ignore.strand = FALSE)
  return(data.s.gr)

}
