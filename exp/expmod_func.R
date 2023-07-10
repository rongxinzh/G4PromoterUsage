### 

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(CAGEfightR))
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
  
  return(G4)
}

GetG4gr <- function(data = NULL) {
  data <- data[, 1:3]
  data$score <- 1
  colnames(data)[1:3] <- c("chr", "start", "end")
  data <- makeGRangesFromDataFrame(data,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE)

  return(data)
}

GetBWObj <- function(cage.data = NULL, bwplist = NULL, bwmlist = NULL) {
  
  bw.plus <- bwplist
  bw.minus <- bwmlist
  bw.plus <- BigWigFileList(bw.plus)
  bw.minus <- BigWigFileList(bw.minus)
  names(bw.plus) <- names(bw.minus) <- cage.data$Name

  return(list(plus = bw.plus, minus = bw.minus))
}

quantifyTSSEh <- function(ctss = NULL, txdb = NULL) {
 
  tcs <- quickTSSs(ctss)
  tss <- tcs %>% calcTPM() %>%
         subsetBySupport(inputAssay = "TPM", 
                         unexpressed = 1,
                         minSamples = 1)

  bcs <- quickEnhancers(ctss)
  bcs <- subsetBySupport(bcs, 
                         inputAssay = "counts", 
                         unexpressed = 0, 
                         minSamples = 1)

  tss <- assignTxType(tss, txModels = txdb, swap = "thick")
  bcs <- assignTxType(bcs, txModels = txdb, swap = "thick")
  enhancers <- subset(bcs, txType %in% c("intergenic", "intron"))

  tss$totalTags <- NULL
  enhancers$totalTags <- NULL
  rowData(tss)$clusterType <- "TSS"
  rowData(enhancers)$clusterType <- "Enhancer"

  rse <- combineClusters(object1 = tss, 
                           object2 = enhancers, 
                           removeIfOverlapping = "object1")
  rowRanges(rse)$clusterType <- factor(rowRanges(rse)$clusterType, levels=c("TSS", "Enhancer"))
  rse <- calcTPM(rse) %>% calcPooled()

  return(rse)

}

extendPrmtr <- function(RSE = NULL) {
  RSE.tss <- rowRanges(RSE) %>% subset(clusterType == "TSS")
  start(RSE.tss) <- ifelse(strand(RSE.tss) == "+", start(RSE.tss) - 300, start(RSE.tss))
  end(RSE.tss) <- ifelse(strand(RSE.tss) == "+", end(RSE.tss), end(RSE.tss) + 300)
  score(RSE.tss) <- score(RSE.tss) / ncol(RSE)

  return(RSE.tss)

}

GetBW <- function(tss.info = NULL, bw.file = NULL) {
  
  tmp.file1 <- tempfile()
  tmp.data <- tss.info[, 1:4]
  fwrite(tmp.data, tmp.file1, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  tmp.file2 <- tempfile()
  command <- system(paste0("bigWigAverageOverBed ", bw.file, " ", tmp.file1, " ", tmp.file2))
  sig <- fread(tmp.file2, sep = '\t', header = FALSE) %>% data.frame
  sig <- sig[, c(1, 6)]
  colnames(sig) <- c("newId", "epi.value")

  unlink(tmp.file1)
  unlink(tmp.file2)

  return(sig)
}
