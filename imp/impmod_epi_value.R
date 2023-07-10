
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))

load(file = "../exp/RSE.tss.Rdata")

set.seed(1)

bw.path <- "../data/roadmap/"
bw.files <- list.files(bw.path, recursive=TRUE)
bw.files <- data.frame(file.name = bw.files, epi = str_split(bw.files, "(-|\\.)", simplify = TRUE)[, 2], code = str_split(bw.files, "(-|\\.)", simplify = TRUE)[, 1])

RSE.k.tss.dt <- data.frame(RSE.k.tss)[, c(1:3, 6)]
RSE.k.tss.dt$newId <- paste0("id_", 1:dim(RSE.k.tss.dt)[1])
RSE.k.tss.dt <- RSE.k.tss.dt %>% dplyr::select(seqnames:end, newId, score)

RSE.h.tss.dt <- data.frame(RSE.h.tss)[, c(1:3, 6)]
RSE.h.tss.dt$newId <- paste0("id_", 1:dim(RSE.h.tss.dt)[1])
RSE.h.tss.dt <- RSE.h.tss.dt %>% dplyr::select(seqnames:end, newId, score)

GetBW <- function(tss.info = NULL, bw.file = NULL) {
  
  tmp.file1 <- tempfile()
  tmp.data <- tss.info[, 1:4]
  fwrite(tmp.data, tmp.file1, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  tmp.file2 <- tempfile()
  command <- system(paste0("bigWigAverageOverBed ../data/roadmap/", bw.file, " ", tmp.file1, " ", tmp.file2))
  sig <- fread(tmp.file2, sep = '\t', header = FALSE) %>% data.frame
  sig <- sig[, c(1, 6)]
  colnames(sig) <- c("newId", "epi.value")

  unlink(tmp.file1)
  unlink(tmp.file2)

  return(sig)
}

mat.k <- NULL
mat.h <- NULL

for (tmp.epi in unique(bw.files$epi)) {
  
  print(tmp.epi)
  tmp.bw <- bw.files %>% filter(epi == tmp.epi) 
 
  sig.k <- GetBW(RSE.k.tss.dt, tmp.bw[tmp.bw$code == "E123", 1])
  sig.k <- RSE.k.tss.dt %>% left_join(sig.k, by = "newId")

  sig.h <- GetBW(RSE.h.tss.dt, tmp.bw[tmp.bw$code == "E118", 1])
  sig.h <- RSE.h.tss.dt %>% left_join(sig.h, by = "newId")

  mat.k <- bind_cols(mat.k, data.frame(sig.k$epi.value))
  mat.h <- bind_cols(mat.h, data.frame(sig.h$epi.value))
}

colnames(mat.k) <- colnames(mat.h) <- unique(bw.files$epi)
mat.k$score <- RSE.k.tss.dt$score
mat.h$score <- RSE.h.tss.dt$score
mat.k.unscaled <- mat.k
mat.h.unscaled <- mat.h

mat.k[, 1:(dim(mat.k)[2]-1)] <- scale(mat.k[, 1:(dim(mat.k)[2]-1)])
mat.h[, 1:(dim(mat.h)[2]-1)] <- scale(mat.h[, 1:(dim(mat.h)[2]-1)])

mat.k$score <- ifelse(mat.k$score > quantile(mat.k$score)[4], "c1", ifelse(mat.k$score < quantile(mat.k$score)[2], "c0", "c2"))
mat.k <- mat.k[mat.k$score != "c2",]
mat.k$score <- as.factor(mat.k$score)

mat.h$score <- ifelse(mat.h$score > quantile(mat.h$score)[4], "c1", ifelse(mat.h$score < quantile(mat.h$score)[2], "c0", "c2"))
mat.h <- mat.h[mat.h$score != "c2",]
mat.h$score <- as.factor(mat.h$score)

save(list=c("mat.k", "mat.h", "mat.k.unscaled", "mat.h.unscaled"), file = "./mat_for_ML.Rdata")


