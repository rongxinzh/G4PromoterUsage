
set.seed(1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(EnrichedHeatmap))

source("./tfmod_func.R")
load(file = "../exp/RSE.tss.Rdata")

RSE.k.tss.dt <- data.frame(RSE.k.tss)[, c(1:3, 6)]
RSE.k.tss.dt$newId <- paste0("id_", 1:dim(RSE.k.tss.dt)[1])
RSE.k.tss.dt <- RSE.k.tss.dt %>% dplyr::select(seqnames:end, newId, score)

RSE.h.tss.dt <- data.frame(RSE.h.tss)[, c(1:3, 6)]
RSE.h.tss.dt$newId <- paste0("id_", 1:dim(RSE.h.tss.dt)[1])
RSE.h.tss.dt <- RSE.h.tss.dt %>% dplyr::select(seqnames:end, newId, score)

tfbs.k <- GetTFBS("K562")
tfbs.h <- GetTFBS("HepG2")
tfbs.k.gr <- Getgr(tfbs.k)
tfbs.h.gr <- Getgr(tfbs.h)

k562.G4 <- GetG4("k562")
hepg2.G4 <- GetG4("hepg2")
k562.G4.gr <- Getgr(k562.G4)
hepg2.G4.gr <- Getgr(hepg2.G4)

RSE.k.tssg4 <- bt.intersect(a = data.frame(RSE.k.tss), b = k562.G4, wa = TRUE, c = TRUE) %>% unique
RSE.k.tssg4$label <- ifelse(RSE.k.tssg4[, ncol(RSE.k.tssg4)] == 0, "non G4-associated", "G4-associated")
RSE.h.tssg4 <- bt.intersect(a = data.frame(RSE.h.tss), b = hepg2.G4, wa = TRUE, c = TRUE) %>% unique
RSE.h.tssg4$label <- ifelse(RSE.h.tssg4[, ncol(RSE.h.tssg4)] == 0, "non G4-associated", "G4-associated")

RSE.k.tssg4.gr <- TSSgr(RSE.k.tssg4)
RSE.h.tssg4.gr <- TSSgr(RSE.h.tssg4)

RSE.k.tssg4Y.gr <- subset(RSE.k.tssg4.gr, label == "G4-associated")
RSE.k.tssg4N.gr <- subset(RSE.k.tssg4.gr, label == "non G4-associated")
RSE.h.tssg4Y.gr <- subset(RSE.h.tssg4.gr, label == "G4-associated")
RSE.h.tssg4N.gr <- subset(RSE.h.tssg4.gr, label == "non G4-associated")

mat.k.Y <- normalizeToMatrix(tfbs.k.gr, RSE.k.tssg4Y.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.k.N <- normalizeToMatrix(tfbs.k.gr, RSE.k.tssg4N.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

mat.h.Y <- normalizeToMatrix(tfbs.h.gr, RSE.h.tssg4Y.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
mat.h.N <- normalizeToMatrix(tfbs.h.gr, RSE.h.tssg4N.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

pdf("../figure/cage/tf/K562_G4_TF_TSS.pdf", width = 5, height = 5)
plot(colMeans(data.frame(mat.k.Y)), type="l", lwd=3, col = "orange", xlab = "", ylab = "Binding density of TFs", xaxt="n")
lines(colMeans(data.frame(mat.k.N)), lwd=3, col="steelblue")
axis(side=1, at=c(1,101,201),labels=c("-1kb","Promoter center","1kb"))
legend(1, 100, legend=c("G4-associated promoters", "Other promoters"),
       col=c("orange", "steelblue"), lty=1, cex=0.8, lwd = 2)
dev.off()

pdf("../figure/cage/tf/HepG2_G4_TF_TSS.pdf", width = 5, height = 5)
plot(colMeans(data.frame(mat.h.Y)), type="l", lwd=3, col = "orange", xlab = "", ylab = "Binding density of TFs", xaxt="n")
lines(colMeans(data.frame(mat.h.N)), lwd=3, col="steelblue")
axis(side=1, at=c(1,101,201),labels=c("-1kb","Promoter center","1kb"))
legend(1, 120, legend=c("G4-associated promoters", "Other promoters"),
       col=c("orange", "steelblue"), lty=1, cex=0.8, lwd = 2)
dev.off()

pG4 <- fread("../data/G4/G4hunter_w25_s1.5.txt", sep = "\t", header = FALSE) %>% data.frame
RSE.k.tsspg4 <- bt.intersect(a = RSE.k.tssg4, b = pG4, wa = TRUE, c = TRUE) %>% unique
RSE.k.tsspg4$label <- ifelse(RSE.k.tsspg4$V14 > 0, "g1", ifelse(RSE.k.tsspg4$V16 > 0, "g2", "g3"))
RSE.h.tsspg4 <- bt.intersect(a = RSE.h.tssg4, b = pG4, wa = TRUE, c = TRUE) %>% unique
RSE.h.tsspg4$label <- ifelse(RSE.h.tsspg4$V14 > 0, "g1", ifelse(RSE.h.tsspg4$V16 > 0, "g2", "g3"))

RSE.k.tsspg4.gr <- TSSgr(RSE.k.tsspg4)
RSE.h.tsspg4.gr <- TSSgr(RSE.h.tsspg4)

RSE.k.tsspg4Y.gr <- subset(RSE.k.tsspg4.gr, label == "g1")
RSE.k.tsspg4N.gr <- subset(RSE.k.tsspg4.gr, label == "g2")
RSE.k.tsspg4NN.gr <- subset(RSE.k.tsspg4.gr, label == "g3")
RSE.h.tsspg4Y.gr <- subset(RSE.h.tsspg4.gr, label == "g1")
RSE.h.tsspg4N.gr <- subset(RSE.h.tsspg4.gr, label == "g2")
RSE.h.tsspg4NN.gr <- subset(RSE.h.tsspg4.gr, label == "g3")

matp.k.Y <- normalizeToMatrix(tfbs.k.gr, RSE.k.tsspg4Y.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
matp.k.N <- normalizeToMatrix(tfbs.k.gr, RSE.k.tsspg4N.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
matp.k.NN <- normalizeToMatrix(tfbs.k.gr, RSE.k.tsspg4NN.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

matp.h.Y <- normalizeToMatrix(tfbs.h.gr, RSE.h.tsspg4Y.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
matp.h.N <- normalizeToMatrix(tfbs.h.gr, RSE.h.tsspg4N.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)
matp.h.NN <- normalizeToMatrix(tfbs.h.gr, RSE.h.tsspg4NN.gr, 
  value_column = "score", extend = 1000, mean_mode = "coverage", w = 10)

pdf("../figure/cage/tf/K562_pG4_TF_TSS.pdf", width = 5, height = 5)
plot(colMeans(data.frame(matp.k.Y)), type="l", lwd=3, col = "orange", 
  xlab = "", ylab = "Binding density of TFs", xaxt="n", ylim = c(10, 110))
lines(colMeans(data.frame(matp.k.N)), lwd=3, col="steelblue")
lines(colMeans(data.frame(matp.k.NN)), lwd=3, col="green")
axis(side=1, at=c(1,101,201),labels=c("-1kb","Promoter center","1kb"))
legend(1, 100, legend=c("G4-associated promoters", "predicted-only G4-associated promoters", "non-G4-associated promoters"),
       col=c("orange", "steelblue", "green"), lty=1, cex=0.8, lwd = 2)
dev.off()

pdf("../figure/cage/tf/HepG2_pG4_TF_TSS.pdf", width = 5, height = 5)
plot(colMeans(data.frame(matp.h.Y)), type="l", lwd=3, col = "orange", 
  xlab = "", ylab = "Binding density of TFs", xaxt="n", ylim = c(10, 130))
lines(colMeans(data.frame(matp.h.N)), lwd=3, col="steelblue")
lines(colMeans(data.frame(matp.h.NN)), lwd=3, col="green")
axis(side=1, at=c(1,101,201),labels=c("-1kb","Promoter center","1kb"))
legend(1, 120, legend=c("G4-associated promoters", "predicted-only G4-associated promoters", "non-G4-associated promoters"),
       col=c("orange", "steelblue", "green"), lty=1, cex=0.8, lwd = 2)
dev.off()

pdf("../figure/cage/tf/K562_aG4_pG4_exp.pdf", height = 5, width = 4)
boxplot(log2(RSE.k.tsspg4[RSE.k.tsspg4$label == "g1", 6]), log2(RSE.k.tsspg4[RSE.k.tsspg4$label == "g2", 6]), col = c("#ff857c", "#99b898"))
dev.off()
# 1.143587e-174

pdf("../figure/cage/tf/HepG2_aG4_pG4_exp.pdf", height = 5, width = 4)
boxplot(log2(RSE.h.tsspg4[RSE.h.tsspg4$label == "g1", 6]), log2(RSE.h.tsspg4[RSE.h.tsspg4$label == "g2", 6]), col = c("#ff857c", "#99b898"))
dev.off()
# 1.331121e-164
