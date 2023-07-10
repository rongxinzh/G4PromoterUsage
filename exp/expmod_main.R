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
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(introdataviz))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(GenomicFeatures))

source("./expmod_func.R")
source("./expmod_plot.R")

gencode.path <- "../data/ref/gencode.v19.annotation.gtf"
txdb <- makeTxDbFromGFF(file = gencode.path, format = "gtf", organism = "Homo sapiens", taxonomyId = 9606)

# hg19 <- SeqinfoForUCSCGenome("hg19")
# convertBED2BigWig("../data/cage/GSE118075/k562/GSM3318239_CNhi10920_biologicalRep1-K562-CAGE_ACC.ctss.bed", "../data/cage/GSE118075/k562/K562_rep1_plus.bw", "../data/cage/GSE118075/k562/K562_rep1_minus.bw", hg19)
# convertBED2BigWig("../data/cage/GSE118075/k562/GSM3318240_CNhi10920_biologicalRep2-K562-CAGE_ATG.ctss.bed", "../data/cage/GSE118075/k562/K562_rep2_plus.bw", "../data/cage/GSE118075/k562/K562_rep2_minus.bw", hg19)
# convertBED2BigWig("../data/cage/GSE118075/hepg2/GSM3318231_CNhi10919_biologicalRep1-HepG2-CAGE_ACC.ctss.bed", "../data/cage/GSE118075/hepg2/HepG2_rep1_plus.bw", "../data/cage/GSE118075/hepg2/HepG2_rep1_minus.bw", hg19)
# convertBED2BigWig("../data/cage/GSE118075/hepg2/GSM3318232_CNhi10919_biologicalRep2-HepG2-CAGE_ATG.ctss.bed", "../data/cage/GSE118075/hepg2/HepG2_rep2_plus.bw", "../data/cage/GSE118075/hepg2/HepG2_rep2_minus.bw", hg19)

bsg <- BSgenome.Hsapiens.UCSC.hg19
odb <- org.Hs.eg.db

k562.G4 <- GetG4("K562")
hepg2.G4 <- GetG4("HepG2")
k562.G4.gr <- GetG4gr(k562.G4)
hepg2.G4.gr <- GetG4gr(hepg2.G4)

cage.data.k <- data.frame(Name = c("kr1", "kr2"), 
                          BigWigPlus = c("K562_r1_plus.bw", "K562_r2_plus.bw"), 
                          BigWigMinus = c("K562_r1_minus.bw", "K562_r2_minus.bw"))
rownames(cage.data.k) <- cage.data.k$Name
k562.bw <- GetBWObj(cage.data.k,
                    bwplist = list.files("../data/cage/GSE118075/k562/", "*plus.bw$", full.names = TRUE), 
                    bwmlist = list.files("../data/cage/GSE118075/k562/", "*minus.bw$", full.names = TRUE))
bw.plus.k <- k562.bw$plus
bw.minus.k <- k562.bw$minus

cage.data.h <- data.frame(Name = c("hr1", "hr2"), 
                          BigWigPlus = c("HepG2_r1_plus.bw", "HepG2_r2_plus.bw"), 
                          BigWigMinus = c("HepG2_r1_minus.bw", "HepG2_r2_minus.bw"))
rownames(cage.data.h) <- cage.data.h$Name
hepg2.bw <- GetBWObj(cage.data.h,
                     bwplist = list.files("../data/cage/GSE118075/hepg2/", "*plus.bw$", full.names = TRUE), 
                     bwmlist = list.files("../data/cage/GSE118075/hepg2/", "*minus.bw$", full.names = TRUE))
bw.plus.h <- hepg2.bw$plus
bw.minus.h <- hepg2.bw$minus

CTSSs.k <- quantifyCTSSs(plusStrand = bw.plus.k,
                         minusStrand = bw.minus.k,
                         genome = seqinfo(bsg),
                         design = cage.data.k) %>% calcTPM() %>% calcPooled()

CTSSs.h <- quantifyCTSSs(plusStrand = bw.plus.h,
                         minusStrand = bw.minus.h,
                         genome = seqinfo(bsg),
                         design = cage.data.h) %>% calcTPM() %>% calcPooled()

RSE.k <- quantifyTSSEh(CTSSs.k, txdb)
RSE.h <- quantifyTSSEh(CTSSs.h, txdb)

PlotTSSWidth("../figure/cage/exp/K562_TSS_width.pdf", 
              data.frame(Length = rowRanges(RSE.k) %>% subset(clusterType == "TSS") %>% width()),
              "#9567bd")

PlotTSSWidth("../figure/cage/exp/HepG2_TSS_width.pdf", 
              data.frame(Length = rowRanges(RSE.h) %>% subset(clusterType == "TSS") %>% width()),
              "#6ba3d6")

RSE.k.tss <- extendPrmtr(RSE.k)
RSE.h.tss <- extendPrmtr(RSE.h)

RSE.k.tssg4 <- bt.intersect(a = RSE.k.tss, b = k562.G4, wa = TRUE, c = TRUE) %>% unique
RSE.k.tssg4$label <- ifelse(RSE.k.tssg4[, ncol(RSE.k.tssg4)] == 0, "non G4-associated", "G4-associated")
RSE.h.tssg4 <- bt.intersect(a = RSE.h.tss, b = hepg2.G4, wa = TRUE, c = TRUE) %>% unique
RSE.h.tssg4$label <- ifelse(RSE.h.tssg4[, ncol(RSE.h.tssg4)] == 0, "non G4-associated", "G4-associated")

PlotG4Prop("../figure/cage/exp/G4_prop_TSS.pdf", RSE.k.tssg4, RSE.h.tssg4)

RSE.k.tssbar <- table(RSE.k.tssg4[, c("V12", "label")]) %>% data.frame() %>% spread(label, Freq) %>% as.data.frame
rownames(RSE.k.tssbar) <- RSE.k.tssbar[, 1]
RSE.k.tssbar <- RSE.k.tssbar[, -1] %>% as.matrix(ncol=2)

RSE.h.tssbar <- table(RSE.h.tssg4[, c("V12", "label")]) %>% data.frame() %>% spread(label, Freq) %>% as.data.frame
rownames(RSE.h.tssbar) <- RSE.h.tssbar[, 1]
RSE.h.tssbar <- RSE.h.tssbar[, -1] %>% as.matrix(ncol=2)

PlotSpineG4Prop("../figure/cage/exp/K562_TSS_G4_prop.pdf", RSE.k.tssbar, c("#c6afd5", "#9567bd"))
PlotSpineG4Prop("../figure/cage/exp/HepG2_TSS_G4_prop.pdf", RSE.h.tssbar, c("#b5dffc", "#6ba3d6"))

RSE.tss.exp <- bind_rows(
                 data.frame(exp = log2(RSE.k.tssg4[RSE.k.tssg4$label == "G4-associated", "V6"]), group = "G4-associated", cell = "K562"),
                 data.frame(exp = log2(RSE.k.tssg4[RSE.k.tssg4$label == "non G4-associated", "V6"]), group = "non G4-associated", cell = "K562"),
                 data.frame(exp = log2(RSE.h.tssg4[RSE.h.tssg4$label == "G4-associated", "V6"]), group = "G4-associated", cell = "HepG2"),
                 data.frame(exp = log2(RSE.h.tssg4[RSE.h.tssg4$label == "non G4-associated", "V6"]), group = "non G4-associated", cell = "HepG2"))

PlotGroupedG4TssExp("../figure/cage/exp/G4group_TSSexp.pdf", RSE.tss.exp)

RSE.k.tssp <- table(RSE.k.tssg4[, "V12"]) %>% data.frame
RSE.k.tssp$perc <- paste0(round(RSE.k.tssp$Freq / sum(RSE.k.tssp$Freq), 3) * 100, "%")

RSE.h.tssp <- table(RSE.h.tssg4[, "V12"]) %>% data.frame
RSE.h.tssp$perc <- paste0(round(RSE.h.tssp$Freq / sum(RSE.h.tssp$Freq), 3) * 100, "%")

PlotTSSFreq("../figure/cage/exp/TSS_region_freq.pdf", RSE.k.tssp, RSE.h.tssp, "#e9c49b")

high.tss.k <- subset(RSE.k, (clusterType == 'TSS') & (score / ncol(RSE.k) >= 10))
high.tss.k <- calcShape(high.tss.k, 
                        pooled = CTSSs.k, 
                        shapeFunction = shapeIQR, 
                        lower = 0.10, 
                        upper = 0.90)

high.tss.h <- subset(RSE.h, (clusterType == 'TSS') & (score / ncol(RSE.h) >= 10))
high.tss.h <- calcShape(high.tss.h, 
                        pooled = CTSSs.h, 
                        shapeFunction = shapeIQR, 
                        lower = 0.10, 
                        upper = 0.90)

PlotLenIQR("../figure/cage/exp/K562_TSS_lenIQR.pdf", high.tss.k, "#9567bd")
PlotLenIQR("../figure/cage/exp/HepG2_TSS_lenIQR.pdf", high.tss.h, "#6ba3d6")

rowData(high.tss.k)$shape <- ifelse(rowData(high.tss.k)$IQR < 10, "Sharp", "Broad")
rowData(high.tss.h)$shape <- ifelse(rowData(high.tss.h)$IQR < 10, "Sharp", "Broad")

PlotShapeProp("../figure/cage/exp/highTSS_shape_prop.pdf", 
              rowData(high.tss.k) %>% data.frame,
              rowData(high.tss.h) %>% data.frame)

prmtr.seqs.k <- high.tss.k %>% rowRanges() %>% swapRanges() %>% promoters(upstream = 100, downstream = 10) %>% getSeq(bsg, .)
prmtr.seqs.h <- high.tss.h %>% rowRanges() %>% swapRanges() %>% promoters(upstream = 100, downstream = 10) %>% getSeq(bsg, .)

PlotTSSMotif("../figure/cage/exp/highTSS_motif_K562.pdf", prmtr.seqs.k, high.tss.k)
PlotTSSMotif("../figure/cage/exp/highTSS_motif_HepG2.pdf", prmtr.seqs.h, high.tss.h)

PlotTSSExp("../figure/cage/exp/example.png", txdb, RSE.k, CTSSs.k, k562.G4.gr)

save(list = c("RSE.k.tss", "RSE.h.tss"), file = "RSE.tss.Rdata")

k562.G4.rand <- bt.shuffle(k562.G4[, 1:3], "../data/ref/hg19.chrom.sizes.txt", noOverlapping = TRUE, chrom = TRUE)
k562.G4.rand.gr <- GetG4gr(k562.G4.rand)
k562.G4.rand.gr <- GRanges(k562.G4.rand.gr, state = "rand")
hepg2.G4.rand <- bt.shuffle(hepg2.G4[, 1:3], "../data/ref/hg19.chrom.sizes.txt", noOverlapping = TRUE, chrom = TRUE)
hepg2.G4.rand.gr <- GetG4gr(hepg2.G4.rand)
hepg2.G4.rand.gr <- GRanges(hepg2.G4.rand.gr, state = "rand")
k562.G4.gr <- GRanges(k562.G4.gr, state = "G4")
hepg2.G4.gr <- GRanges(hepg2.G4.gr, state = "G4")

k562.G4.gr.all <- c(k562.G4.gr, k562.G4.rand.gr)
hepg2.G4.gr.all <- c(hepg2.G4.gr, hepg2.G4.rand.gr)

mat.k <- normalizeToMatrix(k562.G4.gr.all, rowRanges(RSE.k) %>% subset(clusterType == "TSS"),
                           value_column = "state", extend = 1000, mean_mode = "coverage", w = 10)
pdf("../figure/cage/exp/K562_G4_TSS_distrbution.pdf", height = 15, width = 2)
EnrichedHeatmap(mat.k, name = "Group")
dev.off()

mat.h <- normalizeToMatrix(hepg2.G4.gr.all, rowRanges(RSE.h) %>% subset(clusterType == "TSS"),
                           value_column = "state", extend = 1000, mean_mode = "coverage", w = 10)
pdf("../figure/cage/exp/HepG2_G4_TSS_distrbution.pdf", height = 15, width = 2)
EnrichedHeatmap(mat.h, name = "Group")
dev.off()
