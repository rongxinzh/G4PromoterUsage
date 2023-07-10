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
suppressPackageStartupMessages(library(bedtoolsr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(introdataviz))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggsankey))
suppressPackageStartupMessages(library(ggpie))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(GenomicFeatures))

bsg <- BSgenome.Hsapiens.UCSC.hg19
gencode.path <- "../data/ref/gencode.v19.annotation.gtf"
txdb <- makeTxDbFromGFF(file = gencode.path, format = "gtf", organism = "Homo sapiens", taxonomyId = 9606)
odb <- org.Hs.eg.db

k562.G4 <- fread("../data/G4/cellline/K562_hg19.bed", sep = "\t", header = FALSE) %>% data.frame
hepg2.G4 <- fread("../data/G4/cellline/HepG2_hg19.bed", sep = "\t", header = FALSE) %>% data.frame
k562.G4.gr <- k562.G4
k562.G4.gr$score <- 1
colnames(k562.G4.gr)[1:3] <- c("chr", "start", "end")
k562.G4.gr <- makeGRangesFromDataFrame(k562.G4.gr,
                                       keep.extra.columns = TRUE,
                                       ignore.strand = TRUE)
hepg2.G4.gr <- hepg2.G4
hepg2.G4.gr$score <- 1
colnames(hepg2.G4.gr)[1:3] <- c("chr", "start", "end")
hepg2.G4.gr <- makeGRangesFromDataFrame(hepg2.G4.gr,
                                        keep.extra.columns = TRUE,
                                        ignore.strand = TRUE)

# create cage information object
cage.data <- bind_rows(
data.frame(Class = "K562", Name = c("kr1", "kr2"), 
           BigWigPlus = c("K562_r1_plus.bw", "K562_r2_plus.bw"), 
           BigWigMinus = c("K562_r1_minus.bw", "K562_r2_minus.bw")),
data.frame(Class = "HepG2", Name = c("hr1", "hr2"), 
           BigWigPlus = c("HepG2_r1_plus.bw", "HepG2_r2_plus.bw"), 
           BigWigMinus = c("HepG2_r1_minus.bw", "HepG2_r2_minus.bw"))
)
cage.data$Class = factor(cage.data$Class, level = c("K562", "HepG2"))
rownames(cage.data) <- cage.data$Name

# create bw object
bw.plus <- c(list.files("../data/cage/GSE118075/k562/", "*plus.bw$", full.names = TRUE), list.files("../data/cage/GSE118075/hepg2/", "*plus.bw$", full.names = TRUE))
bw.minus <- c(list.files("../data/cage/GSE118075/k562/", "*minus.bw$", full.names = TRUE), list.files("../data/cage/GSE118075/hepg2/", "*minus.bw$", full.names = TRUE))
bw.plus <- BigWigFileList(bw.plus)
bw.minus <- BigWigFileList(bw.minus)
names(bw.plus) <- names(bw.minus) <- cage.data$Name

CTSSs <- quantifyCTSSs(plusStrand = bw.plus,
                       minusStrand = bw.minus,
                       genome = seqinfo(bsg),
                       design = cage.data)
CTSSs <- CTSSs %>% calcTPM() %>% calcPooled()
TCs <- quickTSSs(CTSSs)
TSSs <- TCs %>% calcTPM() %>%
        subsetBySupport(inputAssay = "TPM", 
                        unexpressed = 1, 
                        minSamples = 1)

BCs <- quickEnhancers(CTSSs)
BCs <- subsetBySupport(BCs, 
                       inputAssay = "counts", 
                       unexpressed = 0, 
                       minSamples = 1)

TSSs <- assignTxID(TSSs, txModels = txdb, swap = "thick")
TSSs <- assignTxType(TSSs, txModels = txdb, swap = "thick")

BCs <- assignTxType(BCs, txModels = txdb, swap = "thick")
Enhancers <- subset(BCs, txType %in% c("intergenic", "intron"))

TSSs$totalTags <- NULL
Enhancers$totalTags <- NULL

rowData(TSSs)$balance <- NA
rowData(TSSs)$bidirectionality <- NA
rowData(Enhancers)$txID <- NA

rowData(TSSs)$clusterType <- "TSS"
rowData(Enhancers)$clusterType <- "Enhancer"

RSE <- combineClusters(object1 = TSSs, 
                       object2 = Enhancers, 
                       removeIfOverlapping = "object1")
RSE <- calcTPM(RSE) %>% calcPooled()
RSE <- assignGeneID(RSE, geneModels = txdb)

RSE.filtered <- RSE %>%
    subset(clusterType == "TSS" & !is.na(geneID)) %>%
    subsetByComposition(inputAssay = "counts", 
                        genes = "geneID", 
                        unexpressed = 0.1, 
                        minSamples = 1)

rowData(RSE.filtered)$symbol <- mapIds(odb, 
                                       keys = str_split(rowData(RSE.filtered)$geneID, "\\.", simplify = TRUE)[, 1],
                                       column = "SYMBOL", 
                                       keytype = "ENSEMBL")

RSE.filtered <- RSE.filtered %>% subset(!is.na(symbol))

pdf("../figure/cage/ap/Multi_TSS_freq.pdf", width = 5, height = 2.6)
RSE.filtered %>% rowData() %>% as.data.frame() %>% as_tibble() %>% dplyr::count(geneID) %>% 
    ggplot(aes(x = n, fill = n >= 2)) + 
    geom_bar(alpha = 0.75) +
    scale_fill_colorblind("Multi-TSS") +
    labs(x = "Number of TSSs per gene", y = "Frequency") +
    theme_classic()
dev.off()

TSS.info <- RSE.filtered %>% rowData() %>% subset(select = c(score, txType, geneID, symbol)) %>% as.data.frame()

dge <- DGEList(counts = assay(RSE.filtered, "counts"), genes = TSS.info)

dge <- calcNormFactors(dge)
mod <- model.matrix(~ Class, data = colData(RSE.filtered))

disp <- estimateDisp(dge, design = mod, tagwise = FALSE)
QLfit <- glmQLFit(disp, design = mod, robust = TRUE)

ds <- diffSpliceDGE(QLfit, geneid = "geneID")
ds.tss <- topSpliceDGE(ds, test = "exon", n = 1000000)

ds.tss.info <- data.frame(chr = str_split(rownames(ds.tss), ":", simplify = TRUE)[, 1], 
                          start = as.numeric(str_split(rownames(ds.tss), "(:|-)", simplify = TRUE)[, 2]), 
                          end = as.numeric(str_split(rownames(ds.tss), "(;|-)", simplify = TRUE)[, 2]),
                          strand = str_split(rownames(ds.tss), ";", simplify = TRUE)[, 2],
                          tssInfo = rownames(ds.tss), ds.tss)
ds.tss.info$start <- ifelse(ds.tss.info$strand == "+", ds.tss.info$start - 300, ds.tss.info$start)
ds.tss.info$end <- ifelse(ds.tss.info$strand == "+", ds.tss.info$end, ds.tss.info$end + 300)
ds.tss.info <- bt.intersect(a = ds.tss.info, b = k562.G4.gr, wa = TRUE, c = TRUE) %>% unique
ds.tss.info <- bt.intersect(a = ds.tss.info, b = hepg2.G4.gr, wa = TRUE, c = TRUE) %>% unique
ds.tss.info$anno <- ifelse((ds.tss.info$V14 > 0) & (ds.tss.info$V15 > 0), 2, ifelse((ds.tss.info$V14 == 0) & (ds.tss.info$V15 == 0), 0, 1))
colnames(ds.tss.info) <- c("chr", "start", "end", "strand", "tssInfo", colnames(ds.tss), "K562", "HepG2", "anno")

ds.tss.up <- ds.tss.info %>% filter(logFC > log2(2), P.Value < 0.05)
ds.tss.down <- ds.tss.info %>% filter(logFC < -log2(2), P.Value < 0.05)
ds.tss.dys <- bind_rows(data.frame(ds.tss.up, group = "up"), data.frame(ds.tss.down, group = "down"))

pdf("../figure/cage/ap/AP_count.pdf", width = 5, height = 5)
barplot(c(dim(ds.tss.up)[1], dim(ds.tss.down)[1]), ylim = c(0, 800), width = 1)
dev.off()

venn.diagram(list(Up = unique(ds.tss.up$symbol), Down = unique(ds.tss.down$symbol)), fill = c("#ecc991", "#76bccc"), 
             lwd = 1, filename = "../figure/cage/ap/AP_updown_venn.png", width = 1600, height = 1600)

ds.tss.dys.anno <- ds.tss.dys %>% group_by(tssInfo) %>% summarise(x = sum(anno)) %>% data.frame 
ds.tss.dys.anno$x <- ifelse(ds.tss.dys.anno$x == 0, "None", "Contain >= 1 G4s")
pdf("../figure/cage/ap/TSS_contain_G4TSS_prop.pdf", width = 5, height = 5)
ggdonut(data = ds.tss.dys.anno, group_key = "x", count_type = "full",
        label_info = "all", label_type = "horizon", donut.label.size = 6,
        label_size = 6, label_pos = "in", label_threshold = 10)
dev.off()

ap <- intersect(unique(ds.tss.up$symbol), unique(ds.tss.down$symbol)) %>% unique
ds.tss.ap <- ds.tss.dys %>% filter(symbol %in% ap)

ds.tss.info.adds <- table(ds.tss.info %>% filter(anno > 0) %>% dplyr::select(anno)) %>% as.vector()
ds.tss.info.adds <- ds.tss.info.adds[1] / ds.tss.info.adds[2]

ds.tss.dys.adds <- table(ds.tss.dys %>% filter(anno > 0) %>% dplyr::select(anno)) %>% as.vector()
ds.tss.dys.adds <- ds.tss.dys.adds[1] / ds.tss.dys.adds[2]

ds.tss.prop <- data.frame(prop = round(c(ds.tss.info.adds, ds.tss.dys.adds), 3),
                          group = c("Background", "Alternative TSS"))

pdf("../figure/cage/ap/fold.pdf", width = 3.6, height = 3.6)
ggbarplot(ds.tss.prop, "group", "prop", width = 0.6,
  fill = "group", color = "group", palette = "Paired",
  label = TRUE, lab.pos = "in", lab.col = "white") + 
border() + rremove("xlab") + rremove("legend") + ylab("Fold")
dev.off()

ds.apgene.G4 <- ds.tss.ap %>% group_by(symbol) %>% summarise(x = sum(anno)) %>% data.frame 
ds.apgene.G4$x <- ifelse(ds.apgene.G4$x == 0, "None", "Contain >= 1 G4s")

pdf("../figure/cage/ap/Gene_contain_G4TSS_prop.pdf", width = 5, height = 5)
ggdonut(data = ds.apgene.G4, group_key = "x", count_type = "full",
        label_info = "all", label_type = "horizon", donut.label.size = 6,
        label_size = 6, label_pos = "in", label_threshold = 10)
dev.off()

wilcox.test(ds.tss.info %>% filter(anno == 1) %>% dplyr::select(logFC) %>% unlist() %>% as.vector() %>% abs(),
            ds.tss.info %>% filter(anno == 2) %>% dplyr::select(logFC) %>% unlist() %>% as.vector() %>% abs(), 
            alternative = "greater" )

wilcox.test(ds.tss.dys %>% filter(anno == 1) %>% dplyr::select(logFC) %>% unlist() %>% as.vector() %>% abs(),
            ds.tss.dys %>% filter(anno == 2) %>% dplyr::select(logFC) %>% unlist() %>% as.vector() %>% abs(), 
            alternative = "greater" )

ds.tss.dys.ft <- ds.tss.dys %>% filter(anno > 0)
ds.tss.dys.ft$logFC <- abs(ds.tss.dys.ft$logFC)
pdf("../figure/cage/ap/Gene_contain_G4TSS_prop_dys.pdf", width = 4, height = 4)
ggboxplot(ds.tss.dys.ft, "anno", "logFC",
          fill = "anno", palette = c("#347ebd", "#5bb85c")) + border() +
rremove("xlab") + rremove("legend") +
ylab("absolute log2FC")
dev.off()

