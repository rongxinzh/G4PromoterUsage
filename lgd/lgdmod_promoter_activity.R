
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(proActiv))
suppressPackageStartupMessages(library(EnrichedHeatmap))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(edgeR))

sample.info.path <- "../data/RNA-seq/sample_info.csv"
gtf.path <- "../data/ref/gencode.v19.annotation.gtf"

prmtr.anno <- preparePromoterAnnotation(file = gtf.path, species = 'Homo_sapiens')

sample.info <- fread(sample.info.path, sep = ",", header = TRUE) %>% data.frame
all.ID <- unique(sample.info[, "ID"])

GrG4 <- function() {
  G4 <- data.table::fread("../data/G4/G4hunter_w25_s1.5.txt", sep = "\t", header = FALSE) %>% data.frame
  G4 <- G4[, c(1:3, 5)]
  colnames(G4) <- c("chr", "start", "end", "strand")
  G4$score <- 1
  G4 <- GenomicRanges::makeGRangesFromDataFrame(G4,
                                                keep.extra.columns = TRUE, 
                                                ignore.strand = FALSE)

  return(G4)
}

GetDEP <- function(data = NULL, design = NULL, cont = NULL) {

    rownames(design)=colnames(data)
    cont.mat <- makeContrasts(contrasts = paste0(unique(cont), collapse = "-"), levels = design)
    #cont.mat <- makeContrasts(paste0(unique(cont), collapse = "-"), levels = design)
    data <- log2(data + 1)
    fit <- lmFit(data, design)
    fit2 <- contrasts.fit(fit, cont.mat)
    fit2 <- eBayes(fit2) 
    output = topTable(fit2, coef = 1, n = Inf)
    deg.data = na.omit(output)
    deg.data$promoterId <- as.numeric(rownames(deg.data))

    return(deg.data)
}

pG4 <- GrG4()
up.prmtr <- list()
down.prmtr <- list()
up.pgene <- list()
down.pgene <- list() 

pv <- NULL

for (tmp.round in all.ID) {
  tmp.sample <- sample.info %>% filter(ID == tmp.round)
  tmp.ID <- unique(tmp.sample$AccessionID1)
  tmp.sample$path <- paste0("../data/RNA-seq/", tmp.ID, "/align/", tmp.sample$AccessionID2, "_SJ.out.tab")
  tmp.treat <- setdiff(unique(tmp.sample$Ligand), "Control")

  for (tmp.tt in tmp.treat) {

    tmp.ct.path <- tmp.sample %>% filter(Ligand == "Control") %>% select(path) %>% unlist %>% as.vector
    tmp.tt.path <- tmp.sample %>% filter(Ligand == tmp.tt) %>% select(path) %>% unlist %>% as.vector
    tmp.files <- c(tmp.ct.path, tmp.tt.path)

    tmp.condition <- c(rep("Control", length(tmp.ct.path)), rep(tmp.tt, length(tmp.tt.path))) 

    tmp.activ <- proActiv(files = tmp.files, 
                          promoterAnnotation = prmtr.anno,
                          condition = tmp.condition)
    tmp.activ <- tmp.activ[complete.cases(assays(tmp.activ)$promoterCounts),]
    tmp.activ <- tmp.activ[rowData(tmp.activ)$internalPromoter == FALSE, ]

    tmp.dge <- DGEList(counts=assays(tmp.activ)$promoterCounts, group=colData(tmp.activ)$condition, 
                       genes=rowData(tmp.activ)[, c("geneId", "promoterId")])

    tmp.keep <- rowSums(cpm(tmp.dge) > 1 ) >= 1
    tmp.dge <- tmp.dge[tmp.keep, , keep.lib.sizes = FALSE]
    tmp.dge <- calcNormFactors(tmp.dge)
    
    tmp.design = model.matrix(~colData(tmp.activ)$condition)

    tmp.disp <- estimateDisp(tmp.dge, design = tmp.design, tagwise = FALSE)
    tmp.QLfit <- glmQLFit(tmp.disp, design = tmp.design, robust = TRUE)

    tmp.qlf <- glmQLFTest(tmp.QLfit)
    tmp.qlf.total <- topTags(tmp.qlf, n = 1000000) %>% data.frame

    tmp.ds <- diffSpliceDGE(tmp.QLfit, geneid = "geneId")
    tmp.ds.tss <- topSpliceDGE(tmp.ds, test = "exon", n = 1000000)

    x1 <- tmp.ds.tss %>% filter(logFC > log2(2), P.Value < 0.05)
    x2 <- tmp.ds.tss %>% filter(logFC < -log2(2), P.Value < 0.05)
    x3 <- length(intersect(unique(x1$geneId), unique(x2$geneId)))
    x4 <- length(unique(x1$geneId)) - x3
    x5 <- length(unique(x2$geneId)) - x3

    pv <- bind_rows(pv, 
      data.frame(GSE.ID = tmp.ID, ligand = tmp.tt, up = x4, down = x5, 
        intersection = x3,
        pval = 0.05, FDR = 0
        ))

    pdf(paste0("../figure/rnaseq/lgd/", tmp.round, "_", tmp.ID, "_", tmp.tt, "_Prmtr_Volcano.pdf"))
    print(EnhancedVolcano(tmp.qlf.total, lab = str_split(tmp.qlf.total$geneId, "\\.", simplify = TRUE)[, 1], x = 'logFC', y = 'PValue', pCutoff = 0.05, FCcutoff = log2(2)))
    dev.off()

  }

}

fwrite(pv, "../data/RNA-seq/pv.csv", col.names = TRUE, row.names = FALSE)

