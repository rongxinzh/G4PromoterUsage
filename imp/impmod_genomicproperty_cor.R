
set.seed(1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(tidyr))

load(file = "mat_for_ML.Rdata")

colnames(mat.k.unscaled)[1] <- "DNA Methylation"
colnames(mat.h.unscaled)[1] <- "DNA Methylation"

k.corr <- round(cor(select(mat.k.unscaled, -score), method = "spearman"), 2)
h.corr <- round(cor(select(mat.h.unscaled, -score), method = "spearman"), 2)

k.p.mat <- cor_pmat(select(mat.k.unscaled, -score), method = "spearman")
h.p.mat <- cor_pmat(select(mat.h.unscaled, -score), method = "spearman")

p.k <- ggcorrplot(k.corr,
                  hc.order = TRUE,
                  type = "upper",
                  outline.color = "white",
                  colors = c("#6D9EC1", "white", "#E46726"),
                  p.mat = k.p.mat, tl.srt = 90, tl.cex = 16) + theme(axis.text.x = element_text(vjust = 0.5))


p.h <- ggcorrplot(h.corr,
                  hc.order = TRUE,
                  type = "upper",
                  outline.color = "white",
                  colors = c("#6D9EC1", "white", "#E46726"),
                  p.mat = h.p.mat, tl.srt = 90, tl.cex = 16) + theme(axis.text.x = element_text(vjust = 0.5))

ggsave("../figure/cage/imp/cor_K562.pdf", p.k, width = 10, height = 10)
ggsave("../figure/cage/imp/cor_HepG2.pdf", p.h, width = 10, height = 10)

cor.tpm.k <- NULL
for (i in 1:(ncol(mat.k.unscaled)-1)) {
  d1 <- mat.k.unscaled[, i] %>% unlist %>% as.vector %>% as.numeric
  d2 <- mat.k.unscaled$score %>% unlist %>% as.vector
  ct <- cor.test(d1, d2, method = "spearman")
  cor.tpm.k <- bind_rows(cor.tpm.k, data.frame(rho = ct$estimate, p.val = ct$p.value, property = colnames(mat.k.unscaled)[i]))
}
cor.tpm.k$col <- ifelse(cor.tpm.k$p.val > 0.05, "GC", ifelse(abs(cor.tpm.k$rho) > 0.2, "GA", "GB"))
cor.tpm.k <- cor.tpm.k %>% arrange(desc(rho))
cor.tpm.k$property <- factor(cor.tpm.k$property, levels = cor.tpm.k$property)

fig.k <- ggplot(cor.tpm.k, aes(x = property, y = rho, fill = factor(col))) +
         geom_col() +
         theme_bw() +
         scale_fill_manual(values=c("#ea4c4d", "#3a539a", "#949494")) +
         theme(legend.position = "off", 
               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
               axis.title.x = element_blank())

ggsave("../figure/cage/imp/K562_cor_tpm.pdf", fig.k, width = 8, height = 3.6)

cor.tpm.h <- NULL
for (i in 1:(ncol(mat.h.unscaled)-1)) {
  d1 <- mat.h.unscaled[, i] %>% unlist %>% as.vector %>% as.numeric
  d2 <- mat.h.unscaled$score %>% unlist %>% as.vector
  ct <- cor.test(d1, d2, method = "spearman")
  cor.tpm.h <- bind_rows(cor.tpm.h, data.frame(rho = ct$estimate, p.val = ct$p.value, property = colnames(mat.h.unscaled)[i]))
}
cor.tpm.h$col <- ifelse(cor.tpm.h$p.val > 0.05, "GC", ifelse(abs(cor.tpm.h$rho) > 0.2, "GA", "GB"))
cor.tpm.h <- cor.tpm.h %>% arrange(desc(rho))
cor.tpm.h$property <- factor(cor.tpm.h$property, levels = cor.tpm.h$property)

fig.h <- ggplot(cor.tpm.h, aes(x = property, y = rho, fill = factor(col))) +
         geom_col() +
         theme_bw() +
         scale_fill_manual(values=c("#ea4c4d", "#3a539a", "#949494")) +
         theme(legend.position = "off", 
               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
               axis.title.x = element_blank())

ggsave("../figure/cage/imp/HepG2_cor_tpm.pdf", fig.h, width = 8, height = 3.6)

mat.k.unscaled$score <- ifelse(mat.k.unscaled$score > quantile(mat.k.unscaled$score)[4], "c1", ifelse(mat.k.unscaled$score < quantile(mat.k.unscaled$score)[2], "c0", "c2"))
mat.k.unscaled <- mat.k.unscaled[mat.k.unscaled$score != "c2",]

mat.h.unscaled$score <- ifelse(mat.h.unscaled$score > quantile(mat.h.unscaled$score)[4], "c1", ifelse(mat.h.unscaled$score < quantile(mat.h.unscaled$score)[2], "c0", "c2"))
mat.h.unscaled <- mat.h.unscaled[mat.h.unscaled$score != "c2",]

mat.k.unscaled <- gather(mat.k.unscaled, property, value, `DNA Methylation`:H4K91ac)
mat.h.unscaled <- gather(mat.h.unscaled, property, value, `DNA Methylation`:H4K91ac)

mat.k.unscaled.p <- mat.k.unscaled %>% 
                      group_by(property, score) %>% 
                      summarise(value = list(value)) %>% 
                      spread(score, value) %>% 
                      group_by(property) %>% 
                      mutate(pval = wilcox.test(unlist(c1), unlist(c0))$p.value) %>% select(property, pval) %>% 
                      arrange(pval) %>% data.frame() 

mat.h.unscaled.p <- mat.h.unscaled %>% 
                      group_by(property, score) %>% 
                      summarise(value = list(value)) %>% 
                      spread(score, value) %>% 
                      group_by(property) %>% 
                      mutate(pval = wilcox.test(unlist(c1), unlist(c0))$p.value) %>% select(property, pval) %>% 
                      arrange(pval) %>% data.frame() 

mat.k.unscaled$property <- factor(mat.k.unscaled$property, levels = mat.k.unscaled.p$property)
mat.h.unscaled$property <- factor(mat.h.unscaled$property, levels = mat.h.unscaled.p$property)

f1 <- ggplot(data = mat.k.unscaled, mapping = aes(x = property, y = log2(value+1), fill = score)) +
      geom_boxplot(position = position_dodge(0.5)) + 
      scale_fill_manual(values=c("#0092d1","#62b232")) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.text = element_text(size=10), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14)) + 
      ylab("Signal values (log2)")
ggsave("../figure/cage/imp/K562_c0c1_values.pdf", f1, width = 8, height = 3.6)

f2 <- ggplot(data = mat.h.unscaled, mapping = aes(x = property, y = log2(value+1), fill = score)) +
      geom_boxplot(position = position_dodge(0.5)) + 
      scale_fill_manual(values=c("#0092d1","#62b232")) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.text = element_text(size=10), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14)) + 
      ylab("Signal values (log2)")
ggsave("../figure/cage/imp/HepG2_c0c1_values.pdf", f2, width = 8, height = 3.6)
