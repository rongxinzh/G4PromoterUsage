
set.seed(1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(PRROC))

load(file = "SVM_train.Rdata")

fig.k <- ggplot(varImp(svm.k)$importance %>% data.frame %>% select(c1) %>% mutate(feature = rownames(varImp(svm.k)$importance)), 
           aes(x = reorder(feature, c1, decreasing = TRUE), y = c1, label = round(c1, 1))) +
         geom_point(stat='identity', size=6.6, colour = "#9567bd") +
         geom_text(color="white", size=2) +
         #coord_flip() + 
         theme_bw() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))
         #rotate_x_text(90)
ggsave("../figure/cage/imp/SVM_K562_Imp.pdf", fig.k, height = 3.6, width = 8)

fig.h <- ggplot(varImp(svm.h)$importance %>% data.frame %>% select(c1) %>% mutate(feature = rownames(varImp(svm.h)$importance)), 
           aes(x = reorder(feature, c1, decreasing = TRUE), y = c1, label = round(c1, 1))) +
         geom_point(stat='identity', size=6.6, colour = "#6ba3d6") +
         geom_text(color="white", size=2) +
         #coord_flip() + 
         theme_bw()  +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14))
         #rotate_x_text(90)
ggsave("../figure/cage/imp/SVM_HepG2_Imp.pdf", fig.h, height = 3.6, width = 8)

svm.k.train.prob <- predict(svm.k, newdata = select(data.k$train, -score), type = "prob")
svm.k.valid.prob <- predict(svm.k, newdata = select(data.k$valid, -score), type = "prob")
svm.k.test.prob <- predict(svm.k, newdata = select(data.k$test, -score), type = "prob")

svm.k.tr.roc <- roc.curve(scores.class0 = svm.k.train.prob[, 2], 
                          weights.class0 = ifelse((data.k$train$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.k.tr.prc <- pr.curve(scores.class0 = svm.k.train.prob[, 2], 
                         weights.class0 = ifelse((data.k$train$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.k.vd.roc <- roc.curve(scores.class0 = svm.k.valid.prob[, 2], 
                          weights.class0 = ifelse((data.k$valid$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.k.vd.prc <- pr.curve(scores.class0 = svm.k.valid.prob[, 2], 
                         weights.class0 = ifelse((data.k$valid$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.k.ts.roc <- roc.curve(scores.class0 = svm.k.test.prob[, 2], 
                          weights.class0 = ifelse((data.k$test$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.k.ts.prc <- pr.curve(scores.class0 = svm.k.test.prob[, 2], 
                         weights.class0 = ifelse((data.k$test$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)

svm.h.train.prob <- predict(svm.h, newdata = select(data.h$train, -score), type = "prob")
svm.h.valid.prob <- predict(svm.h, newdata = select(data.h$valid, -score), type = "prob")
svm.h.test.prob <- predict(svm.h, newdata = select(data.h$test, -score), type = "prob")

svm.h.tr.roc <- roc.curve(scores.class0 = svm.h.train.prob[, 2], 
                          weights.class0 = ifelse((data.h$train$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.h.tr.prc <- pr.curve(scores.class0 = svm.h.train.prob[, 2], 
                         weights.class0 = ifelse((data.h$train$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.h.vd.roc <- roc.curve(scores.class0 = svm.h.valid.prob[, 2], 
                          weights.class0 = ifelse((data.h$valid$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.h.vd.prc <- pr.curve(scores.class0 = svm.h.valid.prob[, 2], 
                         weights.class0 = ifelse((data.h$valid$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.h.ts.roc <- roc.curve(scores.class0 = svm.h.test.prob[, 2], 
                          weights.class0 = ifelse((data.h$test$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.h.ts.prc <- pr.curve(scores.class0 = svm.h.test.prob[, 2], 
                         weights.class0 = ifelse((data.h$test$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)

svm.k.invd.prob <- predict(svm.k, newdata = select(mat.h, -score), type = "prob")
svm.h.invd.prob <- predict(svm.h, newdata = select(mat.k, -score), type = "prob")

svm.k.invd.roc <- roc.curve(scores.class0 = svm.k.invd.prob[, 2], 
                            weights.class0 = ifelse((mat.h$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.k.invd.prc <- pr.curve(scores.class0 = svm.k.invd.prob[, 2], 
                           weights.class0 = ifelse((mat.h$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)

svm.h.invd.roc <- roc.curve(scores.class0 = svm.h.invd.prob[, 2], 
                            weights.class0 = ifelse((mat.k$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)
svm.h.invd.prc <- pr.curve(scores.class0 = svm.h.invd.prob[, 2], 
                           weights.class0 = ifelse((mat.k$score %>% unlist %>% as.vector) == "c0", 0, 1), curve=TRUE)

pdf("../figure/cage/imp/SVM_K562_train_ROC.pdf", width = 5.5, height = 5)
plot(svm.k.tr.roc)
dev.off()
pdf("../figure/cage/imp/SVM_K562_valid_ROC.pdf", width = 5.5, height = 5)
plot(svm.k.vd.roc)
dev.off()
pdf("../figure/cage/imp/SVM_K562_test_ROC.pdf", width = 5.5, height = 5)
plot(svm.k.ts.roc)
dev.off()
pdf("../figure/cage/imp/SVM_K562_train_PRC.pdf", width = 5.5, height = 5)
plot(svm.k.tr.prc)
dev.off()
pdf("../figure/cage/imp/SVM_K562_valid_PRC.pdf", width = 5.5, height = 5)
plot(svm.k.vd.prc)
dev.off()
pdf("../figure/cage/imp/SVM_K562_test_PRC.pdf", width = 5.5, height = 5)
plot(svm.k.ts.prc)
dev.off()

pdf("../figure/cage/imp/SVM_HepG2_train_ROC.pdf", width = 5.5, height = 5)
plot(svm.h.tr.roc)
dev.off()
pdf("../figure/cage/imp/SVM_HepG2_valid_ROC.pdf", width = 5.5, height = 5)
plot(svm.h.vd.roc)
dev.off()
pdf("../figure/cage/imp/SVM_HepG2_test_ROC.pdf", width = 5.5, height = 5)
plot(svm.h.ts.roc)
dev.off()
pdf("../figure/cage/imp/SVM_HepG2_train_PRC.pdf", width = 5.5, height = 5)
plot(svm.h.tr.prc)
dev.off()
pdf("../figure/cage/imp/SVM_HepG2_valid_PRC.pdf", width = 5.5, height = 5)
plot(svm.h.vd.prc)
dev.off()
pdf("../figure/cage/imp/SVM_HepG2_test_PRC.pdf", width = 5.5, height = 5)
plot(svm.h.ts.prc)
dev.off()

pdf("../figure/cage/imp/SVM_K562_ind_ROC.pdf", width = 5.5, height = 5)
plot(svm.k.invd.roc)
dev.off()
pdf("../figure/cage/imp/SVM_K562_ind_PRC.pdf", width = 5.5, height = 5)
plot(svm.k.invd.prc)
dev.off()

pdf("../figure/cage/imp/SVM_HepG2_ind_ROC.pdf", width = 5.5, height = 5)
plot(svm.h.invd.roc)
dev.off()
pdf("../figure/cage/imp/SVM_HepG2_ind_PRC.pdf", width = 5.5, height = 5)
plot(svm.h.invd.prc)
dev.off()

