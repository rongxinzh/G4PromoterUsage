
set.seed(1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(doParallel))

load(file = "mat_for_ML.Rdata")

DataSplit <- function(data = NULL) {

  idx <- sample(seq(1, 3), size = nrow(data), replace = TRUE, prob = c(.6, .2, .2))
  train <- data[idx == 1,]
  valid <- data[idx == 2,]
  test <- data[idx == 3,]

  return(list(train = train, valid = valid, test = test))
}

CVsearch <- function(vdata = NULL, C.seq = NULL, sigma.seq = NULL) {

  cv.control <- trainControl(method = 'repeatedcv',
                             number = 10,
                             repeats = 10,
                             search = 'grid',
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary)

  cv.tune.grid <- expand.grid(C = C.seq, sigma = sigma.seq) 
  cv.model.list <- list()

  cl <- makePSOCKcluster(20)
  registerDoParallel(cl)

  fit <- train(score~.,
               data = vdata,
               method = 'svmRadial',
               metric = 'ROC',
               tuneGrid = cv.tune.grid,
               trControl = cv.control)
  roc.res <- data.frame(fit$results)

  stopCluster(cl)

  return(roc.res)
}

data.k <- DataSplit(mat.k)
data.h <- DataSplit(mat.h)

c.seq <- c(0.1, 0.5, seq(1, 50, 1))
sigma.seq <- c(0.01, 0.1, 0.5, seq(1, 10, 1))

roc.res.k <- CVsearch(data.k$valid, c.seq, sigma.seq)
roc.res.h <- CVsearch(data.h$valid, c.seq, sigma.seq)

roc.max.k <- roc.res.k[roc.res.k$ROC == max(roc.res.k$ROC), , drop = FALSE]
roc.max.h <- roc.res.h[roc.res.h$ROC == max(roc.res.h$ROC), , drop = FALSE]

t.control <- trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary)

svm.k <- train(score ~ ., 
               data = data.k$train, 
               method = "svmRadial", 
               trControl = t.control, 
               verbose = FALSE, 
               tuneGrid = data.frame(C = roc.max.k$C, sigma = roc.max.k$sigma),
               metric = "ROC")

svm.h <- train(score ~ ., 
               data = data.h$train, 
               method = "svmRadial", 
               trControl = t.control, 
               verbose = FALSE, 
               tuneGrid = data.frame(C = roc.max.h$C, sigma = roc.max.h$sigma),
               metric = "ROC")

save(list=ls(), file = "./SVM_train.Rdata")
