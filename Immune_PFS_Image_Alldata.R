# packages ----------------------------------------------------------------

rm(list=ls())
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = F)
seed <- 123
# install.packages("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/R_packages/yulab.utils_0.2.0.tgz", repos = NULL, type = "mac.binary")
# 安装和加载相关包
# install.packages(c("clusterProfiler", "limma", "maftools", "ESTIMATE"))
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")
# 然后重新安装clusterProfiler
# 卸载已有的 yulab.utils 包
# 卸载 clusterProfiler 包

# 然后安装相关包
# BiocManager::install("ESTIMATE")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("CIBERSORT")

# install.packages("DESeq2")  # 如果尚未安装 DESeq2
# install.packages('CIBERSORT')
# library('devtools')
# devtools::install_github("https://github.com/PoisonAlien/TCGAmutations")
# 如果没有安装，先安装 IOBR 包
# install.packages("devtools")
## 官方安装教程
## 安装依赖包
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# depens<-c('tibble', 'survival', 'survminer', 'limma', "DESeq2","devtools", 'limSolve', 'GSVA', 'e1071', 'preprocessCore', 
#           "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor",  "timeROC", "pracma", "factoextra", 
#           "FactoMineR", "WGCNA", "patchwork", 'ggplot2', "biomaRt", 'ggpubr')
# 
# for(i in 1:length(depens)){
#   depen<-depens[i]
#   if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
# }

#安装IOBR包
# devtools::install_github("IOBR/IOBR") # 安装IOBR包
# 加载必要的包
# install.packages("meta")
# install.packages("forestploter")
library(survival)    # 用于Cox回归分析
library(meta)        # 用于Meta分析
library(metafor)
library(forestploter) # 用于绘制森林图
library(tidyverse)
library(maftools)
library(BiocManager)
# BiocManager::install("RTCGA.mutations")
# BiocManager::install(c("TCGAutils", "SummarizedExperiment"))
library(TCGAutils)
library(SummarizedExperiment)
# BiocManager::install("PoisonAlien/TCGAmutations")
devtools::install_local("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/R_packages/TCGAmutations-master")
# devtools::install_github("BioinformaticsFMRP/TCGAmutations")
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(RTCGA.mRNA)
library(RTCGA.mutations)
# TCGAmutations::tcga_available()
library(TCGAmutations)
library(visdat)
library(compareGroups)
library(ggstatsplot)
library(IOBR)
library(doParallel)
library(ggprism) ##设置主题私用
library(ggrepel) ##给火山图加标签使用
library(DESeq2)
library(CIBERSORT)
library(edgeR)
# library(org.Hs.eg.db)
# install.packages("ESTIMATE")
# library(R.utils)
library(readr)
library(enrichplot)
library(tidyverse)
library(survminer)
library(survival)

library(clusterProfiler)
library(limma)
library(maftools)
# library(ESTIMATE)

library(openxlsx)
library(timeROC)
library(rms)
library(quadprog)
library(parallel)
library(dplyr)
library(lava)
library(snowfall)
library(caret)
library(randomForestSRC)
library(randomForest)
library(glmnet)
library(pls)
library(plsRcox)
library(plsRglm)
library(superpc)
library(gbm)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(e1071)
library(xgboost)
library(Matrix)
library(pROC)  
library(StepReg)
library(naivebayes)
library(MASS)
library(readxl)
library(dplyr)
library(tidyr)
library(caret)
library(VIM)
library(ggplot2)
library(pROC)
library(randomForest)
library(e1071)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(devtools)
# install_github("binderh/CoxBoost")
library(CoxBoost)

# load file ---------------------------------------------------------------

ctfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Claude4整理结果.xlsx"
purectdata <- read_excel(ctfile_path)
ctidfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Copy of 3_ID.xlsx"
ctiddata <- read_excel(ctidfile_path)
ctdata <- merge(purectdata, ctiddata, by = "文件夹名", all.x = TRUE)
ctdata <- ctdata[,-c(1:2,56:60)]
# 把最后一列放到第一列
ctdata <- ctdata[, c(ncol(ctdata), 1:(ncol(ctdata)-1))]

pfsfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/汇总.xlsx"
pfsfile <- read_excel(pfsfile_path)
colnames(pfsfile) <- pfsfile[1,]
pfsfile <- pfsfile[-1,]
data <- merge(ctdata, pfsfile, by = "Hospitalization_Number", all.x = TRUE)
# 将第 57, 58, 60, 61 列的数字转换为日期
data[, c(56, 57, 59, 60)] <- lapply(data[, c(56, 57, 59, 60)], function(x) as.Date(as.numeric(x), origin = "1899-12-30"))
# 在第58和59列之间插入一列，为58列减去59列的天数差
data <- data %>%
  mutate(Difference_57_56 = as.numeric(data[[57]] - data[[56]])) %>%
  relocate(Difference_57_56, .after = 57)

# 同理
data <- data %>%
  mutate(Difference_61_60 = as.numeric(data[[61]] - data[[60]])) %>%
  relocate(Difference_61_60, .after = 61)


sum(is.na(data))

# KNN填补缺失值
data_filled_ <- kNN(data, k = 5)
# 不知道有些无限接近于0的负值表现出来就是0还是生存分析不能接受0值，反正就有0的数据容易报错
data_filled_ <- data_filled_[-which(data_filled_[,58]<=0),]
sum(is.na(data_filled_))


# 因为它是是否存活和是否进展，所以要反过来，存活是0，不存活才是1
data_filled_[data_filled_[, 59] %in% "否", 59] <- 0
data_filled_[data_filled_[, 59] %in% "是", 59] <- 1
data_filled_[data_filled_[, 63] %in% "是", 63] <- 0
data_filled_[data_filled_[, 63] %in% "否", 63] <- 1


colnames(data_filled_)[58] <- "PFS"
colnames(data_filled_)[59] <- "Progression"
colnames(data_filled_)[62] <- "OS.time"
colnames(data_filled_)[63] <- "OS"



# 将所有变量转为数值型
set.seed(123)
data_filled <- data_filled_[,c(58, 59, 2:53)]


# cols <- c(3:8, 10:34)
data_filled_factor <- data_filled
data_filled_factor[, c(2,5:9)] <- lapply(data_filled_factor[, c(2,5:9)], function(x) as.numeric(x))
str(data_filled_factor)

sum(is.na(data_filled_factor))

# 筛选出最好的大方向模型之后筛选出最好的种子数
set.seed(123)
trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
trainData_factor <- data_filled_factor[trainIndex, ]
testData_factor <- data_filled_factor[-trainIndex, ]

#下面开始进行模型训练
# 设置种子数和节点数，其中节点数可以调整
rf_nodesize <- 5
seed <- 123
# 1-1.RSF -----------------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)

rs <- predict(fit, newdata = data_filled_factor)$predicted
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs, data_filled_factor))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- data.frame()
result <- rbind(result, cc)


# 1-2.RSF + CoxBoost ------------------------------------------------------


set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

pen <- optimCoxBoostPenalty(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1, 2)]), 
                            trace=TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1, 2)]), 
                stepno = cv.res$optimal.step, penalty = pen$penalty)

rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, newdata = data.matrix(data_filled_factor2[, -c(1, 2)]), newtime = data_filled_factor2[, 1],  newstatus = data_filled_factor2[, 2], type = "lp")))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')

cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result, cc)

# 1-3.RSF + Enet ----------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

x1 <- data.matrix(data_filled_factor2[, rid])
x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = data.matrix(data_filled_factor2[, -c(1, 2)]), s = fit$lambda.min)))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

# 1-4.RSF + GBM -----------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor,
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

set.seed(seed)
fit <- gbm(formula = Surv(PFS, Progression)~., data = data_filled_factor2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(PFS,Progression)~., data = data_filled_factor2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, data_filled_factor2, n.trees = best, type = 'link')))
cc <- data.frame(Cindex=as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')

cc$Model <- paste0('RSF + ', 'GBM')
result <- rbind(result, cc)


# 1-5.RSF + Lasso ---------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

x1 <- data.matrix(data_filled_factor2[, rid])
x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor2[, -c(1, 2)]), s = fit$lambda.min)))

cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)


# 1-6.RSF + plsRcox -------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = data.matrix(data_filled_factor2[, rid]), time = data_filled_factor2$PFS, status = data_filled_factor2$Progression), nt = 10, verbose = FALSE)
fit <- plsRcox(data.matrix(data_filled_factor2[, rid]), time = data_filled_factor2$PFS, event = data_filled_factor2$Progression, nt = as.numeric(cv.plsRcox.res[5]))
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = data.matrix(data_filled_factor2[, -c(1, 2)]))))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')

cc$Model <- paste0('RSF + ', 'plsRcox')
result <- rbind(result, cc)

# 1-7.RSF + Ridge ---------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor, 
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

x1 <- data.matrix(data_filled_factor2[, rid])
x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0,
                type.measure = "class")
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor2[, -c(1, 2)]), s = fit$lambda.min)))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)

# 1-8.RSF + StepCox -------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars

data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]


for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(PFS, Progression)~., data_filled_factor2), direction = direction)
  rs <- cbind(data_filled_factor2[, 1:2], RS=predict(fit, type = 'risk', newdata = data_filled_factor2))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

# 1-9.RSF + SuperPC -------------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]

data <- list(x = t(data.matrix(data_filled_factor2[, -c(1, 2)])), y = data_filled_factor2$PFS,
             censoring.status = data_filled_factor2$Progression, 
             featurenames = colnames(data_filled_factor2)[-c(1, 2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10, 
                     n.components = 3, 
                     min.features = 5, 
                     max.features = nrow(data$x), 
                     compute.fullcv = TRUE, 
                     compute.preval = TRUE)

test <- list(x = t(data.matrix(data_filled_factor2[, -c(1, 2)])), y = data_filled_factor2$PFS, censoring.status=data_filled_factor2$Progression, featurenames = colnames(data_filled_factor2)[-c(1, 2)])
ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
rr <- as.numeric(ff$v.pred)
rs <- cbind(data_filled_factor2[, 1:2], RS = rr)

cc <- data.frame(Cindex=as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1]))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'SuperPC')
result <- rbind(result, cc)


# 1-10.RSF + survival-SVM -------------------------------------------------

set.seed(seed)
fit <- rfsrc(Surv(PFS, Progression)~., data = data_filled_factor,
             ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
fit = survivalsvm(Surv(PFS, Progression)~., data= data_filled_factor2, gamma.mu = 1)
rs <- cbind(data_filled_factor2[, 1:2], RS=as.numeric(predict(fit, data_filled_factor2)$predicted))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'survival-SVM')
result <- rbind(result,cc)

# 2.Enet ------------------------------------------------------------------

pre_var <- colnames(data_filled_factor)[3:54]
x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- cbind(data_filled_factor[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = data.matrix(data_filled_factor[,3:54]), s = fit$lambda.min)))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

# 3.StepCox ---------------------------------------------------------------

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(PFS,Progression)~., data_filled_factor), direction = direction)
  rs <- cbind(data_filled_factor[, 1:2], RS = predict(fit, type = 'risk', newdata = data_filled_factor))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(PFS, Progression)~., data_filled_factor), direction = direction)
  rid <- names(coef(fit))#这里不用卡P值，迭代的结果就是可以纳入的基因
  data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
  data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
  set.seed(seed)
  pen <- optimCoxBoostPenalty(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1,2)]),
                              trace=TRUE, start.penalty = 10, parallel = T)
  cv.res <- cv.CoxBoost(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1,2)]),
                        maxstepno = 500, K = 10 , type = "verweij", penalty = pen$penalty)
  fit <- CoxBoost(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty)
  rs <-cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, newdata = data.matrix(data_filled_factor2[, -c(1, 2)]), newtime=data_filled_factor2[, 1], newstatus=data_filled_factor2[,2], type="lp")))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + CoxBoost')
  result <- rbind(result, cc)
  x1 <- data.matrix(data_filled_factor2[, rid])
  x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2, family = "cox",alpha = alpha, nfolds = 10)
    rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = data.matrix(data_filled_factor2[, -c(1, 2)]), s = fit$lambda.min)))
    cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox', '[', direction, ']', ' + Enet', '[α=', alpha, ']')
    result <- rbind(result, cc)
  }
  set.seed(seed)
  fit <- gbm(formula = Surv(PFS, Progression)~., data = data_filled_factor2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(PFS, Progression)~., data = data_filled_factor2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, data_filled_factor2, n.trees = best, type = 'link')))
  cc <- data.frame(Cindex=as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + GBM')
  result <- rbind(result, cc)
  x1 <- data.matrix(data_filled_factor2[, rid])
  x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "binomial", alpha = 1,
                  type.measure = "class")
  rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor2[, -c(1, 2)]), s = fit$lambda.min)))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Lasso')
  result <- rbind(result, cc)
  # set.seed(seed)
  
  # # 确保所有因子变量转换为数值格式
  # 
  # data_filled_factor2_numeric <- data_filled_factor2 %>%
  #   mutate(across(everything(), ~ as.numeric(as.factor(.))))
  # 
  # # 提取所需的变量
  # x <- data_filled_factor2_numeric[, rid]
  # time <- data_filled_factor2$PFS
  # status <- as.integer(data_filled_factor2$Progression)
  # 
  # # 打印变量类型以确认数据格式
  # str(x)
  # str(time)
  # str(status)
  # 
  # 
  # 
  # 使用 cv.plsRcox 进行交叉验证
  # library(rms)
  # dd <- datadist(data_filled_factor2)
  # options(datadist = "dd")
  # set.seed(seed)
  # cv.plsRcox.res <- cv.plsRcox(list(x = data_filled_factor2[, rid], time = data_filled_factor2$PFS, status = data_filled_factor2$Progression), nt = 10, verbose = FALSE)
  # 
  # 
  # fit <- plsRcox(data.matrix(data_filled_factor2[, rid]), time = data_filled_factor2$PFS,
  #                event = data_filled_factor2$Progression, nt = as.numeric(cv.plsRcox.res[5]))
  # rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = data.matrix(data_filled_factor2[, -c(1,2)]))))
  # cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  #   rownames_to_column('ID')
  # cc$Model <- paste0('StepCox', '[', direction, ']', ' + plsRcox')
  # result <- rbind(result, cc)
  x1 <- data.matrix(data_filled_factor2[, rid])
  x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
  set.seed(seed)
  fit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "binomial", alpha = 0,
                  type.measure = "class")
  rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor2[, -c(1, 2)]), s = fit$lambda.min)))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + Ridge')
  result <- rbind(result, cc)
  set.seed(seed)
  fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor2,
               ntree = 1000, nodesize = rf_nodesize, #该值建议多调整
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- cbind(data_filled_factor2[, 1:2], RS = predict(fit, newdata = data_filled_factor2)$predicted)
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + RSF')
  result <- rbind(result, cc)
  data <- list(x = t(data.matrix(data_filled_factor2[, -c(1, 2)])), y = data_filled_factor2$PFS,
               censoring.status = data_filled_factor2$Progression,
               featurenames = colnames(data_filled_factor2)[-c(1,2)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
  cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                       n.fold = 10,
                       n.components = 3,
                       min.features = (nrow(data$x)-1),
                       max.features = nrow(data$x),
                       compute.fullcv = TRUE,
                       compute.preval = TRUE)
  
  test <- list(x = t(data.matrix(data_filled_factor2[, -c(1,2)])), y = data_filled_factor2$PFS, censoring.status = data_filled_factor2$Progression, featurenames = colnames(data_filled_factor2)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rs <- cbind(data_filled_factor2[,1:2], RS = rr)
  
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + SuperPC')
  result <- rbind(result, cc)
  fit = survivalsvm(Surv(PFS,Progression)~., data = data_filled_factor2, gamma.mu = 1)
  rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, data_filled_factor2)$predicted))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']', ' + survival-SVM')
  result <- rbind(result, cc)
}


# 4-1.CoxBoost ------------------------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]), 
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[,3:54]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- cbind(data_filled_factor[,1:2], RS = as.numeric(predict(fit, newdata = data.matrix(data_filled_factor[, 3:54]), newtime = data_filled_factor[,9], newstatus = data_filled_factor[,10], type = "lp")))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result, cc)

# 4-2.CoxBoost + Enet -----------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 0.001, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 1000, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
# 检查所有系数，查看它们的分布
coef_vals <- coef(fit)
print(coef_vals)

rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
rid
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
x1 <- data.matrix(data_filled_factor[, rid])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = data.matrix(data_filled_factor2[, -c(1,2)]), s = fit$lambda.min)))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

# 4-3.CoxBoost + GBM ------------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
set.seed(seed)
fit <- gbm(formula = Surv(PFS,Progression)~., data = data_filled_factor2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(PFS,Progression)~., data = data_filled_factor2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, data_filled_factor2, n.trees = best, type = 'link')))
cc <- data.frame(Cindex= as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'GBM')
result <- rbind(result, cc)

# 4-4.CoxBoost + Lasso ----------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]

x1 <- data.matrix(data_filled_factor2[, rid])
x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor2[, -c(1,2)]), s = fit$lambda.min)))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)

# 4-5.CoxBoost + plsRcox --------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 10, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = data.matrix(data_filled_factor2[,rid]), time = data_filled_factor2$PFS, status = data_filled_factor2$Progression), nt = 10, verbose = FALSE)
fit <- plsRcox(data.matrix(data_filled_factor2[, rid]), time = data_filled_factor2$PFS, event = data_filled_factor2$Progression, nt = as.numeric(cv.plsRcox.res[5]))
rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = data.matrix(data_filled_factor2[, -c(1,2)]))))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'plsRcox')
result <- rbind(result, cc)

# 4-6.CoxBoost + Ridge ----------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
x1 <- data.matrix(data_filled_factor2[, rid])
x2 <- as.matrix(Surv(data_filled_factor2$PFS, data_filled_factor2$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0,
                type.measure = "class")
rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor2[, -c(1,2)]), s = fit$lambda.min)))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Ridge')
result <- rbind(result, cc)

# 4-7.CoxBoost + StepCox --------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(PFS,Progression)~., data_filled_factor2), direction = direction)
  rs <- cbind(data_filled_factor2[,1:2], RS = predict(fit, type = 'risk', newdata = data_filled_factor2))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

# 4-8.CoxBoost + SuperPC --------------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data <- list(x = t(data.matrix(data_filled_factor2[, -c(1,2)])), y = data_filled_factor2$PFS, censoring.status = data_filled_factor2$Progression,
             featurenames = colnames(data_filled_factor2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = (nrow(data$x)-1),
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval =TRUE)
test <- list(x=t(data.matrix(data_filled_factor2[, -c(1,2)])), y = data_filled_factor2$PFS, censoring.status = data_filled_factor2$Progression, featurenames = colnames(data_filled_factor2)[-c(1,2)])
ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
rr <- as.numeric(ff$v.pred)
rs <- cbind(data_filled_factor2[,1:2], RS = rr)

cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'SuperPC')
result <- rbind(result, cc)

# 4-9.CoxBoost + survival-SVM ---------------------------------------------

set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                            trace = TRUE, start.penalty = 10, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor[, 'PFS'], data_filled_factor[, 'Progression'], data.matrix(data_filled_factor[, 3:54]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
fit = survivalsvm(Surv(PFS, Progression)~., data = data_filled_factor2, gamma.mu = 1)
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, data_filled_factor2)$predicted))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
result <- rbind(result, cc)

# 5.plsRcox ---------------------------------------------------------------

# set.seed(seed)
# cv.plsRcox.res = cv.plsRcox(list(x = data.matrix(data_filled_factor[,pre_var]), time = data_filled_factor$PFS, status = data_filled_factor$Progression), nt = 10, verbose = FALSE)
# fit <- plsRcox(data_filled_factor[,pre_var], time = data_filled_factor$PFS, event = data_filled_factor$Progression, nt = as.numeric(cv.plsRcox.res[5]))
# rs <- cbind(data_filled_factor[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = data.matrix(data_filled_factor[, 3:54]))))
# cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS,Progression) ~ RS, rs))$concordance[1])) %>%
#   rownames_to_column('ID')
# cc$Model <- paste0('plsRcox')
# result <- rbind(result, cc)

# 6.superpc ---------------------------------------------------------------

data <- list(x = t(data.matrix(data_filled_factor[, 3:54])), y = data_filled_factor$PFS, censoring.status = data_filled_factor$Progression, featurenames = colnames(data_filled_factor)[3:54])
set.seed(seed) 
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)

test <- list(x = t(data.matrix(data_filled_factor[,3:54])), y = data_filled_factor$PFS, censoring.status = data_filled_factor$Progression, featurenames = colnames(data_filled_factor)[3:54])
ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
rr <- as.numeric(ff$v.pred)
rs <- cbind(data_filled_factor[,1:2], RS = rr)

cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result, cc)

# 7.GBM -------------------------------------------------------------------

set.seed(seed)
fit <- gbm(formula = Surv(PFS,Progression)~., data = data_filled_factor, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(PFS, Progression)~., data = data_filled_factor, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- cbind(data_filled_factor[,1:2], RS = as.numeric(predict(fit, data_filled_factor, n.trees = best, type = 'link')))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)


# 8.survivalsvm -----------------------------------------------------------

fit = survivalsvm(Surv(PFS,Progression)~., data = data_filled_factor, gamma.mu = 1)
rs <- cbind(data_filled_factor[,1:2], RS = as.numeric(predict(fit, data_filled_factor)$predicted))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival - SVM')
result <- rbind(result, cc)

# 9.Ridge -----------------------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = glmnet(x1, x2, family = "binomial", alpha = 0, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, #例文描述：10-fold cross-validation
                  family = "binomial",
                  type.measure = "class"
)

rs <- cbind(data_filled_factor[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor[, 3:54]), s = cvfit$lambda.min)))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS,Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)

# 10.Lasso ----------------------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
rs <- cbind(data_filled_factor[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = data.matrix(data_filled_factor[, 3:54]), s = fit$lambda.min)))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

# 10.1.Lasso + CoxBoost ---------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 0.1,
                type.measure = "class")
fit$lambda
fit$lambda.1se
fit$lambda.min
# 手动调整一下lambda看看有没有更合适的
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[, c('PFS', 'Progression', rid)]
set.seed(seed)
pen <- optimCoxBoostPenalty(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBoost(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBoost(data_filled_factor2[, 'PFS'], data_filled_factor2[, 'Progression'], data.matrix(data_filled_factor2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, newdata = data.matrix(data_filled_factor2[,-c(1,2)]), newtime = data_filled_factor2[,1], newstatus = data_filled_factor2[,2], type = "lp")))
cc <- data.frame(Cindex= as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

# 10.2.Lasso + GBM --------------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
fit$lambda.1se
fit$lambda
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)
fit <- gbm(formula = Surv(PFS,Progression)~., data = data_filled_factor2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(PFS,Progression)~., data = data_filled_factor2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, data_filled_factor2, n.trees = best, type = 'link')))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1]))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

# 10.3.Lasso + plsRcox （这里的cv会卡一下要手动跑）----------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
fit$lambda.1se
fit$lambda
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = data.matrix(data_filled_factor2[, rid]), time = data_filled_factor2$PFS, status = data_filled_factor2$Progression), nt = 10, verbose = FALSE)
fit <- plsRcox(data.matrix(data_filled_factor2[, rid]), time = data_filled_factor2$PFS, event = data_filled_factor2$Progression, nt = as.numeric(cv.plsRcox.res[5]))
rs <- cbind(data_filled_factor2[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = data.matrix(data_filled_factor2[,-c(1,2)]))))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

# 10.4.Lasso + RSF --------------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
fit$lambda.1se
fit$lambda
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)
fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor2,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- cbind(data_filled_factor2[, 1:2], RS = predict(fit, newdata = data_filled_factor2)$predicted)
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

# 10.5.Lasso + stepcox ----------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
fit$lambda.1se
fit$lambda
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(PFS,Progression)~., data_filled_factor2), direction = direction)
  rs <- cbind(data_filled_factor2[, 1:2], RS = predict(fit, type = 'risk', newdata = data_filled_factor2))
  cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

# 10.6.Lasso + superPC ----------------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
fit$lambda.1se
fit$lambda
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)

data <- list(x = t(data.matrix(data_filled_factor2[,-c(1,2)])), y = data_filled_factor2$PFS, censoring.status = data_filled_factor2$Progression,
             featurenames = colnames(data_filled_factor2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = (nrow(data$x)-1),#最小特征数应小于nrow(data$x)
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
test <- list(x = t(data.matrix(data_filled_factor2[,-c(1,2)])), y = data_filled_factor2$PFS, censoring.status = data_filled_factor2$Progression, featurenames = colnames(data_filled_factor2)[-c(1,2)])
ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
rr <- as.numeric(ff$v.pred)
rs <- cbind(data_filled_factor2[, 1:2], RS = rr)

cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)

# 10.7.Lasso + survival-SVM -----------------------------------------------

x1 <- data.matrix(data_filled_factor[, pre_var])
x2 <- as.matrix(Surv(data_filled_factor$PFS, data_filled_factor$Progression))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #例文描述：10-fold cross-validation
                family = "binomial", alpha = 1,
                type.measure = "class")
fit$lambda.min
fit$lambda.1se
fit$lambda
myCoefs <- coef(fit, s = fit$lambda[8]);
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
rid
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
data_filled_factor2 <- data_filled_factor[,c('PFS', 'Progression', rid)]
set.seed(seed)
fit = survivalsvm(Surv(PFS,Progression)~., data = data_filled_factor2, gamma.mu = 1)
rs <- cbind(data_filled_factor2[,1:2], RS = as.numeric(predict(fit, data_filled_factor2)$predicted))
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ RS, rs))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)


# order -------------------------------------------------------------------

result_order <- result[order(result[,2], decreasing = T),]
file_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/结果总结_连续性PFS.xlsx"

# 如果文件已存在，加载它；如果不存在，创建一个新的Workbook
if (file.exists(file_path)) {
  wb <- loadWorkbook(file_path)
} else {
  wb <- createWorkbook()
}

# 添加一个新的sheet
addWorksheet(wb, "Immune_PFS_Image_Alldata")

# 将数据框写入到新的sheet中
writeData(wb, sheet = "Immune_PFS_Image_Alldata", result_order)

# 保存Workbook
saveWorkbook(wb, file = file_path, overwrite = TRUE)
# fit <- step(coxph(Surv(PFS, Progression)~., data_filled_factor2), direction = "forward")
# pred <- predict(fit, type = 'risk', newdata = data_filled_factor2)
# auc_train_10 <- roc(data_filled_factor2[,c(1,2)], as.numeric(pred))

# #将得到的结果赋给result2变量进行操作
# result2 <- result
# 
# ###将结果的长数据转换为宽数据
# dd2 <- pivot_wider(result2, names_from = 'ID', values_from = 'Cindex') %>% as.data.frame()
# #将C指数定义为数值型
# dd2 <- as.numeric(dd2[,-1])
# #求每个模型的C指数在三个数据集的均值
# dd2$All <- apply(dd2[,2:4], 1, mean)
# #求每个模型的C指数在GEO验证集的均值
# dd2$GEO <- apply(dd2[,3:4], 1, mean)
# ###查看每个模型的C指数
# head(dd2)
# #输出C指数结果
# write.table(dd2,"output_C_index.txt", col.names = T, row.names = F, sep = "\t", quote = F)
# 
# # 这里原文仅看了GEO验证集的C指数均值，可以看到StepCox[forward] + lasso的模型组合在GEO验证集的平均C指数最高
# 
# ### StepCox[forward] + lasso构建预测模型
# fit <- step(coxph(Surv(PFS,Progression)~., est_dd), direction = "forward")
# multiCoxSum = summary(fit)
# outTab = data.frame()
# outTab = cbind(
#   coef = multiCoxSum$coefficients[,"coef"],
#   HR = multiCoxSum$conf.int[,"exp(coef)"],
#   HR.95L = multiCoxSum$conf.int[,"lower .95"],
#   HR.95H = multiCoxSum$conf.int[,"upper .95"],
#   pvalue = multiCoxSum$coefficients[,"Pr(>|z|)"])
# outTab <- as.data.frame(outTab)
# outTab = cbind(id = row.names(outTab), outTab)
# #输出多因素cox结果
# write.table(outTab, file = "output_multiCox.txt", sep = "\t", row.names = F, quote = F)
# 
# # 多因素cox并没有基因的剔除，因此用全部基因做lasso
# rid <- names(coef(fit)) # 这里不用卡P值，迭代的结果就是可以纳入的基因
# #训练集
# est_dd2 <- est_data[,c('PFS', 'Progression', rid)]
# #验证集
# val_dd_list2 <- lapply(val_data_list, function(x){x[,c('PFS', 'Progression', rid)]})
# 
# x1 <- as.matrix(est_dd2[,rid])
# x2 <- as.matrix(Surv(est_dd2$PFS, est_dd2$Progression))
# set.seed(seed)
# fit = glmnet(x1, x2, family = "binomial", alpha = 1, lambda = NULL)
# mypal <- pal_npg("nrc")(10)

# RSF_for循环  28分--------------------------------------------------------------

# 设置并行核心数量
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 1)

# 在每个节点上加载randomForest包
clusterEvalQ(cl, {
  library(randomForest)
  library(pROC) 
  library(randomForestSRC)
  library(magrittr)
  library(tibble)
  library(survival)
})
# 将数据传递到每个节点
clusterExport(cl, c("data_filled_factor", "rf_nodesize"))
# 运行并行计算
results <- parLapply(cl, 1:1000, function(seed) {
  set.seed(seed)
  fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor, 
               ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
               splitrule = 'logrank', 
               importance = T, 
               proximity = T, 
               forest = T, 
               seed = seed)
  
  set.seed(seed)
  rs <- predict(fit, newdata = data_filled_factor)$predicted
  set.seed(seed)
  model_result <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs, data_filled_factor))$concordance[1])) %>%
    rownames_to_column('ID')
  return(model_result)
  
})

stopCluster(cl)
# 计算平均AUC并筛选出平均AUC大于0.8的模型
all_Cindex <- sapply(results, function(x) x$Cindex)
all_Cindex
max(all_Cindex)
which.max(all_Cindex)
average_Cindex <- mean(unlist(all_Cindex))
average_Cindex

# 看75的那个
set.seed(464)  
fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = 464)
set.seed(464)  
rs <- predict(fit, newdata = data_filled_factor)$predicted
set.seed(464)  
cc <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs, data_filled_factor))$concordance[1])) %>%
  rownames_to_column('ID')
cc$Cindex



# # 筛选出最好的大方向模型之后筛选出最好的种子数
# trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
# trainData_factor <- data_filled_factor[trainIndex, ]
# testData_factor <- data_filled_factor[-trainIndex, ]

# 设置并行核心数量
library(parallel)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 1)

# 在每个节点上加载randomForest包
clusterEvalQ(cl, {
  library(randomForest)
  library(pROC) 
  library(randomForestSRC)
  library(magrittr)
  library(tibble)
  library(survival)
  library(caret)
})
# 将数据传递到每个节点
clusterExport(cl, c("trainData_factor", "testData_factor",  "rf_nodesize"))
# 运行并行计算
results_test <- parLapply(cl, 1:1000, function(seed) {
  
  set.seed(seed)
  fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor, 
               ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
               splitrule = 'logrank', 
               importance = T, 
               proximity = T, 
               forest = T, 
               seed = seed)
  
  set.seed(seed)
  rs_test <- predict(fit, newdata = testData_factor)$predicted
  set.seed(seed)
  model_result_test <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_test, testData_factor))$concordance[1])) %>%
    rownames_to_column('ID')
  return(model_result_test)
  
})

stopCluster(cl)
# 计算平均AUC并筛选出平均AUC大于0.8的模型
all_Cindex_test <- sapply(results_test, function(x) x$Cindex)
all_Cindex_test
max(all_Cindex_test)
which.max(all_Cindex_test)
average_Cindex_test <- mean(unlist(all_Cindex_test))
average_Cindex_test


# all_Cindex_test <- 1:119
set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = which.max(all_Cindex_test))

# 绘制ROC曲线   28分-----------------------------------------------------------------


# 
# 
# # 算一下原始数据中3，6，9, 12, 18, 24, 36个月的进展率，以决定使用哪个time point
# # 设置时间点
# time_points <- c((1*30+2*31), (3*30+3*31), (4*30+5*31), (6*30+6*31), (9*30+9*31), (12*30+12*31),  (18*30+18*31))
# # 创建一个空的结果数据框
# pfs_progression_rates <- data.frame(Time_Point = time_points, Progression_Rate = NA)
# # 循环计算各时间点的进展率
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算进展率
#   testData_factor_date <- testData_factor
#   testData_factor_date$Progression_Status <- ifelse(testData_factor_date$PFS <= current_time & testData_factor_date$Progression == 1, 1, 0)
#   progression_rate <- mean(testData_factor_date$Progression_Status, na.rm = TRUE)  # 计算进展率
#   
#   pfs_progression_rates$Progression_Rate[i] <- progression_rate
# }
# # 查看结果
# print(pfs_progression_rates)
# 
# # 计算随访数
# # 创建一个空的数据框来存储结果
# followup_counts <- data.frame(Time_Point = time_points, Patients_at_Risk = NA)
# 
# # 计算每个时间点的在随访患者数
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算在随访中的患者数
#   patients_at_risk <- sum(data_filled_factor$PFS >= current_time, na.rm = TRUE)
#   
#   followup_counts$Patients_at_Risk[i] <- patients_at_risk
# }
# 
# # 查看每个时间点的在随访患者数
# print(followup_counts)
# 
# 
# # 所以3年的进展率是最接近50%的
# 
# 
# # 汇报给robin之后
# # 把1年、2年、3年的ROC都画出来
# # 假设predictions是模型预测的风险评分
# # 使用累积风险分数
# data_filled_factor_pre <- data_filled_factor
# set.seed(which.max(all_Cindex_test))
# # 使用RSF模型生成预测结果，type = "risk" 提取风险评分
# risk_scores <- predict(fit, data_filled_factor_pre, type = "risk")$predicted
# 
# # 提取生存时间和状态
# time <- data_filled_factor_pre$PFS
# status <- data_filled_factor_pre$Progression  # 0或1表示无进展或有进展
# 
# library(timeROC)
# riskRoc <- timeROC(T = time,delta = status,
#                    marker = risk_scores,cause = 1,
#                    weighting="marginal",
#                    times = c(365,730,1095))
# 
# multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
#   library(ggsci)
#   color <- pal_lancet()(length(time))
#   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
#        col=color[1],
#        xlab=xlab,
#        ylab=ylab,main=title)
#   #如果直接plot roc对象，无法修改标题和坐标轴标签
#   for(i in 2:length(time)){
#     plot(ROC,time=time[i],add=T,col=color[i])
#   }
#   legend("bottomright",
#          legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
#          col = color,lwd = 2,
#          bty = "n",cex = cex,text.col = color
#   )
# }
# multiTimeplot(riskRoc,time = c(365,730,1095),
#               title="Time dependent ROC curve",
#               xlab="False positive rate",
#               ylab="True positive rate",
#               cex=0.7)
# riskRoc


# 上面是训练集和验证集合起来的，接下来要把训练集和验证集拆开AUC展示

risk_scores_train <- predict(fit, trainData_factor, type = "risk")$predicted

# 提取生存时间和状态
time_train <- trainData_factor$PFS
status_train <- trainData_factor$Progression  # 0或1表示无进展或有进展

riskRoc_train <- timeROC(T = time_train,delta = status_train,
                         marker = risk_scores_train,cause = 1,
                         weighting="marginal",
                         times = c(365,730,1095))


multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  # 设置绘图参数为正方形
  par(pty = "s")
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
       col=color[1],
       xlab=xlab,
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend(x = 0.33, y = 0.4,
         legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
         col = color,lwd = 1, y.intersp = 0.5,
         bty = "n",cex = cex,text.col = color
  )
}
multiTimeplot(riskRoc_train,time = c(365,730,1095),
              title="Time dependent ROC curve(Train set)",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc_train
library(timeROC)
library(boot)

# 定义计算多时间点 AUC 的函数
calc_multi_auc <- function(data, indices, times) {
  sampled_data <- data[indices, ]
  risk_scores <- predict(fit, sampled_data, type = "risk")$predicted
  time_train <- sampled_data$PFS
  status_train <- sampled_data$Progression
  
  roc <- timeROC(T = time_train, delta = status_train,
                 marker = risk_scores, cause = 1,
                 weighting = "marginal",
                 times = times)
  return(roc$AUC)
}

# 设置自助法参数
set.seed(123)
boot_results <- boot(data = trainData_factor, statistic = calc_multi_auc, R = 1000, times = c(365, 730, 1095))

# 提取置信区间
ci_365 <- boot.ci(boot_results, index = 1, type = "perc")
ci_730 <- boot.ci(boot_results, index = 2, type = "perc")
ci_1095 <- boot.ci(boot_results, index = 3, type = "perc")

# 输出置信区间
print(ci_365)
print(ci_730)
print(ci_1095)

risk_scores_test <- predict(fit, testData_factor, type = "risk")$predicted

# 提取生存时间和状态
time_test <- testData_factor$PFS
status_test <- testData_factor$Progression  # 0或1表示无进展或有进展

riskRoc_test <- timeROC(T = time_test,delta = status_test,
                        marker = risk_scores_test,cause = 1,
                        weighting="marginal",
                        times = c(365,730,1095))

multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  # 设置绘图参数为正方形
  par(pty = "s")
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
       col=color[1],
       xlab=xlab,
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend(x = 0.33, y = 0.4,
         legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
         col = color,lwd = 1, y.intersp = 0.5,
         bty = "n",cex = cex,text.col = color
  )
}
multiTimeplot(riskRoc_test,time = c(365,730,1095),
              title="Time dependent ROC curve(Validation set)",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc_test

library(timeROC)
library(boot)

# 定义计算多时间点 AUC 的函数
calc_multi_auc <- function(data, indices, times) {
  sampled_data <- data[indices, ]
  risk_scores <- predict(fit, sampled_data, type = "risk")$predicted
  time_test <- sampled_data$PFS
  status_test <- sampled_data$Progression
  
  roc <- timeROC(T = time_test, delta = status_test,
                 marker = risk_scores, cause = 1,
                 weighting = "marginal",
                 times = times)
  return(roc$AUC)
}

# 设置自助法参数
set.seed(123)
boot_results <- boot(data = testData_factor, statistic = calc_multi_auc, R = 1000, times = c(365, 730, 1095))

# 提取置信区间
ci_365 <- boot.ci(boot_results, index = 1, type = "perc")
ci_730 <- boot.ci(boot_results, index = 2, type = "perc")
ci_1095 <- boot.ci(boot_results, index = 3, type = "perc")

# 输出置信区间
print(ci_365)
print(ci_730)
print(ci_1095)
# 
# # 创建透明背景的PNG图像
# pdf("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/投稿/画图素材/roc_plot.pdf", width = 800, height = 800, bg = "transparent")
# 
# # 绘制ROC图的代码
# multiTimeplot(riskRoc_test,time = c(365,730,1095),
#               title="Time dependent ROC curve(Validation set)",
#               xlab="False positive rate",
#               ylab="True positive rate",
#               cex=0.7)
# # 关闭图形设备，保存文件
# dev.off()
# 导入和处理四川省医的additional validation    28分----------------------------------------------

ctfile_path_add <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/合作课题组数据/GPT4反馈结果.xlsx"
purectdata_add <- read_excel(ctfile_path_add)
ctidfile_path_add <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/合作课题组数据/四川省医进展数据整理.xlsx"
ctiddata_add <- read_excel(ctidfile_path_add)
ctdata_add <- merge(purectdata_add, ctiddata_add, by = "文件夹名", all.x = TRUE)
ctdata_add <- ctdata_add[,-c(1:2,55:57,60:62)]
ctdata_add
# 把最后一列放到第一列
data_add <- ctdata_add[, c((ncol(ctdata_add)-1), ncol(ctdata_add), 1:(ncol(ctdata_add)-2))]
colnames(data_add)[1] <- "PFS"
colnames(data_add)[2] <- "Progression"

sum(is.na(data_add))

# KNN填补缺失值
data_filled_add_ <- kNN(data_add, k = 5)
data_filled_add_[,1] <- data_filled_add_[,1]*30
sum(is.na(data_filled_add_))


# cols <- c(3:8, 10:34)
data_filled_factor_add <- data_filled_add_[,1:54]
str(data_filled_factor_add)

sum(is.na(data_filled_factor_add))

# 看一下训练集和训练集分别的最优种子数、Cindex和平均Cindex    28分 --------------------------------------

# # 筛选出最好的大方向模型之后筛选出最好的种子数
# trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
# trainData_factor <- data_filled_factor[trainIndex, ]
# testData_factor <- data_filled_factor[-trainIndex, ]
# 替模型训练函数

set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor,
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = which.max(all_Cindex_test))

set.seed(which.max(all_Cindex_test))
rs_train <- predict(fit, newdata = trainData_factor)$predicted
set.seed(which.max(all_Cindex_test))
model_result_train <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_train, trainData_factor))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_train$Cindex

set.seed(which.max(all_Cindex_test))
rs_test <- predict(fit, newdata = testData_factor)$predicted
set.seed(which.max(all_Cindex_test))
model_result_test <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_test, testData_factor))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_test$Cindex

set.seed(which.max(all_Cindex_test))
rs_add <- predict(fit, newdata = data_filled_factor_add)$predicted
set.seed(which.max(all_Cindex_test))
model_result_add <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_add, data_filled_factor_add))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_add$Cindex




# 1. 使用surv_cutpoint寻找最佳cut point
trainData_factor1 <- trainData_factor
trainData_factor1$rs_train <- rs_train

cutpoint <- surv_cutpoint(trainData_factor1,
                          time = "PFS", 
                          event = "Progression", 
                          # pmethod = "logrank",
                          variables = "rs_train")
trainData_factor1 <- surv_categorize(cutpoint)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_train, data = trainData_factor1)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = trainData_factor1,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for trainitional Validation Set",
           palette = c("#00AFBB", "#E7B800"))

# 1. 使用surv_cutpoint寻找最佳cut point
testData_factor1 <- testData_factor
testData_factor1$rs_test <- rs_test

cutpoint <- surv_cutpoint(testData_factor1,
                          time = "PFS", 
                          event = "Progression", 
                          # pmethod = "logrank",
                          variables = "rs_test")
testData_factor1 <- surv_categorize(cutpoint)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_test, data = testData_factor1)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = testData_factor1,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for testitional Validation Set",
           palette = c("#00AFBB", "#E7B800"))

# 确保将rs_add添加到data_filled_factor_add数据集
data_filled_factor_add1 <- data_filled_factor_add
data_filled_factor_add1$rs_add <- rs_add

cutpoint <- surv_cutpoint(data_filled_factor_add1,
                          time = "PFS", 
                          event = "Progression", 
                          # pmethod = "logrank",
                          variables = "rs_add")
data_filled_factor_add1 <- surv_categorize(cutpoint)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_add, data = data_filled_factor_add1)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = data_filled_factor_add1,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for Additional Validation Set",
           palette = c("#00AFBB", "#E7B800"))





# 筛选出最好的大方向模型之后筛选出最好的种子数

# 变量降维过程 及nomogram的绘制（已注释掉） --------------------------------------------------------------




set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = which.max(all_Cindex_test))
library(ggRandomForests)
gg_dta_depth <- gg_minimal_depth(fit)
gg_dta_vimp <- gg_minimal_vimp(fit)
plot(gg_dta_depth)
depth_threshold <- 9.7623
# 获取最小深度小于阈值的变量
selected_vars1 <- gg_dta_depth$topvars
selected_vars1

# trainData_factor_ <- trainData_factor[, c(selected_vars1, "PFS", "Progression")]
# testData_factor_ <- testData_factor[, c(selected_vars1, "PFS", "Progression")]
# data_filled_factor_add_ <- data_filled_factor_add[, c(selected_vars1, "PFS", "Progression")]

set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = which.max(all_Cindex_test))

importance_values <- var.select(fit)

varselect_df <- importance_values$varselect

# 提取变量名和重要性值
importance_df <- data.frame(
  Variable = rownames(varselect_df),
  Importance = varselect_df$vimp
)

# 按照vimp值降序排序
importance_df <- importance_df[order(-importance_df$Importance), ]

# 绘制变量重要性图
ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance from RSF Model",
       x = "Variable",
       y = "Variable Importance (vimp)")

# 交叉验证筛选出VIMP的阈值
# 定义一个函数，用于根据VIMP阈值进行选择并计算模型的性能
evaluate_vimp_threshold <- function(vimp_threshold) {
  selected_vars3 <- importance_df$Variable[importance_df$Importance > vimp_threshold]
  selected_vars2 <- intersect(selected_vars1, selected_vars3)
  trainData_factor_ <- trainData_factor[, c(selected_vars2, "PFS", "Progression")]
  testData_factor_ <- testData_factor[, c(selected_vars2, "PFS", "Progression")]
  data_filled_factor_add_ <- data_filled_factor_add[, c(selected_vars2, "PFS", "Progression")]
  fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor_,
               ntree = 1000, nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = which.max(all_Cindex_test))
  set.seed(which.max(all_Cindex_test))
  rs_test_ <- predict(fit, newdata = testData_factor_)$predicted
  testData_factor__ <- testData_factor_
  testData_factor__$rs_test_ <- rs_test_
  
  # 生存分析: testData_factor__
  cutpoint_test_ <- surv_cutpoint(testData_factor__,
                                  time = "PFS",
                                  event = "Progression",
                                  variables = "rs_test_")
  testData_factor__ <- surv_categorize(cutpoint_test_)
  surv_diff_test <- survdiff(Surv(PFS, Progression) ~ rs_test_, data = testData_factor__)
  p_value_test <- 1 - pchisq(surv_diff_test$chisq, df = length(surv_diff_test$n) - 1)
  
  # 如果低风险组的中位生存时间大于高风险组，p值小于0.05，则继续
  if (p_value_test < 0.05) {
    set.seed(which.max(all_Cindex_test))
    rs_add_ <- predict(fit, newdata = data_filled_factor_add_)$predicted
    data_filled_factor_add__ <- data_filled_factor_add_
    data_filled_factor_add__$rs_add_ <- rs_add_
    
    # 生存分析: data_filled_factor_add__
    cutpoint_add_ <- surv_cutpoint(data_filled_factor_add__,
                                   time = "PFS",
                                   event = "Progression",
                                   variables = "rs_add_")
    data_filled_factor_add__ <- surv_categorize(cutpoint_add_)
    surv_diff_add <- survdiff(Surv(PFS, Progression) ~ rs_add_, data = data_filled_factor_add__)
    p_value_add <- 1 - pchisq(surv_diff_add$chisq, df = length(surv_diff_add$n) - 1)
    
    if (p_value_add < 0.05 ) {
      p_value <- data.frame(p_value_test = p_value_test, p_value_add = p_value_add)
      return(p_value)
    } else {
      return(NULL)  # 如果条件不满足，返回NULL
    }
  } else {
    return(NULL)  # 如果条件不满足，返回NULL
  }
}

# 遍历不同的VIMP阈值并计算性能
thresholds <- seq(0, max(importance_df$Importance), by = 0.0002)
p_values <- sapply(thresholds, evaluate_vimp_threshold)

# 筛选出同时满足第一行和第二行都小于0.05且低风险组预后好于高风险组的p_values
valid_indices <- which(sapply(1:length(p_values), function(i) !is.null(p_values[[i]])))

# 从符合条件的p_values中找到第一行最小值的索引
optimal_index <- valid_indices[which.min(sapply(valid_indices, function(i) p_values[[i]]$p_value_test))]

# 根据索引获取相应的阈值
optimal_threshold <- thresholds[optimal_index]
optimal_threshold

selected_vars3 <- importance_df$Variable[importance_df$Importance > optimal_threshold]
selected_vars2 <- intersect(selected_vars3, selected_vars1)
selected_vars2
# set.seed(which.max(all_Cindex_test))
# fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor, 
#              ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
#              splitrule = 'logrank', 
#              importance = T, 
#              proximity = T, 
#              forest = T, 
#              seed = which.max(all_Cindex_test))
# 
# importance_values <- var.select(fit)
# 
# varselect_df <- importance_values$varselect
# 
# # 提取变量名和重要性值
# importance_df <- data.frame(
#   Variable = rownames(varselect_df),
#   Importance = varselect_df$vimp
# )
# 
# # 按照vimp值降序排序
# importance_df <- importance_df[order(-importance_df$Importance), ]
# importance_df[,3] <- order(importance_df$Importance, decreasing = T)
# importance_df[importance_df$Variable%in%selected_vars3,]
# selected_vars2 <- intersect(selected_vars1, importance_df[importance_df$V3<=20,1])
# # 1. 选择vimp值大于0的变量
# selected_vars2 <- importance_df$Variable[importance_df$Importance > 0]
# library(glmnet)
# trainData_factor__ <- trainData_factor_[, c(selected_vars2, "PFS", "Progression")]
# trainData_factor___ <- subset(trainData_factor__, select = -c(PFS, Progression))
# lasso_fit <- cv.glmnet(as.matrix(trainData_factor___), trainData_factor__$PFS, alpha = 1)
# # 获取Lasso回归中选中的变量
# selected_vars_lasso <- rownames(coef(lasso_fit))[coef(lasso_fit) != 0]





# 
# # 1. 选择vimp值大于0的变量
# selected_vars1 <- importance_df$Variable[importance_df$Importance > 0]
# 
# # 2. 试一下VIMP法结合最小深度法
# #最小深度法查看变量重要性
# library(ggRandomForests)
# gg_dta_depth <- gg_minimal_depth(fit)
# plot(gg_dta_depth)
# #两种方法的结合  VIMP+min_depth
# gg_dta_vimp <- gg_minimal_vimp(fit)
# plot(gg_dta_vimp)
# 
# depth_threshold <- 9.7623
# # 获取最小深度小于阈值的变量
# selected_vars2 <- gg_dta_depth$topvars
# 
# selected_vars <- intersect(selected_vars1, selected_vars2)
# # 打印出选择的变量
# print(selected_vars)

# 重新绘制细节较多的筛选后图片
# 加载必要的库
library(ggplot2)
library(dplyr)


# 初始化结果数据框
results_claude <- data.frame(
  Variable = character(),    # 变量名
  Source = character(),      # 固定字符串 "claude"
  Cindex = numeric()         # C-index
)
# 提取特征列的名称（假设从第3列开始是特征）
features_claude <- colnames(trainData_factor)[-c(1, 2)]
# 循环遍历每个特征，拟合单变量 RSF 模型
for (feature in features_claude) {
  # 动态构造公式
  formula <- as.formula(paste("Surv(PFS, Progression) ~", feature))
  
  # 拟合单变量 RSF 模型
  set.seed(which.max(all_Cindex_test))
  rsf_model <- rfsrc(
    formula,
    data = trainData_factor,
    ntree = 1000,
    nodesize = rf_nodesize,
    splitrule = "logrank",
    importance = FALSE,  # 单变量分析无需计算变量重要性
    proximity = FALSE,
    forest = TRUE,
    seed = which.max(all_Cindex_test)
  )
  # 提取 C-index
  set.seed(which.max(all_Cindex_test))
  rs <- predict(rsf_model, newdata = testData_factor[,colnames(testData_factor)%in%c("PFS", "Progression", feature)])$predicted
  set.seed(which.max(all_Cindex_test))
  c_index <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs, testData_factor[,colnames(testData_factor)%in%c("PFS", "Progression", feature)]))$concordance[1])) %>%
    rownames_to_column('ID')
  
  # 将结果添加到数据框
  results_claude <- rbind(
    results_claude,
    data.frame(
      Variable = feature,
      Source = "claude",
      Cindex = c_index
    )
  )
}


# 合并数据 (假设 gg_dta_vimp 和 results_claude 已正确读取并合并)
# gg_dta_vimp <- as.data.frame(gg_dta_vimp)
gg_dta_vimp_data <- gg_dta_vimp %>%
  dplyr::rename(Variable = names) %>% 
  left_join(results_claude %>% 
              dplyr::select(Variable, Cindex = Cindex.Cindex), by = "Variable")

# 检查是否存在空数据
gg_dta_vimp_data <- gg_dta_vimp_data %>% filter(!is.na(Cindex))

gg_dta_vimp_data_plot <- gg_dta_vimp_data
gg_dta_vimp_data_plot$Variable <- gsub("_", " ", gg_dta_vimp_data_plot$Variable)
gg_dta_vimp_data_plot$Variable[which(gg_dta_vimp_data_plot$Variable%in%"Ill defined Margin")] <- "Obscure Margin"
gg_dta_vimp_data_plot$Variable[which(gg_dta_vimp_data_plot$Variable%in%"Well defined Margin")] <- "Well-defined Margin"
gg_dta_vimp_data_plot[gg_dta_vimp_data_plot$vimp<=11,3] <- ">threshold"
gg_dta_vimp_data_plot[gg_dta_vimp_data_plot$vimp>11,3] <- "<threshold"

range(gg_dta_vimp_data_plot$Cindex)
ggplot(gg_dta_vimp_data_plot, aes(x = vimp, y = depth, size = Cindex, color = col)) +
  geom_point(alpha = 0.8) +  # 使用透明度增强可视化效果
  scale_size_continuous(name = "C-index", range = c(4, 7)) +  # 调整圆圈大小范围
  scale_color_manual(values = c(">threshold" = "#E77D72", "<threshold" = "#56BCC2"), name = "VIMP") +  # 自定义颜色
  geom_hline(yintercept = 17.5, linetype = "dashed", color = "#E77D72") +  # 水平阈值线
  geom_vline(xintercept = 11.5, linetype = "dashed", color = "#E77D72") +  # 垂直阈值线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # 添加对角线 y=x
  theme_minimal() +  # 使用简约主题
  labs(
    x = "VIMP Rank",
    y = "Minimal Depth Rank"
    # title = "VIMP Method and Minimal Depth Variable Dimensionality Reduction Plots"
  ) +
  theme(
    legend.position = "right",  # 图例位置
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.text = element_text(size = 6),  # 坐标轴字体大小
    axis.title = element_text(size = 13),  # 坐标轴标题字体大小
    panel.background = element_rect(fill = "#EBEBEB", color = NA),  # 背景改为灰色
    panel.grid = element_line(color = "white"),  # 坐标线改为白色
    plot.margin = ggplot2::margin(t = 110)  # 增大上方的边距，给文字腾出空间
  ) +
  # 限制 y 轴范围，使背景显示范围为 60 处
  coord_cartesian(ylim = c(0, 52), xlim = c(min(gg_dta_vimp_data_plot$vimp), max(gg_dta_vimp_data_plot$vimp) + 10)) +
  # 添加右侧次坐标轴，用于显示变量名
  scale_y_continuous(
    sec.axis = sec_axis(
      ~.,  # 保持与主轴相同的比例
      name = "Variable",  # 设置右侧坐标轴标题
      breaks = gg_dta_vimp_data_plot$depth,  # 使用 depth 列作为刻度
      labels = gg_dta_vimp_data_plot$Variable  # 使用 Variable 列作为标签
    )
  ) +
  theme(
    legend.key = element_rect(fill = "grey90", color = NA),  # 设置图例项的背景为灰色
  ) +
  # 添加上方坐标轴变量名
  geom_text(aes(x = vimp, y = 54, label = Variable),  # 调整 y 值使文字位于背景外部
            size = 2.5,  # 设置字体大小
            angle = 90,  # 旋转文本
            vjust = 0,  # 垂直对齐方式
            hjust = 0,  # 水平对齐方式
            inherit.aes = FALSE,  # 使 text 层独立于其它图层
            color = "black")  












# cox_data <- trainData_factor[, c("PFS", "Progression", selected_vars)]
# # nomogram的绘制
# # 构建多变量Cox模型
# final_formula <- as.formula(paste("Surv(PFS, Progression) ~", paste(selected_vars, collapse = " + ")))
# # head(cox_data)
# str(cox_data)
# 

# 
# 
# # 重新定义单因素回归变量
# univariate_results <- data.frame(
#   Variable = character(),
#   HR = numeric(),
#   CI_lower = numeric(),
#   CI_upper = numeric(),
#   p_value = numeric()
# )
# 
# # 遍历所有变量，包括风险分层
# variables <- c(setdiff(names(cox_data), c("PFS", "Progression")))
# 
# for (var in variables) {
#   formula <- as.formula(paste("Surv(PFS, Progression) ~", var))
#   model <- coxph(formula, data = cox_data)
#   summary_model <- summary(model)
#   univariate_results <- univariate_results %>%
#     add_row(
#       Variable = var,
#       HR = summary_model$conf.int[1],
#       CI_lower = summary_model$conf.int[3],
#       CI_upper = summary_model$conf.int[4],
#       p_value = summary_model$coefficients[5]
#     )
# }
# 
# # 打印单因素 Cox 结果
# print(univariate_results)
# univariate_results[univariate_results$p_value<=0.05,1]
# 
# # 过滤出合法数据
# univariate_results <- univariate_results %>%
#   filter(HR > 0 & CI_lower > 0 & CI_upper > 0)
# 
# # 再次检查
# print(univariate_results)
# rownames(univariate_results)[1] <- "Risk_Group(High_Risk)"
# univariate_results[1,1] <- "Risk_Group(High_Risk)"
# # 格式化数据
# forest_univariate_table <- cbind(
#   univariate_results$Variable,
#   paste0(sprintf("%.2f", univariate_results$HR), " (",
#          sprintf("%.2f", univariate_results$CI_lower), "-",
#          sprintf("%.2f", univariate_results$CI_upper), ")"),
#   ifelse(univariate_results$p_value < 0.0001, "<0.0001", sprintf("%.4f", univariate_results$p_value))
# )
# 
# # 绘制森林图
# 
# forestplot(
#   labeltext = forest_univariate_table,
#   mean = univariate_results$HR,
#   lower = univariate_results$CI_lower,
#   upper = univariate_results$CI_upper,
#   xlog = TRUE,
#   title = "Forest Plot of Univariate Cox Regression for Risk Score Group After Claude Matrix Variable Dimensionality Reduction",
#   zero = 1,
#   boxsize = 0.2,
#   xlab = "Hazard Ratio (log scale)"
# )

# library(rms)
# cox_rms_model <- cph(final_formula, data = cox_data, x = TRUE, y = TRUE, surv = TRUE)
# dd <- datadist(cox_data)
# options(datadist = "dd")
# med <- Quantile(cox_rms_model) #计算中位生存时间
# surv <- Survival(cox_rms_model) #构建生存概率函数
# 
# # 设置全局小数点精度为两位
# options(digits = 2)
# 
# nom <- nomogram(cox_rms_model,
#                 fun = list(function(x) surv(365, x),
#                            function(x) surv(365*2, x),
#                            function(x) surv(365*3, x)),
#                 funlabel = c("1-Year Progression Rate", "2-Year Progression Rate", "3-Year Progression Rate"),
#                 lp = T)
# plot(nom, xfrac = 0.7, cex.axis = 0.6)
# title(main = "Nomogram for Predicting Progression",
#       line = 9)
# # 恢复默认 options
# options(digits = 10)

# # 风险分层后做独立性研究 -------------------------------------------------------------
# 
# # 计算 TPR (灵敏度) 和 FPR (1-特异度)
# roc_obj <- roc(
#   response = data_filled_factor_clinic2$Progression,
#   predictor = risk_scores_filter
# )
# 
# # 手动提取灵敏度、特异度和阈值
# tpr <- roc_obj$sensitivities
# fpr <- 1 - roc_obj$specificities
# thresholds <- roc_obj$thresholds
# 
# # 计算 Youden's Index = TPR - FPR
# youden_index <- tpr - fpr
# 
# # 找到最大 Youden's Index 的阈值
# optimal_index <- which.max(youden_index)
# optimal_threshold <- thresholds[optimal_index]
# 
# # 打印最佳阈值
# print(optimal_threshold)
# 
# data_filled_factor_clinic3 <- data_filled_factor_clinic2
# # 根据阈值划分高低风险组
# data_filled_factor_clinic3$risk_group <- ifelse(
#   risk_scores_filter > optimal_threshold,
#   "High",
#   "Low"
# )
# 
# # 检查分组结果
# table(data_filled_factor_clinic3$risk_group)
# 
# 
# # 筛选高风险组数据
# high_risk_data <- data_filled_factor_clinic3[data_filled_factor_clinic3$risk_group == "High", ]
# 
# # 单因素 Cox 回归
# high_risk_results <- data.frame(
#   Variable = character(),
#   HR = numeric(),
#   CI_lower = numeric(),
#   CI_upper = numeric(),
#   p_value = numeric()
# )
# 
# variables <- names(high_risk_data)[!(names(high_risk_data) %in% c("PFS", "Progression", "risk_group"))]
# 
# for (var in variables) {
#   formula <- as.formula(paste("Surv(PFS, Progression) ~", var))
#   model <- coxph(formula, data = high_risk_data)
#   summary_model <- summary(model)
#   high_risk_results <- high_risk_results %>%
#     add_row(
#       Variable = var,
#       HR = summary_model$conf.int[1],
#       CI_lower = summary_model$conf.int[3],
#       CI_upper = summary_model$conf.int[4],
#       p_value = summary_model$coefficients[5]
#     )
# }
# 
# # 打印高风险组单因素 Cox 结果
# print(high_risk_results)
# 
# 
# # 筛选低风险组数据
# low_risk_data <- data_filled_factor_clinic3[data_filled_factor_clinic3$risk_group == "Low", ]
# 
# # 单因素 Cox 回归
# low_risk_results <- data.frame(
#   Variable = character(),
#   HR = numeric(),
#   CI_lower = numeric(),
#   CI_upper = numeric(),
#   p_value = numeric()
# )
# 
# variables <- names(low_risk_data)[!(names(low_risk_data) %in% c("PFS", "Progression", "risk_group"))]
# 
# for (var in variables) {
#   formula <- as.formula(paste("Surv(PFS, Progression) ~", var))
#   model <- coxph(formula, data = low_risk_data)
#   summary_model <- summary(model)
#   low_risk_results <- low_risk_results %>%
#     add_row(
#       Variable = var,
#       HR = summary_model$conf.int[1],
#       CI_lower = summary_model$conf.int[3],
#       CI_upper = summary_model$conf.int[4],
#       p_value = summary_model$coefficients[5]
#     )
# }
# 
# # 打印低风险组单因素 Cox 结果
# print(low_risk_results)
# 
# # 多因素 Cox 回归
# high_multivariable_formula <- as.formula(paste(
#   "Surv(PFS, Progression) ~",
#   paste(variables, collapse = " + ")
# ))
# high_multivariable_model <- coxph(high_multivariable_formula, data = high_risk_data)
# summary_high_multivariable <- summary(high_multivariable_model)
# 
# # 打印高风险组多因素 Cox 回归结果
# print(summary_high_multivariable)
# 
# # 多因素 Cox 回归
# low_multivariable_formula <- as.formula(paste(
#   "Surv(PFS, Progression) ~",
#   paste(variables, collapse = " + ")
# ))
# low_multivariable_model <- coxph(low_multivariable_formula, data = low_risk_data)
# summary_low_multivariable <- summary(low_multivariable_model)
# 
# # 打印低风险组多因素 Cox 回归结果
# print(summary_low_multivariable)
# 
# 
# 各份数据的风险分数二分类生存曲线 ----------------------------------------------------------

# set.seed(119)
# trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
# trainData_factor <- data_filled_factor[trainIndex, ]
# testData_factor <- data_filled_factor[-trainIndex, ]

# 筛选出最好的大方向模型之后筛选出最好的种子数
# data_filled_factor1 <- data_filled_factor[, c("PFS", "Progression", selected_vars)]
trainData_factor2 <- trainData_factor[, c("PFS", "Progression", selected_vars2)]
testData_factor2 <- testData_factor[, c("PFS", "Progression", selected_vars2)]
data_filled_factor_add2 <- data_filled_factor_add[, c("PFS", "Progression", selected_vars2)]

# 替模型训练函数
set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor2, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = which.max(all_Cindex_test))

set.seed(which.max(all_Cindex_test))
rs_train2 <- predict(fit, newdata = trainData_factor2)$predicted
set.seed(which.max(all_Cindex_test))
model_result_train2 <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_train2, trainData_factor2))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_train2$Cindex

set.seed(which.max(all_Cindex_test))
rs_test2 <- predict(fit, newdata = testData_factor2)$predicted
set.seed(which.max(all_Cindex_test))
model_result_test2 <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_test2, testData_factor2))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_test2$Cindex

set.seed(which.max(all_Cindex_test))
rs_add2 <- predict(fit, newdata = data_filled_factor_add2)$predicted
set.seed(which.max(all_Cindex_test))
model_result_add2 <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_add2, data_filled_factor_add2))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_add2$Cindex


# 1. 使用surv_cutpoint寻找最佳cut point
trainData_factor3 <- trainData_factor2
trainData_factor3$rs_train2 <- rs_train2

cutpoint_train1 <- surv_cutpoint(trainData_factor3,
                          time = "PFS",
                          event = "Progression",
                          # pmethod = "logrank",
                          variables = "rs_train2")
trainData_factor3 <- surv_categorize(cutpoint_train1)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_train2, data = trainData_factor3)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = trainData_factor3,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for trainitional Validation Set",
           palette = c("#00AFBB", "#E7B800"))

# 1. 使用surv_cutpoint寻找最佳cut point
testData_factor3 <- testData_factor2
testData_factor3$rs_test2 <- rs_test2

cutpoint_test1 <- surv_cutpoint(testData_factor3,
                          time = "PFS",
                          event = "Progression",
                          # pmethod = "logrank",
                          variables = "rs_test2")
testData_factor3 <- surv_categorize(cutpoint_test1)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_test2, data = testData_factor3)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = testData_factor3,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for testitional Validation Set",
           palette = c("#00AFBB", "#E7B800"))

# 确保将rs_add添加到data_filled_factor_add数据集
data_filled_factor_add3 <- data_filled_factor_add2
data_filled_factor_add3$rs_add2 <- rs_add2

cutpoint_add1 <- surv_cutpoint(data_filled_factor_add3,
                          time = "PFS", 
                          event = "Progression", 
                          # pmethod = "logrank",
                          variables = "rs_add2")
data_filled_factor_add3 <- surv_categorize(cutpoint_add1)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_add2, data = data_filled_factor_add3)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = data_filled_factor_add3,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for Additional Validation Set",
           palette = c("#00AFBB", "#E7B800"))




# 
# 
# # load("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/cox结果最好的那个data_filled_factor_clinic2.RData")
# colnames(data_filled_factor_clinic2)[13] <- "Classification_Other_NSCLC"
# data_filled_factor_clinic2[c(12, 257),13] <- 1
# data_filled_factor_clinic2[,3] <-  risk_scores_filter
# 
# # 计算 TPR (灵敏度) 和 FPR (1-特异度)
# roc_obj <- roc(
#   response = data_filled_factor_clinic2$Progression,
#   predictor = risk_scores_filter
# )
# 
# # 手动提取灵敏度、特异度和阈值
# tpr <- roc_obj$sensitivities
# fpr <- 1 - roc_obj$specificities
# thresholds <- roc_obj$thresholds
# 
# # 计算 Youden's Index = TPR - FPR
# youden_index <- tpr - fpr
# 
# # 找到最大 Youden's Index 的阈值
# optimal_index <- which.max(youden_index)
# optimal_threshold <- thresholds[optimal_index]
# 
# # 打印最佳阈值
# print(optimal_threshold)
# 
# # 创建风险分层变量
# data_filled_factor_clinic4 <- data_filled_factor_clinic2
# 
# data_filled_factor_clinic4$risk_group <- cut(
#   risk_scores_filter,
#   breaks = c(-Inf, optimal_threshold, Inf),
#   labels = c("Low", "High")
# )
# 
# # 确保新变量生成成功
# table(data_filled_factor_clinic4$risk_group)
# # 训练集=验证集的二分类生存曲线
# # 加载 survminer 包
# library(survminer)
# 
# # 创建生存对象
# surv_object <- Surv(data_filled_factor_clinic4$PFS, data_filled_factor_clinic4$Progression)
# 
# # 绘制 Kaplan-Meier 曲线
# km_fit <- survfit(surv_object ~ risk_group, data = data_filled_factor_clinic4)
# 
# # 绘制曲线
# ggsurvplot(
#   km_fit,
#   data = data_filled_factor_clinic4,
#   risk.table = TRUE,
#   pval = TRUE,  # 显示 Log-rank 检验的 p 值
#   conf.int = TRUE, 
#   legend.title = "Risk Group",
#   legend.labs = c("Low Risk", "High Risk"),
#   title = "Kaplan-Meier Curve of Claude RSF Model Risk Score Group After Claude Matrix Variable Dimensionality Reduction for Overall Set",
#   palette = c("#E7B800", "#00AFBB")
# )
#  
# # 训练集和验证集的生存曲线
# 
# set.seed(464)
# trainIndex_clinic2 <- createDataPartition(data_filled_factor_clinic2$Progression, p = 0.8, list = FALSE)
# trainData_factor_clinic2 <- data_filled_factor_clinic2[trainIndex_clinic2, ]
# validationData_factor_clinic2 <- data_filled_factor_clinic2[-trainIndex_clinic2, ]
# 
# # 计算 TPR (灵敏度) 和 FPR (1-特异度)
# roc_obj_clinic2_training <- roc(
#   response = trainData_factor_clinic2$Progression,
#   predictor = trainData_factor_clinic2$risk_scores_filter
# )
# 
# # 手动提取灵敏度、特异度和阈值
# tpr_clinic2_training <- roc_obj_clinic2_training$sensitivities
# fpr_clinic2_training <- 1 - roc_obj_clinic2_training$specificities
# thresholds_clinic2_training <- roc_obj_clinic2_training$thresholds
# 
# # 计算 Youden's Index = TPR - FPR
# youden_index_clinic2_training <- tpr_clinic2_training - fpr_clinic2_training
# 
# # 找到最大 Youden's Index 的阈值
# optimal_index_clinic2_training <- which.max(youden_index_clinic2_training)
# optimal_threshold_clinic2_training <- thresholds_clinic2_training[optimal_index_clinic2_training]
# 
# # 打印最佳阈值
# print(optimal_threshold_clinic2_training)
# 
# # 创建风险分层变量
# trainData_factor_clinic4 <- trainData_factor_clinic2
# 
# trainData_factor_clinic4$risk_group <- cut(
#   trainData_factor_clinic2$risk_scores_filter,
#   breaks = c(-Inf, optimal_threshold_clinic2_training, Inf),
#   labels = c("Low", "High")
# )
# 
# # 确保新变量生成成功
# table(trainData_factor_clinic4$risk_group)
# 
# 
# # 训练集=验证集的二分类生存曲线——训练集
# # 加载 survminer 包
# library(survminer)
# 
# # 创建生存对象
# surv_object_clinic4_training<- Surv(trainData_factor_clinic4$PFS, trainData_factor_clinic4$Progression)
# 
# # 绘制 Kaplan-Meier 曲线
# km_fit_clinic4_training <- survfit(surv_object_clinic4_training ~ risk_group, data = trainData_factor_clinic4)
# 
# # 绘制曲线
# ggsurvplot(
#   km_fit_clinic4_training,
#   data = trainData_factor_clinic4,
#   risk.table = TRUE,
#   pval = TRUE,  # 显示 Log-rank 检验的 p 值
#   conf.int = TRUE, 
#   legend.title = "Risk Group",
#   legend.labs = c("Low Risk", "High Risk"),
#   title = "Kaplan-Meier Curve of Claude RSF Model Risk Score Group After Claude Matrix Variable Dimensionality Reduction for Training Set",
#   palette = c("#E7B800", "#00AFBB")
# )
# 
# # 计算 TPR (灵敏度) 和 FPR (1-特异度)
# roc_obj_clinic2_validation <- roc(
#   response = validationData_factor_clinic2$Progression,
#   predictor = validationData_factor_clinic2$risk_scores_filter
# )
# 
# # 手动提取灵敏度、特异度和阈值
# tpr_clinic2_validation <- roc_obj_clinic2_validation$sensitivities
# fpr_clinic2_validation <- 1 - roc_obj_clinic2_validation$specificities
# thresholds_clinic2_validation <- roc_obj_clinic2_validation$thresholds
# 
# # 计算 Youden's Index = TPR - FPR
# youden_index_clinic2_validation <- tpr_clinic2_validation - fpr_clinic2_validation
# 
# # 找到最大 Youden's Index 的阈值
# optimal_index_clinic2_validation <- which.max(youden_index_clinic2_validation)
# optimal_threshold_clinic2_validation <- thresholds_clinic2_validation[optimal_index_clinic2_validation]
# 
# # 打印最佳阈值
# print(optimal_threshold_clinic2_validation)
# 
# # 创建风险分层变量
# validationData_factor_clinic4 <- validationData_factor_clinic2
# 
# validationData_factor_clinic4$risk_group <- cut(
#   validationData_factor_clinic2$risk_scores_filter,
#   breaks = c(-Inf, optimal_threshold_clinic2_validation, Inf),
#   labels = c("Low", "High")
# )
# 
# # 确保新变量生成成功
# table(validationData_factor_clinic4$risk_group)
# 
# 
# # 训练集=验证集的二分类生存曲线——验证集
# # 加载 survminer 包
# library(survminer)
# 
# # 创建生存对象
# surv_object_clinic4_validation<- Surv(validationData_factor_clinic4$PFS, validationData_factor_clinic4$Progression)
# 
# # 绘制 Kaplan-Meier 曲线
# km_fit_clinic4_validation <- survfit(surv_object_clinic4_validation ~ risk_group, data = validationData_factor_clinic4)
# 
# # 绘制曲线
# ggsurvplot(
#   km_fit_clinic4_validation,
#   data = validationData_factor_clinic4,
#   risk.table = TRUE,
#   pval = TRUE,  # 显示 Log-rank 检验的 p 值
#   conf.int = TRUE, 
#   legend.title = "Risk Group",
#   legend.labs = c("Low Risk", "High Risk"),
#   title = "Kaplan-Meier Curve of Claude RSF Model Risk Score Group After Claude Matrix Variable Dimensionality Reduction for Validation Set",
#   palette = c("#E7B800", "#00AFBB")
# )
# 
# 
# 
# 
# # 风险分数分为2分类变量的生存曲线——四川省人医的额外验证集
# data_filled_factor_add <- data_filled_add_[,c("PFS", "Progression", selected_vars)]
# addvalidationData_factor_clinic1 <- data_filled_factor_add
# 
# set.seed(which.max(all_Cindex))  
# fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor_filter, 
#               ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
#               splitrule = 'logrank', 
#               importance = T, 
#               proximity = T, 
#               forest = T, 
#               seed = which.max(all_Cindex))
# 
# risk_scores_filter <- predict(fit, addvalidationData_factor_clinic1, type = "risk")$predicted
# 
# # 合并在一起为一个新的数据框
# addvalidationData_factor_clinic2 <- cbind(addvalidationData_factor_clinic1[,1:2], risk_scores_filter)
# 
# # 计算 TPR (灵敏度) 和 FPR (1-特异度)
# roc_obj_clinic2_addvalidation <- roc(
#   response = addvalidationData_factor_clinic2$Progression,
#   predictor = addvalidationData_factor_clinic2$risk_scores_filter
# )
# 
# # 手动提取灵敏度、特异度和阈值
# tpr_clinic2_addvalidation <- roc_obj_clinic2_addvalidation$sensitivities
# fpr_clinic2_addvalidation <- 1 - roc_obj_clinic2_addvalidation$specificities
# thresholds_clinic2_addvalidation <- roc_obj_clinic2_addvalidation$thresholds
# 
# # 计算 Youden's Index = TPR - FPR
# youden_index_clinic2_addvalidation <- tpr_clinic2_addvalidation - fpr_clinic2_addvalidation
# 
# # 找到最大 Youden's Index 的阈值
# optimal_index_clinic2_addvalidation <- which.max(youden_index_clinic2_addvalidation)
# optimal_threshold_clinic2_addvalidation <- thresholds_clinic2_addvalidation[optimal_index_clinic2_addvalidation]
# 
# # 打印最佳阈值
# print(optimal_threshold_clinic2_addvalidation)
# 
# # 创建风险分层变量
# addvalidationData_factor_clinic4 <- addvalidationData_factor_clinic2
# 
# addvalidationData_factor_clinic4$risk_group <- cut(
#   addvalidationData_factor_clinic2$risk_scores_filter,
#   breaks = c(-Inf, optimal_threshold_clinic2_addvalidation, Inf),
#   labels = c("Low", "High")
# )
# 
# # 确保新变量生成成功
# table(addvalidationData_factor_clinic4$risk_group)
# 
# 
# # 训练集=验证集的二分类生存曲线
# # 加载 survminer 包
# library(survminer)
# 
# # 创建生存对象
# surv_object_clinic4_addvalidation<- Surv(addvalidationData_factor_clinic4$PFS, addvalidationData_factor_clinic4$Progression)
# 
# # 绘制 Kaplan-Meier 曲线
# km_fit_clinic4_addvalidation <- survfit(surv_object_clinic4_addvalidation ~ risk_group, data = addvalidationData_factor_clinic4)
# 
# # 绘制曲线
# ggsurvplot(
#   km_fit_clinic4_addvalidation,
#   data = addvalidationData_factor_clinic4,
#   risk.table = TRUE,
#   pval = TRUE,  # 显示 Log-rank 检验的 p 值
#   conf.int = TRUE, 
#   legend.title = "Risk Group",
#   legend.labs = c("Low Risk", "High Risk"),
#   title = "Kaplan-Meier Curve of Claude RSF Model Risk Score Group After Claude Matrix Variable Dimensionality Reduction for Additional Validation Set",
#   palette = c("#E7B800", "#00AFBB")
# )


# 验证影像学变量独立于临床变量 ----------------------------------------------------------

# 因为临床变量在传统观念上对判断肺癌预后上更有帮助
# 所以要证明这个小模型相对临床变量预测是独立的，并不受到临床变量的影像



# 导入临床数据进行处理与合并

ctfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Claude4整理结果.xlsx"
purectdata <- read_excel(ctfile_path)
ctidfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Copy of 3_ID.xlsx"
ctiddata <- read_excel(ctidfile_path)
ctdata <- merge(purectdata, ctiddata, by = "文件夹名", all.x = TRUE)
ctdata <- ctdata[,-c(1:2,56:60)]
# 把最后一列放到第一列
ctdata <- ctdata[, c(ncol(ctdata), 1:(ncol(ctdata)-1))]

pfsfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/汇总.xlsx"
pfsfile <- read_excel(pfsfile_path)
colnames(pfsfile) <- pfsfile[1,]
pfsfile <- pfsfile[-1,]
data <- merge(ctdata, pfsfile, by = "Hospitalization_Number", all.x = TRUE)
# 将第 57, 58, 60, 61 列的数字转换为日期
data[, c(56, 57, 59, 60)] <- lapply(data[, c(56, 57, 59, 60)], function(x) as.Date(as.numeric(x), origin = "1899-12-30"))
# 在第58和59列之间插入一列，为58列减去59列的天数差
data <- data %>%
  mutate(Difference_57_56 = as.numeric(data[[57]] - data[[56]])) %>%
  relocate(Difference_57_56, .after = 57)

# 同理
data <- data %>%
  mutate(Difference_61_60 = as.numeric(data[[61]] - data[[60]])) %>%
  relocate(Difference_61_60, .after = 61)


sum(is.na(data))


clinical_file_path1 <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/汇总的CT_分期分型.xlsx"
clinical_data1 <- read_excel(clinical_file_path1)
clinical_file_path2 <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/肺癌免疫治疗分期分型整理new.xlsx"
clinical_data2 <- read_excel(clinical_file_path2)
# 合并数据
merged_data1 <- merge(data, clinical_data1, by = "Hospitalization_Number", all.x = TRUE)
merged_data1$Sex[merged_data1$Sex == "男"] <- 1
merged_data1$Sex[merged_data1$Sex == "女"] <- 0
merged_data <- merge(merged_data1, clinical_data2, by = "Hospitalization_Number", all.x = TRUE)
merged_data <- merged_data[-67]

# 这些个都不是肺癌，要删掉的
# 4516511
# 5052069
# 4030867
# 4564732
# 5132389
# 2975882
# merged_data <- merged_data[!merged_data[,1]%in%c("4516511","5052069","4030867","4564732","5132389","2975882"),]


# KNN填补缺失值
data_filled_ <- kNN(merged_data, k = 5)
# 不知道有些无限接近于0的负值表现出来就是0还是生存分析不能接受0值，反正就有0的数据容易报错
data_filled_ <- data_filled_[-which(data_filled_[,58]<=0),]
sum(is.na(data_filled_))


# 因为它是是否存活和是否进展，所以要反过来，存活是0，不存活才是1
data_filled_[data_filled_[, 59] %in% "否", 59] <- 0
data_filled_[data_filled_[, 59] %in% "是", 59] <- 1
data_filled_[data_filled_[, 63] %in% "是", 63] <- 0
data_filled_[data_filled_[, 63] %in% "否", 63] <- 1


colnames(data_filled_)[58] <- "PFS"
colnames(data_filled_)[59] <- "Progression"
colnames(data_filled_)[62] <- "OS.time"
colnames(data_filled_)[63] <- "OS"



# 将所有变量转为数值型
set.seed(123)
data_filled <- data_filled_[,c(58, 59, 2:53,67,69,71:81)]

str(data_filled)

data_filled_factor_clinic <- data_filled
data_filled_factor_clinic[, c(2,55,56,60)] <- lapply(data_filled_factor_clinic[, c(2,55,56,60)], function(x) as.numeric(x))
data_filled_factor_clinic[, c(57:60)] <- lapply(data_filled_factor_clinic[, c(57:60)], function(x) as.factor(x))
str(data_filled_factor_clinic)

set.seed(123)
trainIndex <- createDataPartition(data_filled_factor_clinic$Progression, p = 0.8, list = FALSE)
trainData_factor_clinic <- data_filled_factor_clinic[trainIndex, ]
testData_factor_clinic <- data_filled_factor_clinic[-trainIndex, ]

# # 替模型训练函数
# trainData_factor_clinic1 <- trainData_factor_clinic[,c("PFS", "Progression", selected_vars2)]
# testData_factor_clinic1 <- testData_factor_clinic[,c("PFS", "Progression", selected_vars2)]
# set.seed(which.max(all_Cindex))
# fit1 <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor_clinic1,
#               ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
#               splitrule = 'logrank',
#               importance = T,
#               proximity = T,
#               forest = T,
#               seed = which.max(all_Cindex))
# 
# risk_scores_filter <- predict(fit1, data_filled_factor_filter, type = "risk")$predicted

# 合并在一起为一个新的数据框
trainData_factor_clinic1 <- cbind(trainData_factor3, trainData_factor_clinic[,55:67])
testData_factor_clinic1 <- cbind(testData_factor3, testData_factor_clinic[,55:67])

# 删去EC和SC两个变量，因为1的太少0的太多，cox做出来是NA
trainData_factor_clinic1 <- trainData_factor_clinic1[,-c(14,16)]
testData_factor_clinic1 <- testData_factor_clinic1[,-c(14,16)]
colnames(trainData_factor_clinic1)[13] <- "Classification_Other_NSCLC"
colnames(testData_factor_clinic1)[13] <- "Classification_Other_NSCLC"
colnames(trainData_factor_clinic1)[3] <- "Risk_group"
colnames(testData_factor_clinic1)[3] <- "Risk_group"
trainData_factor_clinic1[,3] <- gsub("low", "Low", trainData_factor_clinic1[,3])
trainData_factor_clinic1[,3] <- gsub("high", "High", trainData_factor_clinic1[,3])
testData_factor_clinic1[,3] <- gsub("low", "Low", testData_factor_clinic1[,3])
testData_factor_clinic1[,3] <- gsub("high", "High", testData_factor_clinic1[,3])
data_filled_factor_clinic1 <- rbind(trainData_factor_clinic1,testData_factor_clinic1)
# 单因素cox
# 存储单因素Cox回归结果
# univariate_results <- data.frame(
#   Variable = character(),
#   HR = numeric(),
#   CI_lower = numeric(),
#   CI_upper = numeric(),
#   p_value = numeric()
# )
# variables <- names(data_filled_factor_clinic2)[!(names(data_filled_factor_clinic2) %in% c("PFS", "Progression"))]
# for (var in variables) {
#   formula <- as.formula(paste("Surv(PFS, Progression) ~", var))
#   model <- coxph(formula, data = data_filled_factor_clinic2)
#   summary_model <- summary(model)
#   univariate_results <- univariate_results %>%
#     add_row(
#       Variable = var,
#       HR = summary_model$conf.int[1],
#       CI_lower = summary_model$conf.int[3],
#       CI_upper = summary_model$conf.int[4],
#       p_value = summary_model$coefficients[5]
#     )
# }
# # 打印单因素Cox结果
# print(univariate_results)
# # 多因素cox
# # 创建多因素Cox回归公式
# multivariable_formula <- as.formula(paste(
#   "Surv(PFS, Progression) ~",
#   paste(variables, collapse = " + ")
# ))
#
# # 运行多因素Cox回归
# multivariable_model <- coxph(multivariable_formula, data = data_filled_factor_clinic2)
# summary_multivariable <- summary(multivariable_model)
#
# # 打印多因素Cox回归结果
# print(summary_multivariable)
#
# # 绘制森林图
# library(forestplot)
#
# # 提取多因素Cox回归结果
# forest_data <- data.frame(
#   Variable = rownames(summary_multivariable$coefficients),
#   HR = summary_multivariable$conf.int[, 1],
#   CI_lower = summary_multivariable$conf.int[, 3],
#   CI_upper = summary_multivariable$conf.int[, 4],
#   p_value = summary_multivariable$coefficients[, 5]
# )
#
# # forest_data <- forest_data[-(4:19),]
# # 过滤出合法数据
# forest_data <- forest_data %>%
#   filter(HR > 0 & CI_lower > 0 & CI_upper > 0)
#
# # 再次检查
# print(forest_data)
#
# # 格式化数据
# forest_table <- cbind(
#   forest_data$Variable,
#   paste0(sprintf("%.2f", forest_data$HR), " (",
#          sprintf("%.2f", forest_data$CI_lower), "-",
#          sprintf("%.2f", forest_data$CI_upper), ")"),
#   sprintf("%.4f", forest_data$p_value)
# )
#
# # 绘制森林图
# forestplot(
#   labeltext = forest_table,
#   mean = forest_data$HR,
#   lower = forest_data$CI_lower,
#   upper = forest_data$CI_upper,
#   xlog = TRUE,
#   title = "Forest Plot for Multivariable Cox Regression",
#   zero = 1,
#   boxsize = 0.2,
#   xlab = "Hazard Ratio (log scale)"
# )


# 风险分层二分类变量做独立性研究 ---------------------------------------------------------

# load("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/cox结果最好的那个data_filled_factor_clinic2.RData")
# colnames(data_filled_factor_clinic2)[13] <- "Classification_Other_NSCLC"
# data_filled_factor_clinic2[c(12, 257),13] <- 1
# data_filled_factor_clinic2[,3] <-  risk_scores_filter
# 
# # 计算 TPR (灵敏度) 和 FPR (1-特异度)
# roc_obj <- roc(
#   response = data_filled_factor_clinic2$Progression,
#   predictor = risk_scores_filter
# )
# 
# # 手动提取灵敏度、特异度和阈值
# tpr <- roc_obj$sensitivities
# fpr <- 1 - roc_obj$specificities
# thresholds <- roc_obj$thresholds
# 
# # 计算 Youden's Index = TPR - FPR
# youden_index <- tpr - fpr
# 
# # 找到最大 Youden's Index 的阈值
# optimal_index <- which.max(youden_index)
# optimal_threshold <- thresholds[optimal_index]
# 
# # 打印最佳阈值
# print(optimal_threshold)
# 
# # 创建风险分层变量
# data_filled_factor_clinic4 <- data_filled_factor_clinic2
# 
# data_filled_factor_clinic4$risk_group <- cut(
#   risk_scores_filter,
#   breaks = c(-Inf, optimal_threshold, Inf),
#   labels = c("Low", "High")
# )
# 
# # 确保新变量生成成功
# table(data_filled_factor_clinic4$risk_group)
# 
# 
# # 重新定义单因素回归变量
# univariate_results <- data.frame(
#   Variable = character(),
#   HR = numeric(),
#   CI_lower = numeric(),
#   CI_upper = numeric(),
#   p_value = numeric()
# )

# Robin说换ezcox画个好看点的森林图
# install.packages("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/R_packages/forestmodel_0.6.2.tgz", repos = NULL, type = "mac.binary")
# install.packages("ezcox")

library(ezcox)
# 将列名中的 "_" 替换为空格
# T分期为0的这些个数据老几把异常了
data_filled_factor_clinic2 <- data_filled_factor_clinic1[!data_filled_factor_clinic1$T_for_TNM_Stage == "0",]
# data_filled_factor_clinic2[data_filled_factor_clinic2$Age < 60,5] <- "<60"
# data_filled_factor_clinic2[data_filled_factor_clinic2$Age >= 60,5] <- ">=60"

# 手动指定因子水平
# 设置新的水平，移除不需要的水平
data_filled_factor_clinic2$Clinical_Stage <- droplevels(data_filled_factor_clinic2$Clinical_Stage)
data_filled_factor_clinic2$T_for_TNM_Stage <- as.character(data_filled_factor_clinic2$T_for_TNM_Stage)
data_filled_factor_clinic2$T_for_TNM_Stage <- as.factor(data_filled_factor_clinic2$T_for_TNM_Stage)


# 检查修改后的因子水平
str(data_filled_factor_clinic2$`Clinical Stage`)
colnames(data_filled_factor_clinic2)[3] <- "Risk Group"
colnames(data_filled_factor_clinic2) <- gsub("_", " ", colnames(data_filled_factor_clinic2))
# 将列名中的 "_" 替换为空格
colnames(data_filled_factor_clinic2) <- gsub("_", " ", colnames(data_filled_factor_clinic2))

# 遍历所有变量，包括风险分层
variables <- c("Risk Group", setdiff(names(data_filled_factor_clinic2), c("PFS", "Progression","Risk Group")))
show_forest(
  data_filled_factor_clinic2, 
  covariates = variables, 
  # controls = "age",
  time = "PFS",
  status = "Progression"
)
show_forest(
  data_filled_factor_clinic2, 
  covariates = variables, 
  controls = "Risk Group",
  time = "PFS",
  status = "Progression"
)


# for (var in variables) {
#   formula <- as.formula(paste("Surv(PFS, Progression) ~", var))
#   model <- coxph(formula, data = trainData_factor_clinic1)
#   summary_model <- summary(model)
#   univariate_results <- univariate_results %>%
#     add_row(
#       Variable = var,
#       HR = summary_model$conf.int[1],
#       CI_lower = summary_model$conf.int[3],
#       CI_upper = summary_model$conf.int[4],
#       p_value = summary_model$coefficients[5]
#     )
# }
# 
# # 打印单因素 Cox 结果
# print(univariate_results)
# 
# 
# # 过滤出合法数据
# univariate_results <- univariate_results %>%
#   filter(HR > 0 & CI_lower > 0 & CI_upper > 0)
# 
# # 再次检查
# print(univariate_results)
# rownames(univariate_results)[1] <- "Risk_Group(High_Risk)"
# univariate_results[1,1] <- "Risk_Group(High_Risk)"
# # 格式化数据
# forest_univariate_table <- cbind(
#   univariate_results$Variable,
#   paste0(sprintf("%.2f", univariate_results$HR), " (",
#          sprintf("%.2f", univariate_results$CI_lower), "-",
#          sprintf("%.2f", univariate_results$CI_upper), ")"),
#   ifelse(univariate_results$p_value < 0.0001, "<0.0001", sprintf("%.4f", univariate_results$p_value))
# )
# 
# # 绘制森林图
# 
# forestplot(
#   labeltext = forest_univariate_table,
#   mean = univariate_results$HR,
#   lower = univariate_results$CI_lower,
#   upper = univariate_results$CI_upper,
#   xlog = TRUE,
#   title = "Forest Plot of Univariate Cox Regression for Risk Score Group After Claude Matrix Variable Dimensionality Reduction",
#   zero = 1,
#   boxsize = 0.2,
#   xlab = "Hazard Ratio (log scale)"
# )
# 
# # 定义多因素 Cox 回归公式，加入风险分层变量
# str(data_filled_factor_clinic4)
# multivariable_formula <- as.formula(paste(
#   "Surv(PFS, Progression) ~ risk_group +",
#   paste(setdiff(names(data_filled_factor_clinic4), c("PFS", "Progression", "risk_group", "risk_scores_filter")), collapse = " + ")
# ))
# 
# # 运行多因素 Cox 回归
# multivariable_model <- coxph(multivariable_formula, data = data_filled_factor_clinic4)
# summary_multivariable <- summary(multivariable_model)
# 
# # 打印多因素 Cox 结果
# print(summary_multivariable)
# 
# 
# 
# # 提取多因素Cox回归结果
# forest_data <- data.frame(
#   Variable = rownames(summary_multivariable$coefficients),
#   HR = summary_multivariable$conf.int[, 1],
#   CI_lower = summary_multivariable$conf.int[, 3],
#   CI_upper = summary_multivariable$conf.int[, 4],
#   p_value = summary_multivariable$coefficients[, 5]
# )
# 
# # forest_high_data <- forest_high_data[-(4:19),]
# # 过滤出合法数据
# forest_data <- forest_data %>%
#   filter(HR > 0 & CI_lower > 0 & CI_upper > 0)
# 
# # 再次检查
# print(forest_data)
# rownames(forest_data)[1] <- "Risk_Group(High_Risk)"
# forest_data[1,1] <- "Risk_Group(High_Risk)"
# # 格式化数据
# forest_table <- cbind(
#   forest_data$Variable,
#   paste0(sprintf("%.2f", forest_data$HR), " (",
#          sprintf("%.2f", forest_data$CI_lower), "-",
#          sprintf("%.2f", forest_data$CI_upper), ")"),
#   ifelse(univariate_results$p_value < 0.0001, "<0.0001", sprintf("%.4f", univariate_results$p_value))
# )
# 
# # 绘制森林图
# forestplot(
#   labeltext = forest_table,
#   mean = forest_data$HR,
#   lower = forest_data$CI_lower,
#   upper = forest_data$CI_upper,
#   xlog = TRUE,
#   title = "Forest Plot of Multivariable Cox Regression for Risk Score Group After Claude Matrix Variable Dimensionality Reduction",
#   zero = 1,
#   boxsize = 0.2,
#   xlab = "Hazard Ratio (log scale)"
# )

# 风险二分类的临床变量关系研究 ----------------------------------------------------------

# 加载必要的R包
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(dplyr)
library(tidyr)
data_filled_factor_clinic3 <- data_filled_factor_clinic1[!data_filled_factor_clinic1$T_for_TNM_Stage == "0",]

# colnames(trainData_factor_clinic1)[15] <- "Risk_group"
data_filled_factor_clinic3 <- data_filled_factor_clinic3 %>%
  mutate(
    Progression = factor(Progression, labels = c("No Progression", "Progression")),
    Age_Group = ifelse(Age >= 60, ">=60", "<60"),
    Sex = factor(Sex, labels = c("Male", "Female")),
    T_for_TNM_Stage = factor(T_for_TNM_Stage),
    N_for_TNM_Stage = factor(N_for_TNM_Stage),
    M_for_TNM_Stage = factor(M_for_TNM_Stage),
    Clinical_Stage = factor(Clinical_Stage),
    Classification = case_when(
      Classification_SCLC == 1 ~ "SCLC",
      Classification_SCC == 1 ~ "SCC",
      Classification_ADC == 1 ~ "ADC",
      Classification_PDC == 1 ~ "PDC",
      TRUE ~ "Other_NSCLC"
    ),
    Risk_group = factor(Risk_group, levels = c("High", "Low"))
  )

# 依据风险组排序患者
data_filled_factor_clinic3 <- data_filled_factor_clinic3 %>%
  arrange(Risk_group)
# 确保热图的患者排序一致
annotation_data <- data_filled_factor_clinic3 %>%
  dplyr::select(
    M_for_TNM_Stage,
    Clinical_Stage,
    Classification,
    T_for_TNM_Stage,
    N_for_TNM_Stage,
    Sex,
    Age_Group,
    Risk_group
  ) %>%
  mutate(Patient = rownames(data_filled_factor_clinic3)) %>%
  tidyr::pivot_longer(
    cols = -Patient,
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  mutate(
    Patient = factor(Patient, levels = rownames(data_filled_factor_clinic3)) # 确保顺序一致
  )

# 自定义颜色（与A图相近的蓝色系）
custom_palette <- c("#A6CEE3", "#1F78B4", "#33A02C", "#FB9A99", "#6A3D9A", "#E31A1C", "#FF7F00", "#B15928")

# 定义统一的值顺序（根据需求调整）
value_levels <- list(
  Risk_group = c("High", "Low"),  # red:high, blue:low
  Age_Group = c(">=60", "<60"),  # red:>=60, blue:<60
  Sex = c("Female", "Male")    # 按需调整
  # 其他变量的顺序可在此处添加
)

# 绘制热图
variables <- unique(annotation_data$Feature)
plots_list <- list()


for (var in variables) {
  data_subset <- annotation_data %>% filter(Feature == var)
  
  # 设置变量的值顺序（如果在 value_levels 中定义）
  if (var %in% names(value_levels)) {
    data_subset$Value <- factor(data_subset$Value, levels = value_levels[[var]])
  }
  
  # 获取唯一值并动态生成调色板
  unique_values <- unique(data_subset$Value)
  colors <- custom_palette[1:length(unique_values)]
  
  # 绘制单个变量的热图
  plot <- ggplot(data_subset, aes(x = Patient, y = Feature, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = setNames(colors, unique_values),
      name = var # 设置图例名称为当前变量名
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10, face = "bold"),
      legend.position = "bottom", # 图例放在下方
      legend.direction = "vertical",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    guides(fill = guide_legend(order = ifelse(var == "Age_Group", 2, 1))) # 设置图例顺序
  
  plots_list[[var]] <- plot
}

# 使用 patchwork 合并图像
B_plot_combined <- wrap_plots(plots_list, ncol = 1) + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

# 输出最终组合图
B_plot_combined


# 确保条形图的排序一致
A_plot <- ggplot(data_filled_factor_clinic3, aes(
  x = factor(row.names(data_filled_factor_clinic3), levels = rownames(data_filled_factor_clinic3)), 
  y = PFS, fill = as.factor(Progression))) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4"), name = "Progression") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 10, face = "bold")
  ) +
  labs(y = "PFS Time (days)", fill = "Progression")
A_plot
# 调整条形图和热图拼接
AB_plot <- A_plot / B_plot_combined + 
  plot_layout(heights = c(1, (length(variables)-3))) # 调整条形图和热图高度比例
AB_plot

# 卡方检验和 C 图
chisq_results <- data_filled_factor_clinic3 %>%
  dplyr::select(-PFS) %>%
  summarise(across(
    .cols = c(Age_Group, Sex, T_for_TNM_Stage, N_for_TNM_Stage, 
              M_for_TNM_Stage, Clinical_Stage, Classification),
    .fns = ~ list(chisq.test(.x, data_filled_factor_clinic3$Risk_group))
  )) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Test_Results") %>%
  mutate(
    Statistic = sapply(Test_Results, function(x) x$statistic),  # 提取卡方值
    P_Value = sapply(Test_Results, function(x) x$p.value),     # 提取 P 值
    Significant = ifelse(P_Value < 0.05, "Significant", "Not Significant")
  ) %>%
  dplyr::select(-Test_Results)
    
C_plot <- ggplot(chisq_results, aes(x = Statistic, y = reorder(Variable, Statistic), fill = Significant)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(
    values = c("Not Significant" = "#E0E0E0", "Significant" = "#1F78B4"),  # 灰色和低饱和度蓝色
    name = "Significance"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # 去掉 X 轴标题
    axis.title.y = element_blank(),  # 去掉 Y 轴标题
    axis.text.y = element_blank(),   # 去掉 Y 轴文本
    legend.position = "bottom" ,      # 图例放置底部
    legend.direction = "vertical",    # 图例纵向排布
    legend.title = element_text(size = 10, face = "bold")
  ) +
  labs(x = "Chi-Square Statistic", fill = "Significance")


C_plot

# 添加空白图例占位
empty_plot <- ggplot() + 
  theme_void() 
emptyC_plot <- empty_plot / C_plot + 
  plot_layout(heights = c(1, (length(variables)-5))) 
emptyC_plot
# 整体图形拼接
final_plot <- (AB_plot|emptyC_plot) +
  plot_layout(widths = c(5, 1)) +
  plot_annotation(title = "Patient Characteristics and Chi-square Test Results Plots")

final_plot


contingency_table <- table(data_filled_factor_clinic3$Classification, data_filled_factor_clinic3$Risk_group)
chisq_test <- chisq.test(contingency_table)

# 查看残差
print(chisq_test$residuals)

# TCGA作为补充验证集分析 ----------------------------------------------------------------

# 载入TCGA数据 ----------------------------------------------------------------


# TCGA数据
load("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/LUAD_mean.Rdata")
load("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/LUSC_mean.Rdata")
# 去掉 "patient-" 并添加 "-01A"
LUAD_Claude3.5_mean$TCGA_ID <- sub("^patient-", "", LUAD_Claude3.5_mean[[1]])
LUSC_Claude3.5_mean$TCGA_ID <- sub("^patient-", "", LUSC_Claude3.5_mean[[1]])
LUAD_Claude3.5_mean <- LUAD_Claude3.5_mean[,-1]
LUSC_Claude3.5_mean <- LUSC_Claude3.5_mean[,-1]
head(LUSC_Claude3.5_mean[,1])
str(LUSC_Claude3.5_mean)
library(readxl)
LUAD_Phenotype <- read_excel("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/LUAD.xlsx")
str(LUAD_Phenotype)
colnames(LUAD_Phenotype)[2] <- "TCGA_ID"
LUAD_pfs <- LUAD_Phenotype[,c(2, 42, 43)]
str(LUAD_pfs)
unique(LUAD_pfs[,3])
head(LUAD_pfs)
LUAD_pfs[,2] <- lapply(LUAD_pfs[,2], as.numeric)
head(LUAD_pfs)
str(LUAD_pfs)

str(LUAD_Claude3.5_mean)
LUAD_data <- merge(LUAD_Claude3.5_mean, LUAD_pfs, by = "TCGA_ID", all.x = TRUE)
str(LUAD_data)
LUAD_data <- LUAD_data[,c(1, 54, 55, 2:53)]
# LUAD_data <- kNN(LUAD_data, k = 5)
str(LUAD_data)
LUAD_data <- LUAD_data[!is.na(LUAD_data[,3]),]
LUAD_data[,2] <- LUAD_data[,2]*30
which(LUAD_data[,2]<=0)
str(LUAD_data)


LUSC_Phenotype <- read_excel("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/LUSC.xlsx")
str(LUSC_Phenotype)
colnames(LUSC_Phenotype)[2] <- "TCGA_ID"
LUSC_pfs <- LUSC_Phenotype[,c(2, 41, 42)]
str(LUSC_pfs)
unique(LUSC_pfs[,3])
str(LUSC_pfs)
# LUSC_pfs[,2] <- lapply(LUSC_pfs[,2], as.numeric)
str(LUSC_pfs)
LUSC_data <- merge(LUSC_Claude3.5_mean, LUSC_pfs, by = "TCGA_ID", all.x = TRUE)
str(LUSC_data)
# LUSC_data <- LUSC_data[,-2]
str(LUSC_data)
LUSC_data <- LUSC_data[,c(1, 54, 55, 2:53)]
str(LUSC_data)
# LUSC_data <- kNN(LUSC_data, k = 5)
LUSC_data <- LUSC_data[!is.na(LUSC_data[,3]),]
str(LUSC_data)
LUSC_data[,2] <- LUSC_data[,2]*30
which(LUSC_data[,2]<=0)
str(LUSC_data)


TCGA_data <- rbind(LUAD_data, LUSC_data)
TCGA_claude <- rbind(LUAD_Claude3.5_mean, LUSC_Claude3.5_mean)

# tcga表达数据
luad_count_data <- read_tsv("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/TCGA-LUAD.star_counts.tsv.gz")
lusc_count_data <- read_tsv("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/TCGA-LUSC.star_counts.tsv.gz")
tcga_count_data <- cbind(luad_count_data, lusc_count_data)
colnames(tcga_count_data) <- gsub("-01A", "", colnames(tcga_count_data))
rownames(tcga_count_data) <- tcga_count_data[,1]
tcga_count_data <- tcga_count_data[,colnames(tcga_count_data)%in%TCGA_data$TCGA_ID]
tcga_count_data <- cbind(rownames(tcga_count_data), tcga_count_data)
colnames(tcga_count_data)[1] <- "Ensembl_ID"

probeMap <- read.table("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/gencode.v36.annotation.gtf.gene.probemap", header = T)
head(probeMap)
# 假设您的probeMap包含Ensembl_ID和Gene_Name等信息
# 根据Ensembl_ID进行合并
tcga_count_data <- merge(tcga_count_data, probeMap[,1:2], by.x = "Ensembl_ID", by.y = "id", all.x = TRUE)

# 查看合并后的数据
head(tcga_count_data)
tcga_count_data <- tcga_count_data %>%
  mutate(rowmean = rowMeans(.[,-c(1,61)]), .before = 2) %>% # "."代表传入的数据，计算除Ensembl_ID和gene外的每一行的均值
  arrange(desc(rowmean)) %>% # 按rowmean降序排列,重复基因取平均值最大的那个
  distinct(gene, .keep_all = T) %>% # 保留唯一值
  dplyr::select(-rowmean) %>% # 去除rowmean列
  column_to_rownames(var = "gene") # 将symbol列变成行名

tcga_count_data <- tcga_count_data[,-1]

tcga_count_data <- 2^tcga_count_data - 1 # 去log化

# tpm数据
luad_tpm_data <- read_tsv("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/TCGA-LUAD.star_tpm.tsv.gz")
lusc_tpm_data <- read_tsv("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/TCGA-LUSC.star_tpm.tsv.gz")
tcga_tpm_data <- cbind(luad_tpm_data, lusc_tpm_data)
colnames(tcga_tpm_data) <- gsub("-01A", "", colnames(tcga_tpm_data))
rownames(tcga_tpm_data) <- tcga_tpm_data[,1]
tcga_tpm_data <- tcga_tpm_data[,colnames(tcga_tpm_data)%in%TCGA_data$TCGA_ID]
tcga_tpm_data <- cbind(rownames(tcga_tpm_data), tcga_tpm_data)
colnames(tcga_tpm_data)[1] <- "Ensembl_ID"

probeMap <- read.table("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/gencode.v36.annotation.gtf.gene.probemap", header = T)
head(probeMap)
# 假设您的probeMap包含Ensembl_ID和Gene_Name等信息
# 根据Ensembl_ID进行合并
tcga_tpm_data <- merge(tcga_tpm_data, probeMap[,1:2], by.x = "Ensembl_ID", by.y = "id", all.x = TRUE)

# 查看合并后的数据
head(tcga_tpm_data)
tcga_tpm_data <- tcga_tpm_data %>%
  mutate(rowmean = rowMeans(.[,-c(1,61)]), .before = 2) %>% # "."代表传入的数据，计算除Ensembl_ID和gene外的每一行的均值
  arrange(desc(rowmean)) %>% # 按rowmean降序排列,重复基因取平均值最大的那个
  distinct(gene, .keep_all = T) %>% # 保留唯一值
  dplyr::select(-rowmean) %>% # 去除rowmean列
  column_to_rownames(var = "gene") # 将symbol列变成行名

tcga_tpm_data <- tcga_tpm_data[,-1]





# TCGA生存曲线 --------------------------------------------------------------------

set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = which.max(all_Cindex_test))

set.seed(which.max(all_Cindex_test))
rs_tcga <- predict(fit, newdata = TCGA_data)$predicted
set.seed(which.max(all_Cindex_test))
model_result_tcga <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_tcga, TCGA_data))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_tcga$Cindex

# 确保将rs_add添加到data_filled_factor_add数据集
TCGA_data1 <- TCGA_data
TCGA_data1$rs_tcga <- rs_tcga

cutpoint <- surv_cutpoint(TCGA_data1,
                          time = "PFS", 
                          event = "Progression", 
                          # pmethod = "logrank",
                          variables = "rs_tcga")
TCGA_data1 <- surv_categorize(cutpoint)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_tcga, data = TCGA_data1)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = TCGA_data1,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for Additional Validation Set",
           palette = c("#00AFBB", "#E7B800"))


# 变量筛选
TCGA_data2 <- TCGA_data[, c("TCGA_ID", "PFS", "Progression", selected_vars2)]

# 替模型训练函数
set.seed(which.max(all_Cindex_test))
fit <- rfsrc(Surv(PFS,Progression)~., data = trainData_factor2, 
             ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = which.max(all_Cindex_test))

set.seed(which.max(all_Cindex_test))
rs_tcga2 <- predict(fit, newdata = TCGA_data2)$predicted
set.seed(which.max(all_Cindex_test))
model_result_tcga2 <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs_tcga2, TCGA_data2))$concordance[1])) %>%
  rownames_to_column('ID')
model_result_tcga2$Cindex



# 1. 使用surv_cutpoint寻找最佳cut point
TCGA_data3 <- TCGA_data2
TCGA_data3$rs_tcga2 <- rs_tcga2

cutpoint_train1 <- surv_cutpoint(TCGA_data3,
                                 time = "PFS",
                                 event = "Progression",
                                 # pmethod = "logrank",
                                 variables = "rs_tcga2")
TCGA_data3 <- surv_categorize(cutpoint_train1)

# 2. 拟合生存曲线
fit <- survfit(Surv(PFS, Progression) ~ rs_tcga2, data = TCGA_data3)

# 3. 绘制生存曲线
ggsurvplot(fit,
           data = TCGA_data3,
           risk.table = TRUE, 
           pval = TRUE, 
           # conf.int = TRUE, 
           # test.for.trend = TRUE,
           legend.title = "Risk Group", 
           legend.labs = c("High Risk", "Low Risk"),
           title = "Kaplan-Meier Curve of Claude RSF Model Risk Scores for trainitional Validation Set",
           palette = c("#00AFBB", "#E7B800"))

# HR值meta分析 ---------------------------------------------------------------

# 首先是降维前的数据
colnames(trainData_factor1)[3] <- "Risk Group"
colnames(testData_factor1)[3] <- "Risk Group"
colnames(data_filled_factor_add1)[3] <- "Risk Group"
colnames(TCGA_data1)[3] <- "Risk Group"
# 然后是降维后的数据
colnames(trainData_factor3)[3] <- "Risk Group"
colnames(testData_factor3)[3] <- "Risk Group"
colnames(data_filled_factor_add3)[3] <- "Risk Group"
colnames(TCGA_data3)[3] <- "Risk Group"

trainData_meta <- cbind(trainData_factor1, trainData_factor3)
trainData_meta <- trainData_meta[,-c(4:5)]
colnames(trainData_meta)[c(3,4)] <- c("pre", "after")
testData_meta <- cbind(testData_factor1, testData_factor3)
testData_meta <- testData_meta[,-c(4:5)]
colnames(testData_meta)[c(3,4)] <- c("pre", "after")
data_filled_factor_add_meta <- cbind(data_filled_factor_add1, data_filled_factor_add3)
data_filled_factor_add_meta <- data_filled_factor_add_meta[,-c(4:5)]
colnames(data_filled_factor_add_meta)[c(3,4)] <- c("pre", "after")
TCGA_data_meta <- cbind(TCGA_data1, TCGA_data3)
TCGA_data_meta <- TCGA_data_meta[,-c(4:5)]
colnames(TCGA_data_meta)[c(3,4)] <- c("pre", "after")


# 进行Cox回归分析
train_cox_pre <- coxph(Surv(PFS, Progression) ~ pre, data = trainData_meta)
train_cox_after <- coxph(Surv(PFS, Progression) ~ after, data = trainData_meta)
test_cox_pre <- coxph(Surv(PFS, Progression) ~ pre, data = testData_meta)
test_cox_after <- coxph(Surv(PFS, Progression) ~ after, data = testData_meta)
data_filled_factor_add_cox_pre <- coxph(Surv(PFS, Progression) ~ pre, data = data_filled_factor_add_meta)
data_filled_factor_add_cox_after <- coxph(Surv(PFS, Progression) ~ after, data = data_filled_factor_add_meta)
TCGA_cox_pre <- coxph(Surv(PFS, Progression) ~ pre, data = TCGA_data_meta)
TCGA_cox_after <- coxph(Surv(PFS, Progression) ~ after, data = TCGA_data_meta)


# 提取pre和after的HR值及置信区间
extract_cox_results <- function(cox_model) {
  hr <- summary(cox_model)$conf.int[1]   # HR值
  ci_lower <- summary(cox_model)$conf.int[3]  # 置信区间下限
  ci_upper <- summary(cox_model)$conf.int[4]  # 置信区间上限
  return(data.frame(HR = hr, CI_lower = ci_lower, CI_upper = ci_upper))
}

# 计算HR值及置信区间的pre和after结果
result_pre_train <- extract_cox_results(train_cox_pre)
result_after_train <- extract_cox_results(train_cox_after)
result_pre_test <- extract_cox_results(test_cox_pre)
result_after_test <- extract_cox_results(test_cox_after)
result_pre_data_filled_factor_add <- extract_cox_results(data_filled_factor_add_cox_pre)
result_after_data_filled_factor_add <- extract_cox_results(data_filled_factor_add_cox_after)
result_pre_TCGA <- extract_cox_results(TCGA_cox_pre)
result_after_TCGA <- extract_cox_results(TCGA_cox_after)

# 创建用于pre和after的meta数据
meta_data_pre <- data.frame(
  study = c(
    # "Training set",
    "ZJH-SMU", "SAMSPH", "TCIA"),  # 研究名称
  HR = c(
    # result_pre_train$HR, 
    result_pre_test$HR, result_pre_data_filled_factor_add$HR, result_pre_TCGA$HR),  # HR值
  LCI = c(
    # result_pre_train$CI_lower, 
    result_pre_test$CI_lower, result_pre_data_filled_factor_add$CI_lower, result_pre_TCGA$CI_lower),  # 置信区间下限
  UCI = c(
    # result_pre_train$CI_upper, 
    result_pre_test$CI_upper, result_pre_data_filled_factor_add$CI_upper, result_pre_TCGA$CI_upper)   # 置信区间上限
)

meta_data_after <- data.frame(
  study = c(
    # "Training set", 
    "ZJH-SMU", "SAMSPH", "TCIA"),  # 研究名称
  HR = c(
    # result_after_train$HR, 
    result_after_test$HR, result_after_data_filled_factor_add$HR, result_after_TCGA$HR),  # HR值
  LCI = c(
    # result_after_train$CI_lower, 
    result_after_test$CI_lower, result_after_data_filled_factor_add$CI_lower, result_after_TCGA$CI_lower),  # 置信区间下限
  UCI = c(
    # result_after_train$CI_upper,
    result_after_test$CI_upper, result_after_data_filled_factor_add$CI_upper, result_after_TCGA$CI_upper)   # 置信区间上限
)

# 进行meta分析：pre和after分别进行分析
meta_analysis_pre <- metagen(log(HR), lower = log(LCI), upper = log(UCI), studlab = study, data = meta_data_pre, sm = "HR", 
                             common = TRUE,  # 使用固定效应模型
                             random = TRUE)  # 使用随机效应模型

meta_analysis_after <- metagen(log(HR), lower = log(LCI), upper = log(UCI), studlab = study, data = meta_data_after, sm = "HR", 
                               common = TRUE,  # 使用固定效应模型
                               random = TRUE)  # 使用随机效应模型

# 查看meta分析结果：pre和after的meta分析
summary(meta_analysis_pre)
summary(meta_analysis_after)


meta::forest(meta_analysis_pre, main = "Pre HR Meta-analysis")
meta::forest(meta_analysis_after, main = "After HR Meta-analysis")


# 
# # 进行pre和after的总体比较（可能涉及到两者的差异分析）
# # 可以直接在summary中对比固定效应模型和随机效应模型的结果
# 
# # 提取meta_analysis_pre和meta_analysis_after的HR值、CI和p值
# extract_meta_results <- function(meta_analysis) {
#   hr_fixed <- exp(meta_analysis$TE.fixed)  # 固定效应模型的HR值
#   ci_fixed_lower <- exp(meta_analysis$lower.fixed)  # 固定效应模型的CI下限
#   ci_fixed_upper <- exp(meta_analysis$upper.fixed)  # 固定效应模型的CI上限
#   p_fixed <- 2 * pnorm(-abs(meta_analysis$zval.fixed))  # 固定效应模型的p值
# 
#   hr_random <- exp(meta_analysis$TE.random)  # 随机效应模型的HR值
#   ci_random_lower <- exp(meta_analysis$lower.random)  # 随机效应模型的CI下限
#   ci_random_upper <- exp(meta_analysis$upper.random)  # 随机效应模型的CI上限
#   p_random <- 2 * pnorm(-abs(meta_analysis$zval.random))  # 随机效应模型的p值
# 
#   return(data.frame(HR_fixed = hr_fixed, CI_fixed_lower = ci_fixed_lower, CI_fixed_upper = ci_fixed_upper, P_fixed = p_fixed,
#                     HR_random = hr_random, CI_random_lower = ci_random_lower, CI_random_upper = ci_random_upper, P_random = p_random))
# }
# 
# 
# 
# # 提取结果
# results_pre <- extract_meta_results(meta_analysis_pre)
# results_after <- extract_meta_results(meta_analysis_after)
# 
# 
# # 创建森林图数据
# forest_data <- data.frame(
#   Study = c("Pre", "After"),
#   No_of_Studies = c(4, 4),  # 研究数量
#   HR_fixed = c(results_pre$HR_fixed, results_after$HR_fixed),
#   CI_fixed_lower = c(results_pre$CI_fixed_lower, results_after$CI_fixed_lower),
#   CI_fixed_upper = c(results_pre$CI_fixed_upper, results_after$CI_fixed_upper),
#   P_fixed = c(results_pre$P_fixed, results_after$P_fixed),
#   HR_random = c(results_pre$HR_random, results_after$HR_random),
#   CI_random_lower = c(results_pre$CI_random_lower, results_after$CI_random_lower),
#   CI_random_upper = c(results_pre$CI_random_upper, results_after$HR_random),
#   P_random = c(results_pre$P_random, results_after$P_random)
# )
# 
# 
# #设置森林图的主题，这一步很大程度上决定后面森林图好不好看
# tm <- forest_theme(core = list(bg_params=list(fill = c("white"))),#森林图背景填充
#                    base_size = 10,#主题字符大小
#                    base_family = "Arial",#字体类型
#                    summary_col = "#377EB8",#合并置信区间（菱形）的颜色
#                    refline_col = "black",#无效线颜色
#                    refline_lwd =2,#无效线粗细
#                    ci_col = "#F781BF",#置信区间颜色
#                    arrow_fill ="#F781BF",#箭头填充颜色
#                    arrow_lwd =2,#箭头粗细
#                    ci_lwd =2,#置信区间粗细
#                    ci_Theight =unit(0.2,'cm'),#置信区间两侧竖线高度
#                    arrow_label_just = "end",#箭头文本与箭头对齐方式
#                    arrow_type = "closed")#箭头类型
# 
# # 使用forestploter绘制森林图
# p <- forestploter::forest(forest_data[,1:5],
#                      est = forest_data$HR_fixed,  # 使用固定效应模型的HR值
#                      lower = forest_data$CI_fixed_lower,  # 固定效应模型的置信区间下限
#                      upper = forest_data$CI_fixed_upper,  # 固定效应模型的置信区间上限
#                      # pval = forest_data$P_fixed,  # 固定效应模型的p值
#                      studlab = forest_data$Study,  # 研究名称
#                      weight = forest_data$No_of_Studies,  # 研究数量
#                      # random = TRUE,  # 显示随机效应模型
#                      # random.est = forest_data$HR_random,  # 随机效应模型的HR值
#                      # random.ci.lb = forest_data$CI_random_lower,  # 随机效应模型的置信区间下限
#                      # random.ci.ub = forest_data$CI_random_upper,  # 随机效应模型的置信区间上限
#                      # random.pval = forest_data$P_random,  # 随机效应模型的p值
#                      ci_column = 3,             # 在第三列（随机效应模型）绘制森林图
#                      xlab = "Hazard Ratio (HR)",  # x轴标签
#                      sizes = c(4,4),# 小球球大小
#                      ref_line = 1, # 参考线
#                      xlim = c(-1, 1),
#                      ticks_at = c(0.1, 1, 10, 100), # -1、0、1
#                      arrow_lab = c("Good prognosis in the low-risk group","Bad prognosis in the low-risk group"), #方向键
#                      # columns = c("Study", "No_of_Studies", "HR_fixed", "CI_fixed_lower", "CI_fixed_upper", "P_fixed", "HR_random", "CI_random_lower", "CI_random_upper", "P_random"),  # 指定列的顺序
#                      # header = c("Study", "No. of Studies", "Fixed-effect model", "Fixed-effect model", "Random-effects model", "Random-effects model"),  # 自定义表头
#                      theme = tm,   # 设置脚注文字颜色为蓝色
#                      grid = TRUE)   # 添加网格线
# p
# 
# p <- forest(forest_data[,1:5],
#             est = forest_data$HR_fixed,
#             lower = forest_data$CI_fixed_lower, 
#             upper = forest_data$CI_fixed_upper,
#             sizes = c(4,4),
#             # is_summary = c(rep(F, nrow(dt)-1), T),
#             ci_column = 3,
#             ref_line = 0,
#             # x_trans = "log",
#             arrow_lab = c("Favours caffeine","Favours decaf"),
#             xlim = c(-1, 1),
#             ticks_at = c(-1, 0, 1),
#             theme = tm)
# p
# # 特定列改变背景色（添加底纹）
# g <- edit_plot(p, col = 1:10, 
#                row =1, 
#                which = "background", 
#                gp = gpar(fill = "#D9D9D9",alpha=0.8))
# g
# 
# # 图上方CI列上方添加文本
# g <- add_text(g, text = "IV, Random, 95% CI",
#               part = "header", 
#               col = 7:8,
#               gp = gpar(fontface = "bold"))
# 
# g <- insert_text(g, text = "Odds ratio",
#                  part = "header", 
#                  col = 7:8,
#                  gp = gpar(fontface = "bold"))
# 
# # 图上方添加分组信息
# g <- add_text(g, text = "Caffeine",
#               part = "header",
#               row = 1,
#               col = 2:3,
#               gp = gpar(fontface = "bold"))
# 
# g <- add_text(g, text = "Decaf",
#               part = "header", 
#               row = 1,
#               col = 4:5,
#               gp = gpar(fontface = "bold"))
# 
# # 红绿灯图上方添加文本
# g <- add_text(g, text = "Risk of Bias",
#               part = "header", 
#               row = 1,
#               col = 9:14,
#               gp = gpar(fontface = "bold"))
# 
# # 图下方添加结局事件
# g <- insert_text(g, 
#                  text = c("Total events:"),
#                  row = 9,
#                  col = 1,
#                  before = FALSE,
#                  just = "left")
# 
# # Note: The row counts need to add one to account for 
# # `insert_text` in the previous step
# g <- add_text(g, text = "58",
#               col = 2,
#               row = 10,
#               just = "left")
# 
# g <- add_text(g, text = "46",
#               col = 4,
#               row = 10,
#               just = "left")
# 
# g

# TCGA差异基因分析和火山图绘制 --------------------------------------------------------------------

TCGA_data4 <- cbind(TCGA_data2$TCGA_ID, TCGA_data3)
colnames(TCGA_data4)[4] <- "Risk_Group"
str(TCGA_data4)
TCGA_data4$Risk_Group <- as.factor(TCGA_data4$`Risk_Group`)
str(TCGA_data4)
# 少了两个，如果结果不好记得回来这里看能不能补上这两个少了的患者
intersect(colnames(tcga_count_data), TCGA_data2$TCGA_ID)
common_ids <- intersect(colnames(tcga_count_data), TCGA_data4$`TCGA_data2$TCGA_ID`)
TCGA_data4 <- TCGA_data4[TCGA_data4$`TCGA_data2$TCGA_ID` %in% common_ids, ]
# 创建 DESeq2 数据集对象
str(tcga_count_data)
# 转换 tcga_count_data 中的所有值为整数
tcga_count_data <- round(tcga_count_data)

# 再次尝试创建 DESeq2 数据集对象
str(TCGA_data4)
dds <- DESeqDataSetFromMatrix(countData = tcga_count_data,
                              colData = TCGA_data4,
                              design = ~ `Risk_Group`)

# # 过滤掉低表达基因（例如，表达量在大部分样本中为 0）
# dds <- dds[rowSums(counts(dds)) > 1, ]

# 差异基因分析
dds <- DESeq(dds)

# 获取差异基因分析结果
res_gene <- results(dds, contrast = c("Risk_Group", "high", "low"))
# 去除 res_gene 中差异为 NA 的基因
res_gene <- na.omit(res_gene)

# 查看清理后的结果
head(res_gene)

# 查看结果摘要
summary(res_gene)

# # 筛选显著差异基因（p值 < 0.05，绝对log2FoldChange > 1）
res_gene_sig <- res_gene[which(res_gene$pvalue < 0.05 & abs(res_gene$log2FoldChange) > 1), ]
# head(res_gene_sig)


# 排序显著基因，先按 log2FoldChange 排序（从大到小，即上调基因）
res_sig_upregulated <- res_gene_sig[order(res_gene_sig$log2FoldChange, decreasing = TRUE), ]
# 提取前十个上调基因
uptop <- rownames(head(res_sig_upregulated, 10))

# 排序显著基因，按 log2FoldChange 排序（从小到大，即下调基因）
res_sig_downregulated <- res_gene_sig[order(res_gene_sig$log2FoldChange, decreasing = FALSE), ]
# 提取前十个下调基因
dwtop <- rownames(head(res_sig_downregulated, 10))

res_gene$label <- ifelse(rownames(res_gene) %in% c(uptop,dwtop), rownames(res_gene), "")
res_gene$log10FDR <- -log10(res_gene$padj)
# 添加 DEG 列，标记上调、下调和无显著差异的基因
res_gene$DEG <- ifelse(res_gene$pvalue < 0.05 & res_gene$log2FoldChange > 1, "Up", 
                           ifelse(res_gene$pvalue < 0.05 & res_gene$log2FoldChange < -1, "Down", "None"))
colnames(res_gene)[2] <- "logFC"
# 查看前几行，确认新列添加成功
head(res_gene)
a <- res_gene[res_gene$pvalue<=0.05,5:6]
ggplot(res_gene, aes(x = logFC, y = log10FDR, colour = DEG)) +
  geom_point(alpha = 0.85, size = 1.5) +  # 设置点的透明度和大小
  scale_color_manual(values = c('steelblue', 'gray', 'brown')) +  # 调整点的颜色
  xlim(c(-11, 11)) +  # 调整 x 轴的取值范围，根据最大值确定，再进行取舍
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +  # 添加 x 轴辅助线，lty 函数调整线的类型
  geom_hline(yintercept = -log10(0.555), lty = 4, col = "black", lwd = 0.8) +  # 添加 y 轴辅助线
  labs(x = "logFC", y = "-log10FDR") +  # x、y 轴标签
  ggtitle("TCGA LIHC DEG") +  # 图表标题自定义
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.title = element_blank()) +  # 设置图表标题和图例位置
  geom_label_repel(data = res_gene, aes(label = label),  # 添加标签
                   size = 3, box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black",
                   show.legend = FALSE, max.overlaps = 10000) +  # 标签设置
  theme_prism(border = TRUE)  # 调整主题，设置边框显示


# TCGA通路分析和GSEA图绘制 ----------------------------------------------------------------

# 将DESeq2结果转换为数据框
gene_list <- res_gene_sig$log2FoldChange
names(gene_list) <- rownames(res_gene_sig)

# 排序基因列表（按log2FoldChange排序）
gene_list <- sort(gene_list, decreasing = TRUE)


# 读取基因集文件
gene_hallmark_sets <- read.gmt("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/h.all.v2024.1.Hs.symbols.gmt.txt")
gene_c2_kegg_sets <- read.gmt("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt.txt")
gene_c5_GO_sets <- read.gmt("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/c5.go.v2024.1.Hs.symbols.gmt")
gene_c6_oncogene_sets <- read.gmt("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/c6.all.v2024.1.Hs.symbols.gmt.txt")
gene_c7_immune_sets <- read.gmt("/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/CT_immune_little_test/R_test/TCGA/c7.all.v2024.1.Hs.symbols.gmt.txt")

# 使用clusterProfiler进行GSEA
gsea_hallmark_results <- GSEA(gene_list, TERM2GENE = gene_hallmark_sets, pvalueCutoff = 0.05, verbose = TRUE)
# 查看结果
summary(gsea_hallmark_results)
gsea_c2_kegg_results <- GSEA(gene_list, TERM2GENE = gene_c2_kegg_sets, pvalueCutoff = 0.05, verbose = TRUE)
summary(gsea_c2_kegg_results)
gsea_c5_GO_results <- GSEA(gene_list, TERM2GENE = gene_c5_GO_sets, pvalueCutoff = 0.05, verbose = TRUE)
summary(gsea_c5_GO_results)
gsea_c6_oncogene_results <- GSEA(gene_list, TERM2GENE = gene_c6_oncogene_sets, pvalueCutoff = 0.05, verbose = TRUE)
summary(gsea_c6_oncogene_results)
gsea_c7_immune_results <- GSEA(gene_list, TERM2GENE = gene_c7_immune_sets, pvalueCutoff = 0.05, verbose = TRUE)
summary(gsea_c7_immune_results)


# gsea_c5_GO_results <- gsea_c5_GO_results[c(3,6,27,30,32,33),]
# 看样子只有GO和oncogene能用
gseaplot2(gsea_c5_GO_results, geneSetID = c(1,7,24,26,30,33),
          pvalue_table = TRUE,
          ES_geom = "line")#"dot"将线转换为

# gseaplot2(gsea_c6_oncogene_results, geneSetID = 1:3,
#           pvalue_table = TRUE,
#           ES_geom = "line")#"dot"将线转换为




# # 假设你已经读取了count数据，并存储为tcga_count_data，行名为基因名，列名为TCGA_ID
# # 同时，假设你的TCGA_data4包含了TCGA_ID和Risk Group
# 
# # 设置分组信息
# group <- factor(TCGA_data4$`Risk Group`, levels = c("low", "high"))
# 
# # 构建DGEList对象
# dge <- DGEList(counts = tcga_count_data, group = group)
# 
# # 过滤低表达基因（可以设置合适的阈值）
# dge <- dge[filterByExpr(dge), ]
# 
# # 对count数据进行归一化
# dge <- calcNormFactors(dge)
# 
# # 差异分析
# design <- model.matrix(~0+factor(group))
# rownames(design)<-colnames(dge)
# colnames(design)<-levels(factor(group))
# 
# dge <- estimateGLMCommonDisp(dge,design)
# dge <- estimateGLMTrendedDisp(dge, design)
# dge <- estimateGLMTagwiseDisp(dge, design)
# fit <- glmFit(dge, design) 
# can.vs.nor <- makeContrasts(high-low, levels=design)
# lrt <- glmLRT(fit,contrast=can.vs.nor) #因为是两组所以这么做，按照说明书
# nrDEG=topTags(lrt, n=nrow(dge))
# nrDEG=as.data.frame(nrDEG)
# head(nrDEG)#查看获得的差异基因情况，nrDEG即为与P4HA1相关的差异基因
# #将结果导出来
# out <- nrDEG[nrDEG$PValue<0.05,]
# 
# geneList2 <- out$logFC
# names(geneList2) <- rownames(out)#提取P4HA1相关基因及其logFC
# geneList2 <-  sort(geneList2, decreasing = TRUE)#按logFC进行排序
# #读取下载的背景基因集合
# gsea_c7_immune_results <- GSEA(geneList2,TERM2GENE = gene_hallmark_sets) #GSEA分析
# gseaplot2(gsea_c7_immune_results, geneSetID = 1:8,
#           pvalue_table = TRUE,
#           ES_geom = "line")#"dot"将线转换为
# 
# result=result=KEGG@result
# head(result)#result即为GSEA富集的结果


# TCGA免疫浸润分析 --------------------------------------------------------------


# 设置并行工作环境，使用所有可用的CPU核心
num_cores <- detectCores() - 1  # 选择可用的核心数，避免占用所有资源
cl <- makeCluster(num_cores)    # 创建一个并行集群
registerDoParallel(cl)          # 注册并行计算环境

# 在每个节点上加载randomForest包
clusterEvalQ(cl, {
  library(IOBR)
  
})
# 将数据传递到每个节点
clusterExport(cl, c("tcga_tpm_data"))

# 运行CIBERSORT时，设置perm参数的并行化
ibersort <- deconvo_tme(eset = tcga_tpm_data,
                        method = "cibersort",
                        arrays = TRUE,
                        perm = 100  # 设置排列次数
                        )  

# 运行完毕后，停止并行计算集群
stopCluster(cl)

head(ibersort)#查看结果
#可视化，cell_bar_plot这个函数还能帮忙把数据转换成长数据后面自己改图用
res_cibersort <- cell_bar_plot(input = ibersort,
                     #输入需展示的数据
                     features = colnames(ibersort)[2:23],
                     #输入需要展示的特征列
                     title = 
                       "CIBERSORT Cell Fraction"
)
# 画好看点的箱线图还能做一下差异分析
library(ggpubr)
library(stringr)

# 分组
cibersort_long <- res_cibersort$data
colnames(TCGA_data4)[1] <- 'TCGA_ID'
cibersort_long <- cibersort_long %>%
  mutate(Group = TCGA_data4$Risk_Group[match(ID, TCGA_data4$TCGA_ID)])


p_cibersort <- ggplot(cibersort_long,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test",label.y = 0.4)
p_cibersort

#图的标题

# ESTIMATE
# 设置并行工作环境，使用所有可用的CPU核心
num_cores <- detectCores() - 1  # 选择可用的核心数，避免占用所有资源
cl <- makeCluster(num_cores)    # 创建一个并行集群
registerDoParallel(cl)          # 注册并行计算环境

clusterEvalQ(cl, {
  library(IOBR)
  
})
# 将数据传递到每个节点
clusterExport(cl, c("tcga_tpm_data"))

estimate <- deconvo_tme(eset = tcga_tpm_data, method = "estimate")

# 运行完毕后，停止并行计算集群
stopCluster(cl)
head(estimate)
res_estimate <- cell_bar_plot(input = estimate,
                               #输入需展示的数据
                               features = colnames(estimate)[2:5],
                               #输入需要展示的特征列
                               title = 
                                 "CIBERSORT Cell Fraction"
)
estimate_long <- res_estimate$data
colnames(TCGA_data4)[1] <- 'TCGA_ID'
estimate_long <- estimate_long %>%
  mutate(Group = TCGA_data4$Risk_Group[match(ID, TCGA_data4$TCGA_ID)])


p_estimate <- ggplot(estimate_long,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test",label.y = 0.4)
p_estimate
# 画好看点的箱线图还能做一下差异

# xCELL
# 设置并行工作环境，使用所有可用的CPU核心
num_cores <- detectCores() - 1  # 选择可用的核心数，避免占用所有资源
cl <- makeCluster(num_cores)    # 创建一个并行集群
registerDoParallel(cl)          # 注册并行计算环境

clusterEvalQ(cl, {
  library(IOBR)
  
})
# 将数据传递到每个节点
clusterExport(cl, c("tcga_tpm_data"))

xcell <- deconvo_tme(eset = tcga_tpm_data, method = "xcell", arrays = TRUE)
# 运行完毕后，停止并行计算集群
stopCluster(cl)
head(xcell)
res_xcell <- cell_bar_plot(input = xcell,
                              #输入需展示的数据
                              features = colnames(xcell)[2:6],
                              #输入需要展示的特征列
                              title = 
                                "CIBERSORT Cell Fraction"
)
xcell_long <- res_xcell$data
colnames(TCGA_data4)[1] <- 'TCGA_ID'
xcell_long <- xcell_long %>%
  mutate(Group = TCGA_data4$Risk_Group[match(ID, TCGA_data4$TCGA_ID)])


p_xcell <- ggplot(xcell_long,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test",label.y = 0.4)
p_xcell
# 
# #查看结果
# 
# cibersort_result <- deconvo_cibersort(
#   eset = tcga_tpm_data,  # 你的TPM数据
#   project = "TCGA_Project",  # 项目名称
#   arrays = FALSE,  # 数据是RNA-seq，不是微阵列
#   perm = 1000,  # 排列次数
#   absolute = FALSE  # 不使用绝对模式
# )
# 
# # 查看结果
# head(cibersort_result)
# 
# # 合并免疫浸润结果和Risk Group
# merged_data_immune <- merge(cibersort_result, TCGA_data4, by.x = "ID", by.y = "TCGA_data2$TCGA_ID")
# 
# # 以T细胞为例
# ggplot(merged_data_immune, aes(x = Risk_Group, y = T_cells_CD8, fill = Risk_Group)) +
#   geom_boxplot() +
#   stat_compare_means(method = "t.test") +  # 添加t检验
#   theme_minimal() +
#   labs(title = "T细胞浸润差异", x = "Risk Group", y = "T细胞浸润分数")

# TCGA突变分析 ----------------------------------------------------------------

# 突变数据
TCGAmutations::tcga_available()
#分别获得肺腺癌与肺鳞癌的maf数据
luad.maf <- TCGAmutations::tcga_load(study = "LUAD")
lusc.maf <- TCGAmutations::tcga_load(study = "LUSC")

#合并两个maf文件
all.maf<- maftools::merge_mafs(list(luad.maf,lusc.maf))
clindata <- all.maf@clinical.data
length(unique(all.maf@data$Tumor_Sample_Barcode))

# 合并高低风险组数据
colnames(clindata)[1] <- "TCGA_ID"
intersect(clindata$TCGA_ID, TCGA_data4$TCGA_ID)
clindata <- merge(
  x = TCGA_data4,      # 保留TCGA_data4的所有样本（左连接）
  y = clindata,        # 要合并的临床数据
  by = "TCGA_ID",      # 以TCGA_ID列为键
  all.x = TRUE         # 确保TCGA_data4的所有行保留（即使clindata无匹配）
)
table(clindata$Risk_Group)
# 将更新后的临床数据赋值回 MAF 对象
# 确保 clindata 是 data.table
setDT(clindata)  # 将 data.frame 转为 data.table

# 现在可以安全赋值
all.maf@clinical.data <- clindata

# 提取了每个样本的名称(前者是1：12位的/后者是全称)
samp = data.frame(ID = str_sub(unique(all.maf@data$Tumor_Sample_Barcode),1,12),
                  long = unique(all.maf@data$Tumor_Sample_Barcode))
samp <- samp[samp$ID%in%clindata$TCGA_ID,]
# all.maf@data <- all.maf@data[all.maf@data$Tumor_Sample_Barcode_min%in%clindata$TCGA_ID,]
# length(unique(all.maf@data$Tumor_Sample_Barcode))
# subsetMaf()是maftools包中的一个函数，用于从原始MAF 数据中提取特定的样本或突变类型的子集。
all.maf = subsetMaf(all.maf, tsb = samp$long)
length(unique(all.maf@data$Tumor_Sample_Barcode))
# # 从 MAF 中提取样本 ID 并去除后缀（保留前12字符）
# maf_samples <- substr(all.maf@data$Tumor_Sample_Barcode_min, 1, 12)
# # 获取有 Risk_Group 标注的样本名（排除NA）
# valid_samples <- clindata[, 1]
# 
# # 找出 MAF 中与 valid_samples 匹配的样本（完整 Barcode）
# matched_barcodes <- all.maf@data$Tumor_Sample_Barcode_min[
#   maf_samples %in% valid_samples
# ] %>% unique()
# 
# # 子集化 MAF 对象
# filtered.maf <- subsetMaf(
#   maf = all.maf,
#   tsb = valid_samples,  # tsb = Tumor Sample Barcode
#   tsbCol = "Tumor_Sample_Barcode_min", # 关键参数：指定匹配列
#   mafObj = TRUE         # 返回 MAF 对象
# )


# 再使用subsetMaf
all.maf.high <- subsetMaf(
  maf = all.maf,
  mafObj = TRUE,
  clinQuery = "Risk_Group == 'high'"  # 现在可以正确解析
)

# 低风险组的 MAF
all.maf.low <- subsetMaf(
  maf = all.maf,
  clinQuery = "Risk_Group == 'low'",
  mafObj = TRUE
)

#突变全景分析
plotmafSummary(maf = all.maf, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE, 
               titvRaw = FALSE)

# 获取 all.maf 的前20个基因（按突变频率排序）
gene_summary <- maftools::getGeneSummary(all.maf)
top20_genes <- gene_summary[1:25, Hugo_Symbol]
# top20_genes <- as.factor(top20_genes)
#前20突变基因展示
vc_cols <- c(
  'Missense_Mutation' = "#1F78B4",
  'Multi_Hit' = "#33A02C",  # 确保顺序固定
  'Splice_Site' = "#FDBF6F",
  'Frame_Shift_Ins' = "#FB9A99",
  'Frame_Shift_Del' = "#A6CEE3",
  'Nonsense_Mutation' = "#B2DF8A"
)
oncoplot(maf = all.maf, colors = vc_cols, genes = top20_genes, keepGeneOrder = T         # 指定基因列表
         # geneOrder = top20_genes # 固定顺序)
         )     


oncoplot(maf = all.maf.high, colors = vc_cols, genes = top20_genes, keepGeneOrder = T)
oncoplot(maf = all.maf.low, colors = vc_cols, genes = top20_genes, keepGeneOrder = T)

#对比高低风险突变图谱
pt.vs.rt <- mafCompare(m1 = all.maf.high, 
                       m2 = all.maf.low, 
                       m1Name = 'high', 
                       m2Name = 'low', minMut = 0)
print(pt.vs.rt)

vs.re.sig<- pt.vs.rt$results%>%
  filter(pval<0.05)

vs.re.sig$Hugo_Symbol[vs.re.sig$Hugo_Symbol%in%top20_genes]





# 定义突变参数
cla1 = c('Frame_Shift_Del','Missense_Mutation', 
         'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
         'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
# 上面的信息藏在variant.classification.summary中
cla2 = colnames(all.maf@variant.classification.summary)[-c(1,ncol(all.maf@variant.classification.summary))]
# 使用intersect函数找到cla1和cla2的交集，即共同的突变分类；然后使用setdiff函数找到cla2中不在cla1中的突变分类，将它们合并以形成最终的突变分类列表cla
cla = c(intersect(cla1,cla2),setdiff(cla2,cla1))
col = RColorBrewer::brewer.pal(n = length(cla), name = 'Paired')
fabcol = RColorBrewer::brewer.pal(n = 9,name = 'Spectral')
riskgroupcolor = fabcol[8:9];names(riskgroupcolor) = na.omit(unique(all.maf@clinical.data$Risk_Group))
oncoplot(maf = all.maf, top = 20,
         clinicalFeatures = "Risk_Group", 
         sortByAnnotation = TRUE, 
         annotationColor = list(Risk_Group = riskgroupcolor))


#碱基改变情况
laml.titv = titv(maf = all.maf, plot = T, useSyn = TRUE)


#共突变或者共斥突变
somaticInteractions(maf = all.maf, top = 25, pvalue = c(0.05, 0.1))


#超突变分析
rainfallPlot(maf = all.maf, detectChangePoints = TRUE, pointSize = 0.4)


# 驱动突变
all.sig = oncodrive(maf = all.maf,
                    AACol = 'HGVSp_Short',
                    minMut = 5,
                    pvalMethod = 'zscore')

# plotOncodrive(res =all.sig, 
#               fdrCutOff = 0.01, 
#               useFraction = TRUE, 
#               labelSize = 0.2)



#绘图展示不同年龄分组的差异突变情况
coBarplot(m1 = all.maf.high, 
          m2 =all.maf.low, 
          m1Name = 'high', 
          m2Name = 'low' ,
          genes = vs.re.sig$Hugo_Symbol)  # 关键参数：仅显示显著基因)

dgi = drugInteractions(maf = all.maf, fontSize = 0.75,genes = vs.re.sig$Hugo_Symbol)
dnmt3a.dgi = drugInteractions(genes = vs.re.sig$Hugo_Symbol, drugs = TRUE)


#计算肿瘤异质性指数math并对比青年与老年
high.math<- math.score(all.maf.high)%>%
  mutate(team="high")
low.math<- math.score(all.maf.low)%>%
  mutate(team="low")

all.math <- rbind(high.math,low.math)

ggstatsplot::ggbetweenstats(data = all.math,x = "team",y = "MATH")


# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
# library('NMF')
# 
# all.tnm = trinucleotideMatrix(maf = all.maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
# 
# 
# all.sign = estimateSignatures(mat = all.tnm, nTry = 9)
# plotCophenetic(res = all.sign)
# laml.sig = extractSignatures(mat = all.tnm, n = 5)

# # 绘制ROC曲线 -----------------------------------------------------------------
# 
# # 算一下原始数据中3，6，9, 12, 18, 24, 36个月的进展率，以决定使用哪个time point
# # 设置时间点
# time_points <- c((1*30+2*31), (3*30+3*31), (4*30+5*31), (6*30+6*31), (9*30+9*31), (12*30+12*31),  (18*30+18*31))
# # 创建一个空的结果数据框
# pfs_progression_rates <- data.frame(Time_Point = time_points, Progression_Rate = NA)
# # 循环计算各时间点的进展率
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算进展率
#   data_filled_factor_date <- data_filled_factor
#   data_filled_factor_date$Progression_Status <- ifelse(data_filled_factor_date$PFS <= current_time & data_filled_factor_date$Progression == 1, 1, 0)
#   progression_rate <- mean(data_filled_factor_date$Progression_Status, na.rm = TRUE)  # 计算进展率
#   
#   pfs_progression_rates$Progression_Rate[i] <- progression_rate
# }
# # 查看结果
# print(pfs_progression_rates)
# 
# # 计算随访数
# # 创建一个空的数据框来存储结果
# followup_counts <- data.frame(Time_Point = time_points, Patients_at_Risk = NA)
# 
# # 计算每个时间点的在随访患者数
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算在随访中的患者数
#   patients_at_risk <- sum(data_filled_factor$PFS >= current_time, na.rm = TRUE)
#   
#   followup_counts$Patients_at_Risk[i] <- patients_at_risk
# }
# 
# # 查看每个时间点的在随访患者数
# print(followup_counts)
# 
# 
# # 所以3年的进展率是最接近50%的
# 
# 
# # 汇报给robin之后
# # 把1年、2年、3年的ROC都画出来
# # 假设predictions是模型预测的风险评分
# # 使用累积风险分数
# 
# 
# set.seed(which.max(all_Cindex))
# fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor, 
#              ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
#              splitrule = 'logrank', 
#              importance = T, 
#              proximity = T, 
#              forest = T, 
#              seed = which.max(all_Cindex))
# risk_scores <- predict(fit, data_filled_factor, type = "risk")$predicted
# 
# # 提取生存时间和状态
# time <- data_filled_factor$PFS
# status <- data_filled_factor$Progression  # 0或1表示无进展或有进展
# 
# library(timeROC)
# riskRoc <- timeROC(T = time,delta = status,
#                    marker = risk_scores,cause = 1,
#                    weighting="marginal",
#                    times = c(365,730,1095))
# 
# multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
#   library(ggsci)
#   color <- pal_lancet()(length(time))
#   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
#        col=color[1],
#        xlab=xlab,
#        ylab=ylab,main=title)
#   #如果直接plot roc对象，无法修改标题和坐标轴标签
#   for(i in 2:length(time)){
#     plot(ROC,time=time[i],add=T,col=color[i])
#   }
#   legend("bottomright",
#          legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
#          col = color,lwd = 1,
#          bty = "n",cex = cex,text.col = color
#   )
# }
# multiTimeplot(riskRoc,time = c(365,730,1095),
#               title="Time dependent ROC curve",
#               xlab="False positive rate",
#               ylab="True positive rate",
#               cex=0.7)
# riskRoc
# 
# # 上面是训练集和验证集合起来的，接下来要把训练集和验证集拆开AUC展示
# set.seed(464)
# trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
# trainData_factor <- data_filled_factor[trainIndex, ]
# testData_factor <- data_filled_factor[-trainIndex, ]
# 
# risk_scores_train <- predict(fit, trainData_factor, type = "risk")$predicted
# 
# # 提取生存时间和状态
# time_train <- trainData_factor$PFS
# status_train <- trainData_factor$Progression  # 0或1表示无进展或有进展
# 
# riskRoc_train <- timeROC(T = time_train,delta = status_train,
#                          marker = risk_scores_train,cause = 1,
#                          weighting="marginal",
#                          times = c(365,730,1095))
# 
# multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
#   library(ggsci)
#   color <- pal_lancet()(length(time))
#   # 设置绘图参数为正方形
#   par(pty = "s")
#   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
#        col=color[1],
#        xlab=xlab,
#        ylab=ylab,main=title)
#   #如果直接plot roc对象，无法修改标题和坐标轴标签
#   for(i in 2:length(time)){
#     plot(ROC,time=time[i],add=T,col=color[i])
#   }
#   legend(x = 0.33, y = 0.4,
#          legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
#          col = color,lwd = 1, y.intersp = 0.5,
#          bty = "n",cex = cex,text.col = color
#   )
# }
# multiTimeplot(riskRoc_train,time = c(365,730,1095),
#               title="Time Dependent ROC Curve for the Claude Matrix Training Set",
#               xlab="False Positive Rate",
#               ylab="True Positive Rate",
#               cex=0.7)
# riskRoc_train
# 
# 
# risk_scores_test <- predict(fit, testData_factor, type = "risk")$predicted
# 
# # 提取生存时间和状态
# time_test <- testData_factor$PFS
# status_test <- testData_factor$Progression  # 0或1表示无进展或有进展
# 
# riskRoc_test <- timeROC(T = time_test,delta = status_test,
#                         marker = risk_scores_test,cause = 1,
#                         weighting="marginal",
#                         times = c(365,730,1095))
# 
# multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
#   library(ggsci)
#   color <- pal_lancet()(length(time))
#   # 设置绘图参数为正方形
#   par(pty = "s")
#   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
#        col=color[1],
#        xlab=xlab,
#        ylab=ylab,main=title)
#   #如果直接plot roc对象，无法修改标题和坐标轴标签
#   for(i in 2:length(time)){
#     plot(ROC,time=time[i],add=T,col=color[i])
#   }
#   legend(x = 0.33, y = 0.4,
#          legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
#          col = color,lwd = 1, y.intersp = 0.5,
#          bty = "n",cex = cex,text.col = color
#   )
# }
# multiTimeplot(riskRoc_test,time = c(365,730,1095),
#               title="Time Dependent ROC Curve for the Claude Matrix Validation Set",
#               xlab="False Positive Rate",
#               ylab="True Positive Rate",
#               cex=0.7)
# riskRoc_test
# 
# # 
# # ##########导入GPT的数据跑ROC
# # ctfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/GPT4整理结果.xlsx"
# # purectdata <- read_excel(ctfile_path)
# # ctidfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Copy of 3_ID.xlsx"
# # ctiddata <- read_excel(ctidfile_path)
# # ctdata <- merge(purectdata, ctiddata, by = "文件夹名", all.x = TRUE)
# # ctdata <- ctdata[,-c(1:2,56:60)]
# # # 把最后一列放到第一列
# # ctdata <- ctdata[, c(ncol(ctdata), 1:(ncol(ctdata)-1))]
# # 
# # pfsfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/汇总.xlsx"
# # pfsfile <- read_excel(pfsfile_path)
# # colnames(pfsfile) <- pfsfile[1,]
# # pfsfile <- pfsfile[-1,]
# # data <- merge(ctdata, pfsfile, by = "Hospitalization_Number", all.x = TRUE)
# # # 将第 57, 58, 60, 61 列的数字转换为日期
# # data[, c(56, 57, 59, 60)] <- lapply(data[, c(56, 57, 59, 60)], function(x) as.Date(as.numeric(x), origin = "1899-12-30"))
# # # 在第58和59列之间插入一列，为58列减去59列的天数差
# # data <- data %>%
# #   mutate(Difference_57_56 = as.numeric(data[[57]] - data[[56]])) %>%
# #   relocate(Difference_57_56, .after = 57)
# # 
# # # 同理
# # data <- data %>%
# #   mutate(Difference_61_60 = as.numeric(data[[61]] - data[[60]])) %>%
# #   relocate(Difference_61_60, .after = 61)
# # 
# # 
# # sum(is.na(data))
# # 
# # # KNN填补缺失值
# # data_filled_ <- kNN(data, k = 5)
# # # 不知道有些无限接近于0的负值表现出来就是0还是生存分析不能接受0值，反正就有0的数据容易报错
# # data_filled_ <- data_filled_[-which(data_filled_[,58]<=0),]
# # sum(is.na(data_filled_))
# # 
# # 
# # # 因为它是是否存活和是否进展，所以要反过来，存活是0，不存活才是1
# # data_filled_[data_filled_[, 59] %in% "否", 59] <- 0
# # data_filled_[data_filled_[, 59] %in% "是", 59] <- 1
# # data_filled_[data_filled_[, 63] %in% "是", 63] <- 0
# # data_filled_[data_filled_[, 63] %in% "否", 63] <- 1
# # 
# # 
# # colnames(data_filled_)[58] <- "PFS"
# # colnames(data_filled_)[59] <- "Progression"
# # colnames(data_filled_)[62] <- "OS.time"
# # colnames(data_filled_)[63] <- "OS"
# # 
# # 
# # 
# # # 将所有变量转为数值型
# # set.seed(123)
# # data_filled <- data_filled_[,c(58, 59, 2:53)]
# # 
# # 
# # # cols <- c(3:8, 10:34)
# # data_filled_factor <- data_filled
# # data_filled_factor[, 2] <- as.numeric(data_filled_factor[, 2])
# # str(data_filled_factor)
# # 
# # set.seed(464)
# # trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
# # trainData_factor <- data_filled_factor[trainIndex, ]
# # testData_factor <- data_filled_factor[-trainIndex, ]
# # 
# # risk_scores_train <- predict(fit, trainData_factor, type = "risk")$predicted
# # 
# # # 提取生存时间和状态
# # time_train <- trainData_factor$PFS
# # status_train <- trainData_factor$Progression  # 0或1表示无进展或有进展
# # 
# # riskRoc_train <- timeROC(T = time_train,delta = status_train,
# #                          marker = risk_scores_train,cause = 1,
# #                          weighting="marginal",
# #                          times = c(365,730,1095))
# # 
# # multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
# #   library(ggsci)
# #   color <- pal_lancet()(length(time))
# #   # 设置绘图参数为正方形
# #   par(pty = "s")
# #   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
# #        col=color[1],
# #        xlab=xlab,
# #        ylab=ylab,main=title)
# #   #如果直接plot roc对象，无法修改标题和坐标轴标签
# #   for(i in 2:length(time)){
# #     plot(ROC,time=time[i],add=T,col=color[i])
# #   }
# #   legend(x = 0.33, y = 0.4,
# #          legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
# #          col = color,lwd = 1, y.intersp = 0.5,
# #          bty = "n",cex = cex,text.col = color
# #   )
# # }
# # multiTimeplot(riskRoc_train,time = c(365,730,1095),
# #               title="Time Dependent ROC Curve for the GPT Matrix Training Set(n=270)",
# #               xlab="False Positive Rate",
# #               ylab="True Positive Rate",
# #               cex=0.7)
# # riskRoc_train
# # 
# # 
# # risk_scores_test <- predict(fit, testData_factor, type = "risk")$predicted
# # 
# # # 提取生存时间和状态
# # time_test <- testData_factor$PFS
# # status_test <- testData_factor$Progression  # 0或1表示无进展或有进展
# # 
# # riskRoc_test <- timeROC(T = time_test,delta = status_test,
# #                         marker = risk_scores_test,cause = 1,
# #                         weighting="marginal",
# #                         times = c(365,730,1095))
# # 
# # multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
# #   library(ggsci)
# #   color <- pal_lancet()(length(time))
# #   # 设置绘图参数为正方形
# #   par(pty = "s")
# #   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
# #        col=color[1],
# #        xlab=xlab,
# #        ylab=ylab,main=title)
# #   #如果直接plot roc对象，无法修改标题和坐标轴标签
# #   for(i in 2:length(time)){
# #     plot(ROC,time=time[i],add=T,col=color[i])
# #   }
# #   legend(x = 0.33, y = 0.4,
# #          legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
# #          col = color,lwd = 1, y.intersp = 0.5,
# #          bty = "n",cex = cex,text.col = color
# #   )
# # }
# # multiTimeplot(riskRoc_test,time = c(365,730,1095),
# #               title="Time Dependent ROC Curve for the GPT Matrix Validation Set(n=67)",
# #               xlab="False Positive Rate",
# #               ylab="True Positive Rate",
# #               cex=0.7)
# # riskRoc_test
# # 
# # #########接下来轮到Gemini的
# # 
# # ctfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Gemini整理结果.xlsx"
# # purectdata <- read_excel(ctfile_path)
# # ctidfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/Copy of 3_ID.xlsx"
# # ctiddata <- read_excel(ctidfile_path)
# # ctdata <- merge(purectdata, ctiddata, by = "文件夹名", all.x = TRUE)
# # ctdata <- ctdata[,-c(1:2,56:60)]
# # # 把最后一列放到第一列
# # ctdata <- ctdata[, c(ncol(ctdata), 1:(ncol(ctdata)-1))]
# # 
# # pfsfile_path <- "/Users/yetaojun/Documents/学习/硕士/硕士课题/GPT+CT+免疫标志物课题/小测试/本地数据/免疫随访汇总/GenAI数据/汇总.xlsx"
# # pfsfile <- read_excel(pfsfile_path)
# # colnames(pfsfile) <- pfsfile[1,]
# # pfsfile <- pfsfile[-1,]
# # data <- merge(ctdata, pfsfile, by = "Hospitalization_Number", all.x = TRUE)
# # # 将第 57, 58, 60, 61 列的数字转换为日期
# # data[, c(56, 57, 59, 60)] <- lapply(data[, c(56, 57, 59, 60)], function(x) as.Date(as.numeric(x), origin = "1899-12-30"))
# # # 在第58和59列之间插入一列，为58列减去59列的天数差
# # data <- data %>%
# #   mutate(Difference_57_56 = as.numeric(data[[57]] - data[[56]])) %>%
# #   relocate(Difference_57_56, .after = 57)
# # 
# # # 同理
# # data <- data %>%
# #   mutate(Difference_61_60 = as.numeric(data[[61]] - data[[60]])) %>%
# #   relocate(Difference_61_60, .after = 61)
# # 
# # 
# # sum(is.na(data))
# # 
# # # KNN填补缺失值
# # data_filled_ <- kNN(data, k = 5)
# # # 不知道有些无限接近于0的负值表现出来就是0还是生存分析不能接受0值，反正就有0的数据容易报错
# # data_filled_ <- data_filled_[-which(data_filled_[,58]<=0),]
# # sum(is.na(data_filled_))
# # 
# # 
# # # 因为它是是否存活和是否进展，所以要反过来，存活是0，不存活才是1
# # data_filled_[data_filled_[, 59] %in% "否", 59] <- 0
# # data_filled_[data_filled_[, 59] %in% "是", 59] <- 1
# # data_filled_[data_filled_[, 63] %in% "是", 63] <- 0
# # data_filled_[data_filled_[, 63] %in% "否", 63] <- 1
# # 
# # 
# # colnames(data_filled_)[58] <- "PFS"
# # colnames(data_filled_)[59] <- "Progression"
# # colnames(data_filled_)[62] <- "OS.time"
# # colnames(data_filled_)[63] <- "OS"
# # 
# # 
# # 
# # # 将所有变量转为数值型
# # set.seed(123)
# # data_filled <- data_filled_[,c(58, 59, 2:53)]
# # 
# # 
# # # cols <- c(3:8, 10:34)
# # data_filled_factor <- data_filled
# # data_filled_factor[, c(2,5:9)] <- lapply(data_filled_factor[, c(2,5:9)], function(x) as.numeric(x))
# # str(data_filled_factor)
# # 
# # set.seed(464)
# # trainIndex <- createDataPartition(data_filled_factor$Progression, p = 0.8, list = FALSE)
# # trainData_factor <- data_filled_factor[trainIndex, ]
# # testData_factor <- data_filled_factor[-trainIndex, ]
# # 
# # risk_scores_train <- predict(fit, trainData_factor, type = "risk")$predicted
# # 
# # # 提取生存时间和状态
# # time_train <- trainData_factor$PFS
# # status_train <- trainData_factor$Progression  # 0或1表示无进展或有进展
# # 
# # riskRoc_train <- timeROC(T = time_train,delta = status_train,
# #                          marker = risk_scores_train,cause = 1,
# #                          weighting="marginal",
# #                          times = c(365,730,1095))
# # 
# # multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
# #   library(ggsci)
# #   color <- pal_lancet()(length(time))
# #   # 设置绘图参数为正方形
# #   par(pty = "s")
# #   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
# #        col=color[1],
# #        xlab=xlab,
# #        ylab=ylab,main=title)
# #   #如果直接plot roc对象，无法修改标题和坐标轴标签
# #   for(i in 2:length(time)){
# #     plot(ROC,time=time[i],add=T,col=color[i])
# #   }
# #   legend(x = 0.33, y = 0.4,
# #          legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
# #          col = color,lwd = 1, y.intersp = 0.5,
# #          bty = "n",cex = cex,text.col = color
# #   )
# # }
# # multiTimeplot(riskRoc_train,time = c(365,730,1095),
# #               title="Time Dependent ROC Curve for the Gemini Matrix Training Set(n=270)",
# #               xlab="False Positive Rate",
# #               ylab="True Positive Rate",
# #               cex=0.7)
# # riskRoc_train
# # 
# # 
# # risk_scores_test <- predict(fit, testData_factor, type = "risk")$predicted
# # 
# # # 提取生存时间和状态
# # time_test <- testData_factor$PFS
# # status_test <- testData_factor$Progression  # 0或1表示无进展或有进展
# # 
# # riskRoc_test <- timeROC(T = time_test,delta = status_test,
# #                         marker = risk_scores_test,cause = 1,
# #                         weighting="marginal",
# #                         times = c(365,730,1095))
# # 
# # multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
# #   library(ggsci)
# #   color <- pal_lancet()(length(time))
# #   # 设置绘图参数为正方形
# #   par(pty = "s")
# #   plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1),
# #        col=color[1],
# #        xlab=xlab,
# #        ylab=ylab,main=title)
# #   #如果直接plot roc对象，无法修改标题和坐标轴标签
# #   for(i in 2:length(time)){
# #     plot(ROC,time=time[i],add=T,col=color[i])
# #   }
# #   legend(x = 0.33, y = 0.4,
# #          legend =paste("AUC at",c("1","2","3"),"year:",round(ROC$AUC,digits = 3)),
# #          col = color,lwd = 1, y.intersp = 0.5,
# #          bty = "n",cex = cex,text.col = color
# #   )
# # }
# # multiTimeplot(riskRoc_test,time = c(365,730,1095),
# #               title="Time Dependent ROC Curve for the Gemini Matrix Validation Set(n=67)",
# #               xlab="False Positive Rate",
# #               ylab="True Positive Rate",
# #               cex=0.7)
# # riskRoc_test
# # nomogram绘制 --------------------------------------------------------------
# 
# 
# 
# # 算一下原始数据中3，6，9, 12, 18, 24, 36个月的进展率，以决定使用哪个time point
# # 设置时间点
# time_points <- c((1*30+2*31), (3*30+3*31), (4*30+5*31), (6*30+6*31), (9*30+9*31), (12*30+12*31),  (18*30+18*31))
# # 创建一个空的结果数据框
# pfs_progression_rates <- data.frame(Time_Point = time_points, Progression_Rate = NA)
# # 循环计算各时间点的进展率
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算进展率
#   data_filled_factor_date <- data_filled_factor
#   data_filled_factor_date$Progression_Status <- ifelse(data_filled_factor_date$PFS <= current_time & data_filled_factor_date$Progression == 1, 1, 0)
#   progression_rate <- mean(data_filled_factor_date$Progression_Status, na.rm = TRUE)  # 计算进展率
#   
#   pfs_progression_rates$Progression_Rate[i] <- progression_rate
# }
# # 查看结果
# print(pfs_progression_rates)
# 
# # 计算随访数
# # 创建一个空的数据框来存储结果
# followup_counts <- data.frame(Time_Point = time_points, Patients_at_Risk = NA)
# 
# # 计算每个时间点的在随访患者数
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算在随访中的患者数
#   patients_at_risk <- sum(data_filled_factor$PFS >= current_time, na.rm = TRUE)
#   
#   followup_counts$Patients_at_Risk[i] <- patients_at_risk
# }
# 
# # 查看每个时间点的在随访患者数
# print(followup_counts)
# 
# 
# # 所以3年的进展率是最接近50%的
# 
# 
# # 汇报给robin之后
# # 把1年、2年、3年的ROC都画出来
# # 假设predictions是模型预测的风险评分
# # 使用累积风险分数
# data_filled_factor_pre <- data_filled_factor
# set.seed(which.max(all_Cindex))
# # 使用RSF模型生成预测结果，type = "risk" 提取风险评分
# risk_scores <- predict(fit, data_filled_factor_pre, type = "risk")$predicted
# 
# # 提取生存时间和状态
# time <- data_filled_factor_pre$PFS
# status <- data_filled_factor_pre$Progression  # 0或1表示无进展或有进展
# 
# # 计算365天、730天和1095天的时间依赖性ROC
# roc_1_year <- timeROC(T = time, delta = status, marker = risk_scores, 
#                       cause = 1, times = 365, iid = TRUE)
# roc_2_years <- timeROC(T = time, delta = status, marker = risk_scores, 
#                        cause = 1, times = 730, iid = TRUE)
# roc_3_years <- timeROC(T = time, delta = status, marker = risk_scores, 
#                        cause = 1, times = 1095, iid = TRUE)
# 
# # 绘制ROC曲线
# plot(roc_1_year, time = 365, col = "blue", lwd = 2, main = "Time-dependent ROC Curves")
# plot(roc_2_years, time = 730, add = TRUE, col = "green", lwd = 2)
# plot(roc_3_years, time = 1095, add = TRUE, col = "red", lwd = 2)
# 
# # 获取AUC值
# auc_1_year <- round(roc_1_year$AUC[2], 3)
# auc_2_years <- round(roc_2_years$AUC[2], 3)
# auc_3_years <- round(roc_3_years$AUC[2], 3)
# # 添加AUC值到图中
# text(0.4, 0.4, paste("1-year AUC =", auc_1_year), col = "blue", cex = 1)
# text(0.4, 0.35, paste("2-year AUC =", auc_2_years), col = "green", cex = 1)
# text(0.4, 0.3, paste("3-year AUC =", auc_3_years), col = "red", cex = 1)
# 
# # 添加图例
# legend("bottomright", legend = c("1-year ROC", "2-year ROC", "3-year ROC"),
#        col = c("blue", "green", "red"), lwd = 2)
# 
# # 提取AUC值
# cat("1-year AUC:", roc_1_year$AUC[2], "\n")
# cat("2-year AUC:", roc_2_years$AUC[2], "\n")
# cat("3-year AUC:", roc_3_years$AUC[2], "\n")
# 
# 
# # # 设定需要观察的时间点
# # time_points <- c(36)  # 单位为月
# # 
# # # 找到每个时间点最接近的索引
# # model_times <- fit$time  # 模型中的实际时间点
# # closest_time_points <- which.min(abs(model_times - time_points))
# # 
# # # 计算每个时间点的累积风险
# # risk_at_times <- scores[, closest_time_points]
# # 
# # # 计算进展率（1 - 生存概率）每个时间点的进展率
# # progression_rate <- 1 - exp(-risk_at_times)  # 累积风险转换为进展率
# # 
# # # 查找进展率为50%的时间点
# # progression_rate_50 <- which.max(progression_rate >= 0.5)
# # cat("进展率达到50%的时间点：", time_points[progression_rate_50], "个月\n")
# # 
# # # 打印每个时间点的进展率
# # print(progression_rate)
# 
# 
# # 根据depth threshold 这个值，来确定筛选出的变量
# # 当小于阈值时就得到了在VIMP法和最小深度法结合的情况下筛选出的变量
# # depth threshold    : 8.2095
# 
# # 接下来要取importance变量画nomogram了
# # 提取变量重要性
# # 假设varselect已经存储在importance_values$varselect中
# # 提取varselect信息
# importance_values <- var.select(fit)
# 
# varselect_df <- importance_values$varselect
# 
# # 提取变量名和重要性值
# importance_df <- data.frame(
#   Variable = rownames(varselect_df),
#   Importance = varselect_df$vimp
# )
# 
# # 按照vimp值降序排序
# importance_df <- importance_df[order(-importance_df$Importance), ]
# 
# # 绘制变量重要性图
# ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "Variable Importance from RSF Model",
#        x = "Variable",
#        y = "Variable Importance (vimp)")
# 
# # 1. 选择vimp值大于0的变量
# selected_vars1 <- importance_df$Variable[importance_df$Importance > 0]
# 
# # 2. 试一下VIMP法结合最小深度法
# #最小深度法查看变量重要性
# library(ggRandomForests)
# gg_dta_depth <- gg_minimal_depth(fit)
# plot(gg_dta_depth)
# #两种方法的结合  VIMP+min_depth
# gg_dta_vimp <- gg_minimal_vimp(fit)
# plot(gg_dta_vimp)
# 
# depth_threshold <- 8.2095
# # 获取最小深度小于阈值的变量
# selected_vars2 <- gg_dta_depth$topvars
# 
# selected_vars <- intersect(selected_vars1, selected_vars2)
# 
# # 打印出选择的变量
# print(selected_vars)
# cox_data <- data_filled_factor[, c("PFS", "Progression", selected_vars)]
# # cox_data$SurvObj <- with(cox_data, Surv(PFS, Progression))
# # 准备一个空的数据框，用于存储Cox回归分析的结果
# cox_results <- data.frame(Variable = character(), P_value = numeric(), stringsAsFactors = FALSE)
# # 循环遍历每个变量，进行单变量Cox回归分析
# # 循环遍历每个变量，进行单变量Cox回归分析
# for (var in selected_vars) {
#   # 动态构造公式
#   formula <- as.formula(paste("Surv(PFS, Progression) ~", var))
#   
#   # 执行Cox回归分析
#   cox_model <- coxph(formula, data = cox_data)
#   
#   # 提取p值
#   p_value <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]
#   
#   # 将结果存入数据框
#   cox_results <- rbind(cox_results, data.frame(Variable = var, P_value = p_value))
# }
# # 筛选出p值显著的变量，例如 p < 0.05
# significant_vars <- cox_results$Variable[cox_results$P_value < 0.05]
# # 打印结果
# print(significant_vars)
# # 构建多变量Cox模型
# final_formula <- as.formula(paste("Surv(PFS, Progression) ~", paste(significant_vars, collapse = " + ")))
# # head(cox_data)
# # str(cox_data)
# library(rms)
# cox_rms_model <- cph(final_formula, data = cox_data, x = TRUE, y = TRUE, surv = TRUE)
# dd <- datadist(cox_data)
# options(datadist = "dd")
# med <- Quantile(cox_rms_model) #计算中位生存时间
# surv <- Survival(cox_rms_model) #构建生存概率函数
# 
# # 设置全局小数点精度为两位
# options(digits = 2)
# 
# nom <- nomogram(cox_rms_model, 
#                 fun = list(function(x) surv(365, x),
#                            function(x) surv(365*2, x),
#                            function(x) surv(365*3, x)),
#                 funlabel = c("1-Year Progression Rate", "2-Year Progression Rate", "3-Year Progression Rate"),
#                 lp = T)
# plot(nom)
# 
# 
# # RSF_for循环 ---------------------------------------------------------------
# 
# # 设置并行核心数量
# num_cores <- detectCores() - 1
# cl <- makeCluster(num_cores)
# clusterSetRNGStream(cl, 1)
# 
# # 在每个节点上加载randomForest包
# clusterEvalQ(cl, {
#   library(randomForest)
#   library(pROC) 
#   library(randomForestSRC)
#   library(magrittr)
#   library(tibble)
#   library(survival)
# })
# # 将数据传递到每个节点
# clusterExport(cl, c("data_filled_factor","rf_nodesize"))
# # 运行并行计算
# results <- parLapply(cl, 1:1000, function(seed) {
#   set.seed(seed)
#   fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor, 
#                ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
#                splitrule = 'logrank', 
#                importance = T, 
#                proximity = T, 
#                forest = T, 
#                seed = seed)
#   
#   set.seed(seed)
#   rs <- predict(fit, newdata = data_filled_factor)$predicted
#   set.seed(seed)
#   model_result <- data.frame(Cindex = as.numeric(summary(coxph(Surv(PFS, Progression) ~ rs, data_filled_factor))$concordance[1])) %>%
#     rownames_to_column('ID')
#   return(model_result)
#   
# })
# 
# stopCluster(cl)
# # 计算平均AUC并筛选出平均AUC大于0.8的模型
# all_Cindex <- sapply(results, function(x) x$Cindex)
# all_Cindex
# max(all_Cindex)
# which.max(all_Cindex)
# average_Cindex <- mean(unlist(all_Cindex))
# average_Cindex
# 
# 
# # 替模型训练函数
# set.seed(which.max(all_Cindex))  
# fit <- rfsrc(Surv(PFS,Progression)~., data = data_filled_factor, 
#              ntree = 1000, nodesize = rf_nodesize,  #该值建议多调整
#              splitrule = 'logrank', 
#              importance = T, 
#              proximity = T, 
#              forest = T, 
#              seed = which.max(all_Cindex))
# set.seed(which.max(all_Cindex))  
# 
# 
# 
# # 使用累积风险分数
# set.seed(which.max(all_Cindex))  
# scores <- predict(fit, newdata = data_filled_factor)$chf
# 
# 
# # 算一下6, 12, 18, 24, 36个月的进展率
# 
# # 设置时间点
# time_points <- c((3*30+3*31), (6*30+6*31), (9*30+9*31), (12*30+12*31),  (18*30+18*31))
# # 创建一个空的结果数据框
# pfs_progression_rates <- data.frame(Time_Point = time_points, Progression_Rate = NA)
# # 循环计算各时间点的进展率
# for (i in seq_along(time_points)) {
#   # 当前时间点
#   current_time <- time_points[i]
#   
#   # 计算进展率
#   data_filled_factor$Progression_Status <- ifelse(data_filled_factor$PFS <= current_time & data_filled_factor$Progression == 1, 1, 0)
#   progression_rate <- mean(data_filled_factor$Progression_Status, na.rm = TRUE)  # 计算进展率
#   
#   pfs_progression_rates$Progression_Rate[i] <- progression_rate
# }
# # 查看结果
# print(pfs_progression_rates)
# 
# 
# # 汇报给robin之后
# # 把1年、2年、3年的ROC都画出来
# # 假设predictions是模型预测的风险评分
# # 使用累积风险分数
# data_filled_factor_pre <- data_filled_factor
# set.seed(which.max(all_Cindex))
# # 使用RSF模型生成预测结果，type = "risk" 提取风险评分
# risk_scores <- predict(fit, data_filled_factor_pre, type = "risk")$predicted
# 
# # 提取生存时间和状态
# time <- data_filled_factor_pre$PFS
# status <- data_filled_factor_pre$Progression  # 0或1表示无进展或有进展
# 
# # 计算365天、730天和1095天的时间依赖性ROC
# roc_1_year <- timeROC(T = time, delta = status, marker = risk_scores, 
#                       cause = 1, times = 365, iid = TRUE)
# roc_2_years <- timeROC(T = time, delta = status, marker = risk_scores, 
#                        cause = 1, times = 730, iid = TRUE)
# roc_3_years <- timeROC(T = time, delta = status, marker = risk_scores, 
#                        cause = 1, times = 1095, iid = TRUE)
# 
# # 绘制ROC曲线
# plot(roc_1_year, time = 365, col = "blue", lwd = 2, main = "Time-dependent ROC Curves")
# plot(roc_2_years, time = 730, add = TRUE, col = "green", lwd = 2)
# plot(roc_3_years, time = 1095, add = TRUE, col = "red", lwd = 2)
# 
# # 获取AUC值
# auc_1_year <- round(roc_1_year$AUC[2], 3)
# auc_2_years <- round(roc_2_years$AUC[2], 3)
# auc_3_years <- round(roc_3_years$AUC[2], 3)
# # 添加AUC值到图中
# text(0.4, 0.4, paste("1-year AUC =", auc_1_year), col = "blue", cex = 1)
# text(0.4, 0.35, paste("2-year AUC =", auc_2_years), col = "green", cex = 1)
# text(0.4, 0.3, paste("3-year AUC =", auc_3_years), col = "red", cex = 1)
# 
# # 添加图例
# legend("bottomright", legend = c("1-year ROC", "2-year ROC", "3-year ROC"),
#        col = c("blue", "green", "red"), lwd = 2)
# 
# # 提取AUC值
# cat("1-year AUC:", roc_1_year$AUC[2], "\n")
# cat("2-year AUC:", roc_2_years$AUC[2], "\n")
# cat("3-year AUC:", roc_3_years$AUC[2], "\n")
# # # 获取模型的时间点
# # model_times <- fit$time
# # 
# # # 查看pred_risk的第一列对应的时间点
# # first_time_point <- model_times[1]
# # # 提取模型中的时间点
# # model_times <- fit$time
# # 
# # # 计算中位生存时间
# # fit_pfs <- survfit(Surv(PFS, Progression) ~ 1, data = data_filled_factor)
# # summary_fit_pfs <- summary(fit_pfs)
# # median_pfs <- summary_fit_pfs$table["median"]
# # print(median_pfs)
# # 
# # # 找到模型中最接近中位生存时间的时间点索引
# # closest_time_point <- which.min(abs(model_times - median_pfs))
# # 
# # # 提取中位生存时间点的累积风险
# # risk_at_median_time <- scores[, closest_time_point]
# 
# 
# # 选择1年的累积风险（假设PFS单位是天，2年约为730天）
# timepoint <- 1095
# set.seed(which.max(all_Cindex))
# scores_at_timepoint <- scores[, which.min(abs(fit$time.interest - timepoint))]
# 
# #####这里是PFS的
# 
# data_filled_factor4 <- data_filled_factor
# # 将累积风险分数放入数据框中
# data_filled_factor4$scores <- scores_at_timepoint  # 假设选择第一个时间点（例如12个月）的累积风险分数
# 
# # 寻找最佳cutoff值来划分高低风险组
# cutpoint <- surv_cutpoint(data_filled_factor4, 
#                           time = "PFS", 
#                           event = "Progression", 
#                           variables = "scores")
# 
# # 显示最佳cutoff值
# print(cutpoint$cutpoint)
# 
# # 使用找到的cutoff值，将患者划分为高低风险组
# data_filled_factor4$risk_group <- ifelse(data_filled_factor4$scores >= as.numeric(cutpoint$cutpoint[1]), "high", "low")
# 
# # 绘制Kaplan-Meier生存曲线
# # fit_PFS <- survfit(Surv(PFS, Progression) ~ risk_group, data = data_filled_factor4)
# # ggsurvplot(fit_PFS, data = data_filled_factor4, 
# #            risk.table = TRUE, 
# #            pval = TRUE, 
# #            conf.int = TRUE,
# #            legend.labs = c("Low Risk", "High Risk"), 
# #            title = "PFS based on Risk Groups with Optimal Cutoff")
# 
# # 创建一个二分类的PFS变量（例如，一年是否无进展）
# data_filled_factor4$PFS_1yr <- ifelse(data_filled_factor4$PFS >= 365 & data_filled_factor4$Progression == 0, 1, 0)
# 
# # 使用logistic回归模型计算不同风险组之间的OR
# logistic_model <- glm(PFS_1yr ~ risk_group, family = binomial, data = data_filled_factor4)
# summary(logistic_model)
# 
# # 提取OR和置信区间
# exp(cbind(OR = coef(logistic_model), confint(logistic_model)))
# 
# 
# 
# ######下面是OS的
# # wok突然想起来前面为了方便把OS的数据都删掉了，只能在这里再来一遍数据整理了
# cols <- c(3:8, 10:34)
# data_filled_factor5 <- data_filled_[,c(42, 43, 2:33)]
# data_filled_factor5[, cols] <- lapply(data_filled_factor5[, cols], factor)
# 
# # 数值型数据就是数值型数据
# data_filled_factor5$`Max_Diameter` <- as.numeric(data_filled_factor5$`Max_Diameter`)
# data_filled_factor5$OS.time <- as.numeric(data_filled_factor5$OS.time)
# data_filled_factor5$OS <- as.numeric(data_filled_factor5$OS)
# sum(is.na(data_filled_factor5))
# 
# # data_filled_factor5 <- data_filled_factor5[,9:49]
# # 将累积风险分数放入数据框中
# data_filled_factor5$scores <- scores_at_timepoint  # 假设选择第一个时间点（例如12个月）的累积风险分数
# 
# # 寻找最佳cutoff值来划分高低风险组
# cutpoint <- surv_cutpoint(data_filled_factor5, 
#                           time = "OS.time", 
#                           event = "OS", 
#                           variables = "scores")
# 
# # 显示最佳cutoff值
# print(cutpoint$cutpoint)
# 
# # 使用找到的cutoff值，将患者划分为高低风险组
# data_filled_factor5$risk_group <- ifelse(data_filled_factor5$scores >= as.numeric(cutpoint$cutpoint[1]), "high", "low")
# 
# # 绘制Kaplan-Meier生存曲线
# # fit_PFS <- survfit(Surv(PFS, Progression) ~ risk_group, data = data_filled_factor5)
# # ggsurvplot(fit_PFS, data = data_filled_factor5, 
# #            risk.table = TRUE, 
# #            pval = TRUE, 
# #            conf.int = TRUE,
# #            legend.labs = c("Low Risk", "High Risk"), 
# #            title = "PFS based on Risk Groups with Optimal Cutoff")
# 
# # 创建一个二分类的PFS变量（例如，一年是否无死亡）
# data_filled_factor5$OS_1yr <- ifelse(data_filled_factor5$OS.time >= 365 & data_filled_factor5$OS == 0, 1, 0)
# 
# # 使用logistic回归模型计算不同风险组之间的OR
# logistic_model <- glm(OS_1yr ~ risk_group, family = binomial, data = data_filled_factor5)
# summary(logistic_model)
# 
# # 提取OR和置信区间
# exp(cbind(OR = coef(logistic_model), confint(logistic_model)))