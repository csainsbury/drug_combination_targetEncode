## cbg data for xgboost classification

## CBG spar plots
library(data.table)
library(tidyverse)
library(tidyquant)
library(timetk)
library(padr)
library(zoo)
library(viridis)
library(imputeTS)
library(xgboost)
library(caret)
library(pROC)

# prediction period (15 min increment values)
hours <- 8
prediction_bin <- hours * 4
gap_hours <- 0.5
gap_bins <- gap_hours * 4
downsample <- 1

sup_path <- c('~/Documents/data/libre_as_supervised/')
file_list <- list.files(paste0(sup_path))
n_ids = 15 # need automated way of determining number of IDs data

# n_vec <- c('al1', 'ak1', 'at1','cd1','hm1','mv1','ps1','se1','ek2')
load_train <- 384 # 192 # 96 # 384
load_test <- 96

for (n in c(1:n_ids)) {
  a <- fread(paste0('~/Documents/data/libre_as_supervised/supervised_', n,'_', load_train, '_', load_test, '_vals.csv'))
  print(nrow(a))
  if (n == 1) {
    x <- a
  } else {
    x <- rbind(x, a)
  }
}

t_col = (ncol(x) - load_test)

train <- x[, 1:t_col]

if (gap_hours > 0) {
  test <- x[, (t_col + 1 + gap_bins):(t_col + prediction_bin + gap_bins)]
} else {
  test <- x[, (t_col + 1):(t_col + prediction_bin)]
}

min_test <- apply(test,1,min)
label <- ifelse(min_test <3, 1, 0)

train$label <- label
train$id <- c(1:nrow(train))

# take every nth row (to prevent over similarity)
n = prediction_bin
train <- train[id%%n==0]

# down sample n
case_n = nrow(train[label ==  1]) # n cases
ratio = 1

cases    <- train[label == 1]
controls <- train[label == 0]

set.seed(42)
select_case <- cases[sample(nrow(cases), case_n), ]
select_control <- controls[sample(nrow(controls), case_n * ratio, replace = F), ]

sc <- rbind(cases, select_control)
sc <- sc[sample(nrow(sc), round(nrow(sc) * downsample, 0)), ]

write.table(sc, file = paste0('~/Documents/data/libreTSdata/libre_', load_train, '_', prediction_bin,'.csv'),
            sep = ',', row.names = F)


## generate data to produce per day features for xgb model

## convert to long
u_ids <- unique(sc$id)
for (i in c(1:length(u_ids))) {
  
  id = u_ids[i]
  sub <- sc[id == u_ids[i]]
  label <- sub$label
  sub <- sub[, -c('id', 'label')]
  
  new <- data.table(t(sub))
  colnames(new) <- c('glu')
  
  new$timestep <- c(1:nrow(new))
  new$id <- id
  new$label <- label
  
  for (j in c(1:(ncol(sub) / 96))) {
    c <- rep(j, 96)
    if (j == 1) {
      cc = c
    } else {
      cc = c(cc, c)
    }
  }
  
  new$day_label = cc
  
  if (i == 1) {
    newcat = new
  } else {
    newcat = rbind(newcat, new)
  }
  
}

s <- newcat

s[, 'max_by_day' := max(glu), by=.(id, day_label)]
s[, 'min_by_day' := min(glu), by=.(id, day_label)]
s[, 'median_by_day' := median(glu), by=.(id, day_label)]
s[, 'iqr_by_day' := quantile(glu)[4]-quantile(glu)[2], by=.(id, day_label)]
s[, 'gradient_by_day' := as.numeric(lm(glu ~ timestep)$coefficients[2]), by=.(id, day_label)]

s[, 'day_n' := c(1:.N), by=.(id, day_label)]
s[, 'gradient' := as.numeric(lm(glu ~ timestep)$coefficients[2]), by=.(id)]
s[, 'sd_cbg' := sd(glu) ,by=.(id)]
s[, 'mean_cbg' := mean(glu) ,by=.(id)]
s[, 'cV' := sd_cbg/mean_cbg ,by=.(id)]

s[, 'ID_n' := c(1:.N), by=.(id)]

hist(s[ID_n==1]$cV, 100)

## per ID generate wide data
ids <-unique(s$id)
for (i in c(1:length(ids))) {
  
  if (i %% 100 == 0) {print(i/length(ids))}
  
  sub <- s[id == ids[i]]
  # sub <- sub[day_n == 1]
  
  sub <- sub %>% select(id, glu, day_label, timestep, max_by_day, min_by_day, median_by_day, iqr_by_day, gradient_by_day, cV, gradient, label)
    
    static <- sub %>% select(id, gradient, cV, label)
    
    byday <- sub %>% select(day_label, max_by_day, min_by_day, median_by_day, iqr_by_day, gradient_by_day)
    byday[, 'day_n' := c(1:.N), by=.(day_label)]
    byday <- byday[day_n == 1]
    byday <- byday[, -c('day_n')]
    
      dayMelt <- melt(byday, id.vars='day_label')
      dm <- as.data.frame(t(dayMelt))
      dm <- dm[-1, ]
      colnames(dm) <- c(dm[1,])
      dm <- dm[-1, ]
    
    dm_add_static <- cbind(dm, static[1,])
    
    glud <- sub %>% select(glu)
    glud$index <- c(1:nrow(glud))
    
    glud$dg <- 0
    for (k in c(1:(nrow(glud) - 1))) {
      grad_vec <- glud$glu[k:nrow(glud)]
      steps <- c(k:nrow(glud))
      glud$dg[k] <- as.numeric(lm(grad_vec ~ steps)$coefficients[2])
    }
    
      gluMelt <- melt(glud, id.vars='index')
      gm <- as.data.frame(t(gluMelt))
      gm <- gm[-1, ]
      colnames(gm) <- c(gm[1,])
      gm <- gm[-1, ]
      
    final <- cbind(gm, dm_add_static)
    
    if (i == 1) {
      export <- final
    } else {
      export <- rbind(export, final)
    }
    
}

coln = colnames(export)
glu_names <- coln[coln == 'glu']
dg_names <- coln[coln == 'dg']
max_names <- coln[coln == 'max_by_day']
min_names <- coln[coln == 'min_by_day']
med_names <- coln[coln == 'median_by_day']
iqr_names <- coln[coln == 'iqr_by_day']
gradientByDay_names <- coln[coln == 'gradient_by_day']

other_names <- coln[coln != 'glu' & coln!='dg' & coln != 'max_by_day'
                    & coln != 'min_by_day' & coln != 'iqr_by_day' & coln != 'median_by_day' &
                      coln != 'gradient_by_day']
for (c in c(1:length(glu_names))) {
  glu_names[c] <- paste0(glu_names[c], '_', c)
}
for (c in c(1:length(dg_names))) {
  dg_names[c] <- paste0(dg_names[c], '_', c)
}
for (c in c(1:length(max_names))) {
  max_names[c] <- paste0(max_names[c], '_', c)
}
for (c in c(1:length(min_names))) {
  min_names[c] <- paste0(min_names[c], '_', c)
}
for (c in c(1:length(iqr_names))) {
  iqr_names[c] <- paste0(iqr_names[c], '_', c)
}
for (c in c(1:length(med_names))) {
  med_names[c] <- paste0(med_names[c], '_', c)
}
for (c in c(1:length(gradientByDay_names))) {
  gradientByDay_names[c] <- paste0(gradientByDay_names[c], '_', c)
}

coln_n <- c(glu_names, dg_names, max_names, min_names, med_names, iqr_names, gradientByDay_names, other_names)
colnames(export) <- coln_n

for (j in c(1:ncol(export))) {
  export[, j] <- as.numeric(export[, j])
}
# 
write.table(export, file = paste0('Documents/data/libre_data/libre_export/v2_', n_ids,'_IDs_trainBins_', load_train, '_predictionBins_', prediction_bin, '_gapBins_', gap_bins,'.csv'), sep = ',', row.names = F)
# 
# export1 <- fread(paste0('Documents/data/CBGdata/abstract_exports/export_admissionDuration_', days_n, '_days.csv'))

##### build model
##
set.seed(42)
export = data.table(export)
k = 4
flds <- createFolds(export$label, k, list = T, returnTrain = F)

auc_vec <- rep(0, k)
n_vec <- rep(0, k)
threshold_frame <- as.data.frame(matrix(nrow = k, ncol = 3))

accuracy_v <- rep(0, k)
balanced_accuracy_v <- rep(0, k)
f1 <- rep(0, k)

maxd = 8

for (kfold in c(1:k)) {
  
  train = export[-flds[[kfold]], ]
  test  = export[flds[[kfold]], ]
  
  #define predictor and response variables in training set
  train_x = data.matrix(train[, -c('id', 'label')])
  train_y = train[,'label']
  
  #define predictor and response variables in testing set
  test_x = data.matrix(test[, -c('id', 'label')])
  test_y = test[,'label']
  
  #define final training and testing sets
  xgb_train = xgb.DMatrix(data = as.matrix(train_x), label = train_y$label)
  xgb_test = xgb.DMatrix(data = as.matrix(test_x), label = test_y$label)
  
  #define watchlist
  watchlist = list(train=xgb_train, test=xgb_test)
  
  #fit XGBoost model and display training and testing data at each round
  param <- list(max.depth = maxd, eta = 0.1, nthread = 64, gamma = 2, min_child_weight = 1, objective = 'binary:logistic')
  model = xgb.train(param, data = xgb_train, watchlist=watchlist,
                    # lambda = 1000,
                    nrounds = 100, verbose = 1)
  
  # model = xgb.train(data = xgb_train, max.depth = maxd, watchlist=watchlist, nrounds = 100, nthread = 16)
  # plot(model$evaluation_log$test_rmse)
  
  # n = which(model$evaluation_log$test_rmse == min(model$evaluation_log$test_rmse))
  n = which(model$evaluation_log[, 3] == min(model$evaluation_log[, 3]))
  print(n)
  n_vec[kfold] <- n
  
  final = xgboost(data = xgb_train, max.depth = maxd, nrounds = n, verbose = 1)
  
  pred_y = predict(final, xgb_test)
  # hist(pred_y, 100)
  
  ##
  pROC_obj <- roc(test_y$label,pred_y,
                  # arguments for ci
                  ci=TRUE,
                  # arguments for plot
                  plot=FALSE, grid=TRUE,
                  print.auc=TRUE)
  
  auc_vec[kfold] <- pROC_obj$auc
  
  c = coords(pROC_obj, "best", "threshold")
  threshold_frame[kfold, ] = c
  
  print(pROC_obj$auc)
  
  if (kfold == 1) {
    plot(pROC_obj$specificities, pROC_obj$sensitivities, xlim = c(1, 0), cex = 0)
    lines(pROC_obj$specificities, pROC_obj$sensitivities, col = rgb(0,0,0,0.6, maxColorValue = 1), lwd=2)
  } else{
    points(pROC_obj$specificities, pROC_obj$sensitivities, cex = 0)
    lines(pROC_obj$specificities, pROC_obj$sensitivities, col = rgb(0,0,0,0.6, maxColorValue = 1), lwd=2)
  }
  
  ## generate validation accuracy metric
  pred_acc <- ifelse(pred_y >= c$threshold, 1, 0)
  cm <- confusionMatrix(as.factor(pred_acc), as.factor(test_y$label))
  
  accuracy_v[kfold] <- cm$overall[1]
  balanced_accuracy_v[kfold] <- cm$byClass[11]
  f1[kfold] <- cm$byClass[7]
  
  outputTable <- data.table('pred' = pred_y, 'ref' = test_y$label)
  
  importance_matrix = xgb.importance(colnames(xgb_train), model = model)
  
  if (kfold == 1) {
    outputTable_cat <- outputTable
  } else {
    outputTable_cat <- rbind(outputTable_cat, outputTable)
  }
  
  # if (kfold == k) {
  #   legend(0.2, 0.8, legend=paste0('kf auroc: ', round(auc_vec, 4)), cex=0.8)
  #   legend(0.2, 0.2, legend=paste0('t/sen/spec: ', round(apply(threshold_frame, 2, median), 2)), cex=0.8)
  # }
  
  
}

print(quantile(auc_vec))
print(quantile(n_vec))
print(threshold_frame)

print(quantile(accuracy_v))
print(quantile(balanced_accuracy_v))
print(quantile(f1))

print(quantile(threshold_frame[,2]))
print(quantile(threshold_frame[,3]))

importance_matrix = xgb.importance(colnames(xgb_train), model = model)
xgb.plot.importance(importance_matrix)

imp <- importance_matrix[order(-importance_matrix$Importance), ]
imp[1:20,]

hist(outputTable_cat$pred, 400)


##################################
# standard metrics


base_predictions_prevHypo <- function(ID, Glu, label, hypothresh) {
  
  # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
  d <- data.table(ID, Glu, label)
  d[, 'prior_below_3' := ifelse(min(Glu) < hypothresh, 1, 0) ,by=.(ID)]
  d[, 'n' := c(1: .N), by=.(ID)]
  single_d <- d[n==1]
  print(confusionMatrix(as.factor(single_d$prior_below_3),
                        as.factor(single_d$label),
                        positive = '1',
                        prevalence = 0.33))
  
  cm = confusionMatrix(as.factor(single_d$prior_below_3),
                       as.factor(single_d$label),
                       positive = '1',
                       prevalence = 0.33)
  
  return(list(as.numeric(cm$overall[1]), # accuracy
              as.numeric(cm$byClass[11]), # balanced accuracy
              as.numeric(cm$byClass[7]),# f1
              as.numeric(cm$byClass[1]),# sens
              as.numeric(cm$byClass[2]))) # spec
  
}

base_predictions_cv <- function(ID, Glu, label) {
  
  # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
  d <- data.table(ID, Glu, label)
  
  d[, 'sd_cbg' := sd(Glu) ,by=.(ID)]
  d[, 'mean_cbg' := mean(Glu) ,by=.(ID)]
  d[, 'cV' := sd_cbg/mean_cbg ,by=.(ID)]
  
  d[, 'n' := c(1: .N), by=.(ID)]
  
  single_d <- d[n==1]
  
  require(pROC)
  roc_obj <- roc(single_d$label, single_d$cV)
  print(auc(roc_obj))
  
  bestthresh <- as.numeric(coords(roc_obj, "best", ret="threshold", transpose = FALSE)[1,])
  print(confusionMatrix(as.factor(ifelse(single_d$cV > bestthresh, 1, 0)),
                        as.factor(single_d$label),
                        positive = '1',
                        prevalence = 0.33))
  cm = confusionMatrix(as.factor(ifelse(single_d$cV > bestthresh, 1, 0)),
                       as.factor(single_d$label),
                       positive = '1',
                       prevalence = 0.33)
  
  return(list(as.numeric(auc(roc_obj)),
              as.numeric(cm$overall[1]), # accuracy
              as.numeric(cm$byClass[11]), # balanced accuracy
              as.numeric(cm$byClass[7]),# f1
              as.numeric(cm$byClass[1]),# sens
              as.numeric(cm$byClass[2]))) # spec
  
  
}

#########
## test by bootstrapping

single_n <- x[n==1]
ids_for_split <- data.table('id' = single_n$ID,
                            'label' = single_n$label)

bootst <- 1000
bootst_n <- 100

predhV_acc <- rep(0, bootst)
predhV_bal_acc <- rep(0, bootst)
predhV_f1 <- rep(0, bootst)
predhV_sens <- rep(0, bootst)
predhV_spec <- rep(0, bootst)
aucV <- rep(0, bootst)
predcvV_acc <- rep(0, bootst)
predcvV_bal_acc <- rep(0, bootst)
predcvV_f1 <- rep(0, bootst)
predcvV_sens <- rep(0, bootst)
predcvV_spec <- rep(0, bootst)

for (j in c(1:bootst)) {
  
  sub <- sample(nrow(ids_for_split), bootst_n)
  bootstrap_ids <- ids_for_split[sub]
  
  dat <- merge(bootstrap_ids, x, by.x = 'id', by.y = 'ID')
  
  predh <- base_predictions_prevHypo(dat$id, dat$Glu, dat$label.x, 3)
  predcv<- base_predictions_cv(dat$id, dat$Glu, dat$label.x)
  
  predhV_acc[j] <- predh[[1]]
  predhV_bal_acc[j] <- predh[[2]]
  predhV_f1[j] <- predh[[3]]
  predhV_sens[j] <- predh[[4]]
  predhV_spec[j] <- predh[[5]]
  
  aucV[j] <- predcv[[1]]
  predcvV_acc[j] <- predcv[[2]]
  predcvV_bal_acc[j] <- predcv[[3]]
  predcvV_f1[j] <- predcv[[4]]
  predcvV_sens[j] <- predcv[[5]]
  predcvV_spec[j] <- predcv[[6]]
}

quantile(predhV_acc)
quantile(predhV_bal_acc)
quantile(predhV_f1)
quantile(predhV_sens)
quantile(predhV_spec)

quantile(aucV)

quantile(predcvV_acc)
quantile(predcvV_bal_acc)
quantile(predcvV_f1)
quantile(predcvV_sens)
quantile(predcvV_spec)


