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

# 101355572

dataset <- 'huge'
days_n  <- 2
minimum_n_cbgs <- days_n + 2
hypo_threshold <- 3

x <- fread(paste0('~/Documents/data/CBGdata/', dataset, '_unipoc_time_series_cohort_first_', days_n,'_days.csv'))
x <- x[order(x$uID, x$admission_vec)]

exclude <- c('dialysis', 'renal', 'abbott', 'Ketones', 'X Ray', 'X-ray', 'training', 'staff', 'returned', 'obsolete', 'old', 'Meterx', 'dead', 'BioChemistry', 'weldon', 'radiology', 'outpatient', 'rcas', 'Labour Ward', 'Neonatal Unit Top Floor')
x_ex <- x
for (e in c(1:length(exclude))) {
  x_ex <- x_ex[-grep(exclude[e], x_ex$loc, ignore.case = TRUE), ]
}
x <- x_ex

#library(dplyr)
# add unique ID by 2 cols - uID and admission vec
x <- x %>% mutate(ID = group_indices(x, .dots=c("uID", "admission_vec"))) 
# x[, 'N' := .N, by=.(uID)]

## here need to truncate all admissions to the right length!!
x$day <- as.Date(x$dateTime)
x[, 'correct_duration_days' := ifelse(day <= (min(day) + (days_n - 1)), 1, 0), by=.(ID)]
x <- x[correct_duration_days == 1]

# limit to those with minimum n CBGs
x[, 'N_truncated' := .N, by=.(ID)]
x <- x[N_truncated >= minimum_n_cbgs]

# split out last day and label
#t <- x[V1 == 203446038]
#x$day <- as.Date(x$dateTime)
x[, 'flag_last_day' := ifelse(day == max(day), 1, 0), by=.(ID)]
x[, 'label' := ifelse(min(Glu[flag_last_day == 1]) <= hypo_threshold, 1, 0), by=.(ID)]
# remove the last day from the training set
x <- x[flag_last_day == 0]

x[, 'N' := .N, by=.(ID)]
x <- x[N>days_n]

x[, 'n' := c(1 : .N), by=.(ID)]
print(sum(x[n==1]$label))
print(nrow(x[n==1]))

# plot(x$dateTime, x$Glu, pch = 16, cex = 0.6, col = ifelse(x$label == 0,
#                                                           rgb(0,0,0,0.4, maxColorValue = 1),
#                                                           rgb(1,0,0,0.4, maxColorValue = 1)))

# produce fixed ratio of case to no case
ratio = 5
event_ids <- unique(x[label==1]$ID)
no_event_ids <- unique(x[label==0]$ID)
id_sample <- no_event_ids[sample(length(no_event_ids), round(length(event_ids) * ratio), 0)]

neg <- x[ID %in% id_sample]
pos <- x[ID %in% event_ids]

x <- rbind(neg, pos)

print(sum(x[n==1]$label))
print(nrow(x[n==1]))

s <- x

s[, 'max_by_day' := max(Glu), by=.(ID, day)]
s[, 'min_by_day' := min(Glu), by=.(ID, day)]
s[, 'median_by_day' := median(Glu), by=.(ID, day)]
s[, 'iqr_by_day' := quantile(Glu)[4]-quantile(Glu)[2], by=.(ID, day)]

s[, 'day_n' := c(1:.N), by=.(ID, day)]
s[, 'day_N' := .N, by=.(ID, day)]
s[, 'gradient' := as.numeric(lm(Glu ~ dateTime)$coefficients[2]), by=.(ID)]
    s[, 'sd_cbg' := sd(Glu) ,by=.(ID)]
    s[, 'mean_cbg' := mean(Glu) ,by=.(ID)]
    s[, 'cV' := sd_cbg/mean_cbg ,by=.(ID)]

s[, 'ID_n' := c(1:.N), by=.(ID)]
s[, 'ID_N' := .N, by=.(ID)]

hist(s[ID_n==1]$cV, 1000)

s[, 'sum_diff_day' := as.numeric(sum(diff(day))), by=.(ID)]

## per ID generate wide data
ids <-unique(s$ID)
for (i in c(1:length(ids))) {
  
  if (i %% 100 == 0) {print(i/length(ids))}
  
  sub <- s[ID == ids[i]]
  sub <- sub[day_n == 1]
  
  sub <- sub %>% select(ID, dateTime, day, Glu, max_by_day, min_by_day, median_by_day, day_N, iqr_by_day, cV, gradient, label)
  
  if ((nrow(sub) > 2) & (as.numeric(min(diff(as.Date(sub$dateTime)))) <= 1)) { # second arguement needed to ensure that minimum time step is not greater than 1 day (thicken will no work if it is)
    
    #sub$day_numeric <- c(1:nrow(sub))
    
    if ((as.Date(sub$dateTime[nrow(sub)]) - as.Date(sub$dateTime[1]) < (days_n - 1))) {
      sub <- as.data.frame(sub)
      sub[nrow(sub)+1,] <- NA
      sub <- as.data.table(sub)
      
      sub$dateTime[nrow(sub)] <- (sub$dateTime[1] + ((days_n - 1) * 60*60*24))
      
    }
    
    cV <- sub$cV[1]
    gradient <- sub$gradient[1]
    id <- sub$ID[1]
    label <- sub$label[1]
    
    sub <- sub %>% select(dateTime, Glu, max_by_day, min_by_day, median_by_day, day_N, iqr_by_day)
    sub <- sub %>% thicken("1 day") %>% select(-dateTime) %>% pad()
    
    # sub.locf <- na.locf(sub)
    sub.locf <- sub
    sub.locf$Glu <- NULL
    
    require(reshape2)
    melt.sub <- melt(sub.locf, id.vars='dateTime_day')
    m <- as.data.frame(t(melt.sub))
    m$cV <- cV
    m$gradient <- gradient
    m$id <- id
    m$label <- label
    m <- m[-1, ]
    colnames(m) <- c(m[1,1:(ncol(m) - 4)],'cV' ,'gradient' ,'id', 'label')
    m <- m[-1, ]
    
    if (i == 1) {
      export <- m
    } else {
      export <- rbind(export, m)
    }
    
  }
  
  
}

exp <- export
coln = colnames(exp)
max_names <- coln[coln == 'max_by_day']
min_names <- coln[coln == 'min_by_day']
med_names <- coln[coln == 'median_by_day']
N_names <- coln[coln == 'day_N']
iqr_names <- coln[coln == 'iqr_by_day']

other_names <- coln[coln != 'max_by_day' & coln != 'min_by_day' & coln != 'iqr_by_day' &
                      coln != 'median_by_day' & coln != 'day_N' ]

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
for (c in c(1:length(N_names))) {
  N_names[c] <- paste0(N_names[c], '_', c)
}

coln_n <- c(max_names, min_names, med_names, N_names, iqr_names, other_names)
colnames(exp) <- coln_n

export <- exp

for (j in c(1:ncol(export))) {
  export[, j] <- as.numeric(export[, j])
}

write.table(export, file = paste0('Documents/data/CBGdata/abstract_exports/export_admissionDuration_', days_n, '_days_hypothresh_NAs_included_', hypo_threshold, '.csv'), sep = ',', row.names = F)

export1 <- fread(paste0('Documents/data/CBGdata/abstract_exports/export_admissionDuration_', days_n, '_days.csv'))

##### build model
##
export = data.table(export)
k = 10
flds <- createFolds(export$label, k, list = T, returnTrain = F)

auc_vec <- rep(0, k)
n_vec <- rep(0, k)
threshold_frame <- as.data.frame(matrix(nrow = k, ncol = 3))

accuracy_v <- rep(0, k)
balanced_accuracy_v <- rep(0, k)
f1 <- rep(0, k)

maxd = 2

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
  param <- list(max.depth = maxd, eta = 0.1, nthread = 64, min_child_weight = 1)
  model = xgb.train(param, data = xgb_train, watchlist=watchlist, nrounds = 100, verbose = 0)
  
  # model = xgb.train(data = xgb_train, max.depth = maxd, watchlist=watchlist, nrounds = 100, nthread = 16)
  # plot(model$evaluation_log$test_rmse)
  
  n = which(model$evaluation_log$test_rmse == min(model$evaluation_log$test_rmse))
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
      importance_matrix_merge <- importance_matrix
    } else {
      outputTable_cat <- rbind(outputTable_cat, outputTable)
      # importance_matrix_merge = merge(importance_matrix_merge, importance_matrix, by.x = 'Feature', by.y = 'Feature')
      
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
importance_matrix <- importance_matrix[order(-importance_matrix$Importance), ]
importance_matrix[1:20,]

hist(outputTable_cat$pred, 100)


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
bootst_n <- 1000

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


