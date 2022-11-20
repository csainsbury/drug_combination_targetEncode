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
library(ParBayesianOptimization)

# 101355572
x <- fread('~/Documents/data/glucose_steroid/new_final.csv')

export = x

##### build model - bayes opt
splits <- rsample::initial_split(export, prop = 0.8)
# Training set
train_df <- rsample::training(splits)
# Test set
test_df <- rsample::testing(splits)
dim(train_df)
## [1] 354  23
dim(test_df)
## [1] 152  23

# The xgboost interface accepts matrices 
X <- train_df %>%
  # Remove the target variable
  select(!label_max) %>%
  as.matrix()

# Get the target variable
y <- train_df %>%
  pull(label_max)

# Cross validation folds
folds <- list(fold1 = as.integer(seq(1, nrow(X), by = 3)),
              fold2 = as.integer(seq(2, nrow(X), by = 3)))

# Function must take the hyper-parameters as inputs
obj_func <- function(eta, max_depth, min_child_weight, subsample, lambda, alpha) {
  
  param <- list(
    
    # Hyter parameters 
    eta = eta,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    subsample = subsample,
    lambda = lambda,
    alpha = alpha,
    
    # Tree model 
    booster = "gbtree",
    
    # Regression problem 
    objective = "binary:logistic",
    
    # Use the Mean Absolute Percentage Error
    eval_metric = "logloss")
  
  xgbcv <- xgb.cv(params = param,
                  data = X,
                  label = y,
                  nround = 100,
                  folds = folds,
                  prediction = TRUE,
                  early_stopping_rounds = 5,
                  verbose = 1,
                  maximize = F)
  
  lst <- list(

    # First argument must be named as "Score"
    # Function finds maxima so inverting the output
    Score = -min(xgbcv$evaluation_log$test_logloss_mean),

    # Get number of trees for the best performing model
    nrounds = xgbcv$best_iteration
  )
  
  return(lst)
}

bounds <- list(eta = c(0.001, 0.6),
               max_depth = c(1L, 20L),
               min_child_weight = c(1, 100),
               subsample = c(0.01, 1),
               lambda = c(1, 100),
               alpha = c(1, 100))
set.seed(1234)
bayes_out <- bayesOpt(FUN = obj_func, bounds = bounds, initPoints = length(bounds) + 2, iters.n = 80)

bayes_out$scoreSummary
data.frame(getBestPars(bayes_out))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# export = data.table(export)
# k = 3
# flds <- createFolds(export$y, k, list = T, returnTrain = F)
# 
# auc_vec <- rep(0, k)
# n_vec <- rep(0, k)
# threshold_frame <- as.data.frame(matrix(nrow = k, ncol = 3))
# 
# accuracy_v <- rep(0, k)
# balanced_accuracy_v <- rep(0, k)
# f1 <- rep(0, k)
# 
# maxd = 2
# 
# for (kfold in c(1:k)) {
#   
#   train = export[-flds[[kfold]], ]
#   test  = export[flds[[kfold]], ]
#   
#   #define predictor and response variables in training set
#   train_x = data.matrix(train[, -c('V1', 'y')])
#   train_y = train[,'y']
#   
#   #define predictor and response variables in testing set
#   test_x = data.matrix(test[, -c('V1', 'y')])
#   test_y = test[,'y']
#   
#   #define final training and testing sets
#   xgb_train = xgb.DMatrix(data = as.matrix(train_x), label = train_y$y)
#   xgb_test = xgb.DMatrix(data = as.matrix(test_x), label = test_y$y)
#   
#   #define watchlist
#   watchlist = list(train=xgb_train, test=xgb_test)
#   
#   #fit XGBoost model and display training and testing data at each round
#   param <- list(max.depth = maxd, eta = 0.1, nthread = 64, min_child_weight = 1)
#   model = xgb.train(param, data = xgb_train, watchlist=watchlist, nrounds = 40, verbose = 0,
#                     subsample = 1, lambda = 1, alpha = 1)
#   
#   # model = xgb.train(data = xgb_train, max.depth = maxd, watchlist=watchlist, nrounds = 100, nthread = 16)
#   # plot(model$evaluation_log$test_rmse)
#   
#   n = which(model$evaluation_log$test_rmse == min(model$evaluation_log$test_rmse))
#   print(n)
#   n_vec[kfold] <- n
#   
#   final = xgboost(data = xgb_train, max.depth = maxd, nrounds = n, verbose = 1)
#   
#   pred_y = predict(final, xgb_test)
#   # hist(pred_y, 100)
#   
#   ##
#   pROC_obj <- roc(test_y$y,pred_y,
#                   # arguments for ci
#                   ci=TRUE,
#                   # arguments for plot
#                   plot=FALSE, grid=TRUE,
#                   print.auc=TRUE)
#   
#   auc_vec[kfold] <- pROC_obj$auc
#   
#   c = coords(pROC_obj, "best", "threshold")
#   threshold_frame[kfold, ] = c
#   
#   print(pROC_obj$auc)
#   
#   if (kfold == 1) {
#     plot(pROC_obj$specificities, pROC_obj$sensitivities, xlim = c(1, 0), cex = 0)
#     lines(pROC_obj$specificities, pROC_obj$sensitivities, col = rgb(0,0,0,0.6, maxColorValue = 1), lwd=2)
#   } else{
#     points(pROC_obj$specificities, pROC_obj$sensitivities, cex = 0)
#     lines(pROC_obj$specificities, pROC_obj$sensitivities, col = rgb(0,0,0,0.6, maxColorValue = 1), lwd=2)
#   }
#   
#   ## generate validation accuracy metric
#   pred_acc <- ifelse(pred_y >= c$threshold, 1, 0)
#   cm <- confusionMatrix(as.factor(pred_acc), as.factor(test_y$y))
#   
#   accuracy_v[kfold] <- cm$overall[1]
#   balanced_accuracy_v[kfold] <- cm$byClass[11]
#   f1[kfold] <- cm$byClass[7]
#   
#     outputTable <- data.table('pred' = pred_y, 'ref' = test_y$label)
#     importance_matrix = xgb.importance(colnames(xgb_train), model = model)
#     importance_matrix = data.table('Feature' = importance_matrix$Feature,
#                                    'Gain' = importance_matrix$Gain)
#     
#     if (kfold == 1) {
#       outputTable_cat <- outputTable
#       colnames(importance_matrix)[2:ncol(importance_matrix)] <- paste0(colnames(importance_matrix)[2:ncol(importance_matrix)], '_', kfold)
#       imp_merge <- importance_matrix
#     } else {
#       outputTable_cat <- rbind(outputTable_cat, outputTable)
#       colnames(importance_matrix)[2:ncol(importance_matrix)] <- paste0(colnames(importance_matrix)[2:ncol(importance_matrix)], '_', kfold)
#       #imp_merge <- merge(imp_merge, importance_matrix, by.x = 'Feature', by.y = 'Feature', all.x = T, allow.cartesian = T)
#       # importance_matrix_merge = merge(importance_matrix_merge, importance_matrix, by.x = 'Feature', by.y = 'Feature')
#       
#     }
#   
#   # if (kfold == k) {
#   #   legend(0.2, 0.8, legend=paste0('kf auroc: ', round(auc_vec, 4)), cex=0.8)
#   #   legend(0.2, 0.2, legend=paste0('t/sen/spec: ', round(apply(threshold_frame, 2, median), 2)), cex=0.8)
#   # }
#   
#   
# }
# 
# print(quantile(auc_vec))
# print(quantile(n_vec))
# print(threshold_frame)
# 
# print(quantile(accuracy_v))
# print(quantile(balanced_accuracy_v))
# print(quantile(f1))
# 
# print(quantile(threshold_frame[,2]))
# print(quantile(threshold_frame[,3]))
# 
# importance_matrix = xgb.importance(colnames(xgb_train), model = model)
# xgb.plot.importance(importance_matrix)
# importance_matrix <- importance_matrix[order(-importance_matrix$Importance), ]
# importance_matrix[1:20,]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## leave one out

export = data.table(export)

k = nrow(export)

prediction <- rep(0, k)
ground_t <- rep(0, k)

n_vec <- rep(0, k)

maxd = 2

for (kfold in c(1:nrow(export))) {
  
  train = export[-kfold,]
  test  = export[kfold,]
  
  #define predictor and response variables in training set
  train_x = data.matrix(train[, -c('uID', 'label_max')])
  train_y = train[,'label_max']
  
  #define predictor and response variables in testing set
  test_x = data.matrix(test[, -c('uID', 'label_max')])
  test_y = test[,'label_max']
  
  #define final training and testing sets
  xgb_train = xgb.DMatrix(data = as.matrix(train_x), label = train_y$label_max)
  xgb_test = xgb.DMatrix(data = as.matrix(test_x), label = test_y$label_max)
  
  #define watchlist
  watchlist = list(train=xgb_train, test=xgb_test)
  
  #fit XGBoost model and display training and testing data at each round
  param <- list(max.depth = maxd, eta = 0.3, nthread = 64, min_child_weight = 1)
  model = xgb.train(param, data = xgb_train, watchlist=watchlist, nrounds = 40, verbose = 0,
                    subsample = 1, lambda = 1, alpha = 1)
  
  # model = xgb.train(data = xgb_train, max.depth = maxd, watchlist=watchlist, nrounds = 100, nthread = 16)
  # plot(model$evaluation_log$test_rmse)
  
  n = which(model$evaluation_log$test_rmse == min(model$evaluation_log$test_rmse))
  print(n)
  n_vec[kfold] <- n
  
  final = xgboost(data = xgb_train, max.depth = maxd, nrounds = n, verbose = 1)
  
  pred_y = predict(final, xgb_test)
  
  prediction[kfold] = pred_y
  ground_t[kfold] = test_y$label_max
  # hist(pred_y, 100)
  
}
  
print(quantile(n_vec))
  ##
  pROC_obj <- roc(ground_t,prediction,
                  # arguments for ci
                  ci=TRUE,
                  # arguments for plot
                  plot=FALSE, grid=TRUE,
                  print.auc=TRUE)
  
  #auc_vec[kfold] <- pROC_obj$auc
  
  c = coords(pROC_obj, "best", "threshold")
  threshold = c$threshold
  
  print(pROC_obj$auc)
  
    plot(pROC_obj$specificities, pROC_obj$sensitivities, xlim = c(1, 0), cex = 0)
    lines(pROC_obj$specificities, pROC_obj$sensitivities, col = rgb(0,0,0,0.6, maxColorValue = 1), lwd=2)

  ## generate validation accuracy metric
  pred_acc <- ifelse(prediction >= c$threshold, 1, 0)
  cm <- confusionMatrix(as.factor(pred_acc), as.factor(ground_t), positive = '1')
  
  importance_matrix = xgb.importance(colnames(xgb_train), model = model)
  xgb.plot.importance(importance_matrix)
  importance_matrix <- importance_matrix[order(-importance_matrix$Importance), ]
  importance_matrix[1:20,]
  
#     #imp_merge <- merge(imp_merge, importance_matrix, by.x = 'Feature', by.y = 'Feature', all.x = T, allow.cartesian = T)
#     # importance_matrix_merge = merge(importance_matrix_merge, importance_matrix, by.x = 'Feature', by.y = 'Feature')
#     
#   }
#   
#   # if (kfold == k) {
#   #   legend(0.2, 0.8, legend=paste0('kf auroc: ', round(auc_vec, 4)), cex=0.8)
#   #   legend(0.2, 0.2, legend=paste0('t/sen/spec: ', round(apply(threshold_frame, 2, median), 2)), cex=0.8)
#   # }
#   
#   
# }
# 
# print(quantile(auc_vec))
# print(quantile(n_vec))
# print(threshold_frame)
# 
# print(quantile(accuracy_v))
# print(quantile(balanced_accuracy_v))
# print(quantile(f1))
# 
# print(quantile(threshold_frame[,2]))
# print(quantile(threshold_frame[,3]))
# 
# importance_matrix = xgb.importance(colnames(xgb_train), model = model)
# xgb.plot.importance(importance_matrix)
# importance_matrix <- importance_matrix[order(-importance_matrix$Importance), ]
# importance_matrix[1:20,]
# # 
# # timp <- t(imp_merge)
# # timp <- data.frame(timp)
# # colnames(timp) <- timp[1,]
# # timp = timp[-1,]
# # for (jj in c(1:ncol(timp))) {
# #   timp[, jj] <- as.numeric(timp[, jj])
# # }
# # 
# # medians <- data.frame(sapply(timp, function(x) median(as.numeric(x), na.rm = T) ))
# # colnames(medians) <- c('medians')
# # medians$name <- row.names(medians)
# # medians <- medians[order(-medians$medians), ]
# # medians
# 
# ## visualise trees
# require(DiagrammeR)
# xgb.plot.tree(model = model, trees = 1)
# 
# hist(outputTable_cat$pred, 100)
# 
# 
# ##################################
# # standard metrics
# 
# base_predictions_last_cbg <- function(ID, Glu, label, hypothresh) {
#   
#   # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
#   d <- data.table(ID, Glu, label)
#   d[, 'last_below_thresh' := ifelse(Glu[length(Glu)] < hypothresh, 1, 0) ,by=.(ID)]
#   d[, 'n' := c(1: .N), by=.(ID)]
#   single_d <- d[n==1]
#   # print(confusionMatrix(as.factor(single_d$last_below_thresh),
#                         # as.factor(single_d$label),
#                         # positive = '1',
#                         # prevalence = (sum(single_d$label) / nrow(single_d))))
#   
#   cm = confusionMatrix(as.factor(single_d$last_below_thresh),
#                        as.factor(single_d$label),
#                        positive = '1',
#                        prevalence = (sum(single_d$label) / nrow(single_d)))
#   
#   return(list(as.numeric(cm$overall[1]), # accuracy
#               as.numeric(cm$byClass[11]), # balanced accuracy
#               as.numeric(cm$byClass[7]),# f1
#               as.numeric(cm$byClass[1]),# sens
#               as.numeric(cm$byClass[2]))) # spec
#   
# }
# 
# base_predictions_prevHypo <- function(ID, Glu, label, hypothresh) {
#   
#   # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
#   d <- data.table(ID, Glu, label)
#   d[, 'prior_below_3' := ifelse(min(Glu) < hypothresh, 1, 0) ,by=.(ID)]
#   d[, 'n' := c(1: .N), by=.(ID)]
#   single_d <- d[n==1]
#   # print(confusionMatrix(as.factor(single_d$prior_below_3),
#   #                       as.factor(single_d$label),
#   #                       positive = '1',
#   #                       prevalence = (sum(single_d$label) / nrow(single_d))))
#   
#   cm = confusionMatrix(as.factor(single_d$prior_below_3),
#                        as.factor(single_d$label),
#                        positive = '1',
#                        prevalence = (sum(single_d$label) / nrow(single_d)))
#   
#   return(list(as.numeric(cm$overall[1]), # accuracy
#               as.numeric(cm$byClass[11]), # balanced accuracy
#               as.numeric(cm$byClass[7]),# f1
#               as.numeric(cm$byClass[1]),# sens
#               as.numeric(cm$byClass[2]))) # spec
#   
# }
# 
# base_predictions_cv <- function(ID, Glu, label) {
#   
#   # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
#   d <- data.table(ID, Glu, label)
#   
#   d[, 'sd_cbg' := sd(Glu) ,by=.(ID)]
#   d[, 'mean_cbg' := mean(Glu) ,by=.(ID)]
#   d[, 'cV' := sd_cbg/mean_cbg ,by=.(ID)]
#   
#   d[, 'n' := c(1: .N), by=.(ID)]
#   
#   single_d <- d[n==1]
#   
#   require(pROC)
#   roc_obj <- roc(single_d$label, single_d$cV)
#   print(auc(roc_obj))
#   
#   bestthresh <- as.numeric(coords(roc_obj, "best", ret="threshold", transpose = FALSE)[1,])
#   # print(confusionMatrix(as.factor(ifelse(single_d$cV > bestthresh, 1, 0)),
#   #                       as.factor(single_d$label),
#   #                       positive = '1',
#   #                       prevalence = (sum(single_d$label) / nrow(single_d))))
#   cm = confusionMatrix(as.factor(ifelse(single_d$cV > bestthresh, 1, 0)),
#                        as.factor(single_d$label),
#                        positive = '1',
#                        prevalence = (sum(single_d$label) / nrow(single_d)))
#   
#   return(list(as.numeric(auc(roc_obj)),
#               as.numeric(cm$overall[1]), # accuracy
#               as.numeric(cm$byClass[11]), # balanced accuracy
#               as.numeric(cm$byClass[7]),# f1
#               as.numeric(cm$byClass[1]),# sens
#               as.numeric(cm$byClass[2]))) # spec
#   
#   
# }
# 
# #########
# ## test by bootstrapping
# 
# single_n <- x[n==1]
# ids_for_split <- data.table('id' = single_n$ID,
#                             'label' = single_n$label)
# 
# bootst <- 200
# bootst_n <- 1000
# 
# # last CBG
# lasthV_acc <- rep(0, bootst)
# lasthV_bal_acc <- rep(0, bootst)
# lasthV_f1 <- rep(0, bootst)
# lasthV_sens <- rep(0, bootst)
# lasthV_spec <- rep(0, bootst)
# 
# # any hypo
# predhV_acc <- rep(0, bootst)
# predhV_bal_acc <- rep(0, bootst)
# predhV_f1 <- rep(0, bootst)
# predhV_sens <- rep(0, bootst)
# predhV_spec <- rep(0, bootst)
# # cV
# aucV <- rep(0, bootst)
# predcvV_acc <- rep(0, bootst)
# predcvV_bal_acc <- rep(0, bootst)
# predcvV_f1 <- rep(0, bootst)
# predcvV_sens <- rep(0, bootst)
# predcvV_spec <- rep(0, bootst)
# 
# for (j in c(1:bootst)) {
#   
#   sub <- sample(nrow(ids_for_split), bootst_n)
#   bootstrap_ids <- ids_for_split[sub]
#   
#   dat <- merge(bootstrap_ids, x, by.x = 'id', by.y = 'ID')
#   
#   lasth <- base_predictions_last_cbg(dat$id, dat$Glu, dat$label.x, 3)
#   predh <- base_predictions_prevHypo(dat$id, dat$Glu, dat$label.x, 3)
#   predcv<- base_predictions_cv(dat$id, dat$Glu, dat$label.x)
#   
#   lasthV_acc[j] <- lasth[[1]]
#   lasthV_bal_acc[j] <- lasth[[2]]
#   lasthV_f1[j] <- lasth[[3]]
#   lasthV_sens[j] <- lasth[[4]]
#   lasthV_spec[j] <- lasth[[5]]
#   
#   predhV_acc[j] <- predh[[1]]
#   predhV_bal_acc[j] <- predh[[2]]
#   predhV_f1[j] <- predh[[3]]
#   predhV_sens[j] <- predh[[4]]
#   predhV_spec[j] <- predh[[5]]
#   
#   aucV[j] <- predcv[[1]]
#   predcvV_acc[j] <- predcv[[2]]
#   predcvV_bal_acc[j] <- predcv[[3]]
#   predcvV_f1[j] <- predcv[[4]]
#   predcvV_sens[j] <- predcv[[5]]
#   predcvV_spec[j] <- predcv[[6]]
# }
# 
# 
# quantile(lasthV_acc)
# quantile(lasthV_bal_acc)
# quantile(lasthV_f1, na.rm = T)
# quantile(lasthV_sens)
# quantile(lasthV_spec)
# 
# quantile(predhV_acc)
# quantile(predhV_bal_acc)
# quantile(predhV_f1)
# quantile(predhV_sens)
# quantile(predhV_spec)
# 
# quantile(aucV)
# 
# quantile(predcvV_acc)
# quantile(predcvV_bal_acc)
# quantile(predcvV_f1)
# quantile(predcvV_sens)
# quantile(predcvV_spec)
# 
# 
