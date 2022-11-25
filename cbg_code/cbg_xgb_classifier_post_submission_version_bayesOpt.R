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
library(dplyr)

source('~/Documents/code/cbg_code/cbg_functions.R')

# install.packages(c('data.table', 'tidyverse', 'tidyquant', 'timetk', 'padr', 'zoo', 'imputeTS', 'xgboost', 'caret', 'pROC'))

dataset <- 'huge'
days_n  <- 7
minimum_n_cbgs <- days_n + 2
hypo_threshold <- 3
lockout = 7

ratio = 2

exclude_locations <- function(data, exclude_vector) {
  
  x_ex <- data
  
  for (e in c(1:length(exclude_vector))) {
    x_ex <- x_ex[-grep(exclude_vector[e], x_ex$loc, ignore.case = TRUE), ]
  }
  
  return(x_ex)
}

initial_process_label <- function(x, days_n, minimum_n_cbgs, hypo_threshold) {
  
  # add unique ID by 2 cols - uID and admission vec
  x <- x %>% mutate(ID = group_indices(x, .dots=c("uID", "admission_vec"))) 

  ## truncate all admissions to the right length
  x$day <- as.Date(x$dateTime)
  x[, 'correct_duration_days' := ifelse(day <= (min(day) + (days_n - 1)), 1, 0), by=.(ID)]
  x <- x[correct_duration_days == 1]
  
  # limit to those with minimum n CBGs
  x[, 'N_truncated' := .N, by=.(ID)]
  x <- x[N_truncated >= minimum_n_cbgs]
  
  # split out last day and add label
  x[, 'flag_last_day' := ifelse(day == max(day), 1, 0), by=.(ID)]
  x[, 'label' := ifelse(min(Glu[flag_last_day == 1]) <= hypo_threshold, 1, 0), by=.(ID)]
  
  # remove the last day from the training set
  x <- x[flag_last_day == 0]
  
  x[, 'N' := .N, by=.(ID)]
  x <- x[N>days_n]
  
  x[, 'n' := c(1 : .N), by=.(ID)]
  print(sum(x[n==1]$label))
  print(nrow(x[n==1]))
  
  return(x)
  
}

prior_process <- function(x, days_n, minimum_n_cbgs) {
  
  # add unique ID by 2 cols - uID and admission vec
  x <- x %>% mutate(ID = group_indices(x, .dots=c("uID", "prior_admission_vec"))) 
  
  ## truncate all admissions to the right length
  x$day <- as.Date(x$dateTime)
  x[, 'correct_duration_days' := ifelse(day <= (min(day) + (days_n - 1)), 1, 0), by=.(ID)]
  x <- x[correct_duration_days == 1]
  
  # limit to those with minimum n CBGs
  x[, 'N_truncated' := .N, by=.(ID)]
  x <- x[N_truncated >= minimum_n_cbgs]
  
  return(x)
  
}

case_control_ratio <- function(x, ratio) {
  
  # produce fixed ratio of case to no case
  event_ids <- unique(x[label==1]$ID)
  no_event_ids <- unique(x[label==0]$ID)
  id_sample <- no_event_ids[sample(length(no_event_ids), round(length(event_ids) * ratio), 0)]
  
  neg <- x[ID %in% id_sample]
  pos <- x[ID %in% event_ids]
  
  x <- rbind(neg, pos)
  
  print(sum(x[n==1]$label))
  print(nrow(x[n==1]))
  
  return(x)
  
}

# process by admission function
admission_process <- function(dt) {
  
  s <- dt
  
  s[, 'max_by_day' := max(Glu), by=.(ID, day)]
  s[, 'min_by_day' := min(Glu), by=.(ID, day)]
  s[, 'median_by_day' := median(Glu), by=.(ID, day)]
  s[, 'iqr_by_day' := quantile(Glu)[4]-quantile(Glu)[2], by=.(ID, day)]
  s[, 'sd_by_day' := sd(Glu) ,by=.(ID, day)]
  s[, 'mean_by_day' := mean(Glu) ,by=.(ID, day)]
  s[, 'cV_by_day' := sd_by_day/mean_by_day ,by=.(ID, day)]
  
  s[, 'training_min' := min(Glu) ,by=.(ID)]
  s[, 'training_max' := max(Glu) ,by=.(ID)]
  
  s[, 'day_n' := c(1:.N), by=.(ID, day)]
  s[, 'day_N' := .N, by=.(ID, day)]
  s[, 'gradient' := as.numeric(lm(Glu ~ dateTime)$coefficients[2]), by=.(ID)]
  s[, 'sd_cbg' := sd(Glu) ,by=.(ID)]
  s[, 'mean_cbg' := mean(Glu) ,by=.(ID)]
  s[, 'cV' := sd_cbg/mean_cbg ,by=.(ID)]
  
  s[, 'ID_n' := c(1:.N), by=.(ID)]
  s[, 'ID_N' := .N, by=.(ID)]
  
  return(s)
}

prior_admission_process <- function(dt) {
  
  s <- dt
  
  s[, 'prior_overall_min' := min(Glu) ,by=.(ID)]
  s[, 'prior_overall_max' := max(Glu) ,by=.(ID)]
  s[, 'prior_duration' := max(unix_dateTime) - min(unix_dateTime) ,by=.(ID)]
  s[, 'prior_gradient' := as.numeric(lm(Glu ~ dateTime)$coefficients[2]), by=.(ID)]
  s[, 'prior_sd_cbg' := sd(Glu) ,by=.(ID)]
  s[, 'prior_mean_cbg' := mean(Glu) ,by=.(ID)]
  s[, 'prior_cV' := prior_sd_cbg/prior_mean_cbg ,by=.(ID)]
  s[, 'prior_N_tests' := .N ,by=.(ID)]
  
  return(s)
  
}

x <- fread(paste0('~/Documents/data/CBGdata/', dataset, '_unipoc_time_series_cohort_first_', days_n,'_days.csv'))
x <- x[order(x$uID, x$admission_vec)]

x <- exclude_locations(x, c('dialysis', 'renal', 'abbott', 'Ketones', 'X Ray', 'X-ray',
                            'training', 'staff', 'returned', 'obsolete', 'old', 'Meterx',
                            'dead', 'BioChemistry', 'weldon', 'radiology', 'outpatient',
                            'rcas', 'Labour Ward', 'Neonatal Unit Top Floor'))

data.orig <- x

x <- initial_process_label(x, days_n, minimum_n_cbgs, hypo_threshold)

x <- case_control_ratio(x, ratio)

s <- admission_process(x)

# get data from most recent admission, and all prior admission N
s[, 'qualifying_admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(uID)]
s <- s %>% select(-admission_vector, -n_admissions, -admission_vec)

# identify the most recent admission - to ensure each ID appears only once
s[, 'last_admission_flag' := ifelse(qualifying_admission_vec == max(qualifying_admission_vec), 1, 0), by=.(uID)]

    # prior data from ids with multiple admissions
    ma <- s[last_admission_flag==0]
    # ids from prior admissions:
    ids_old_admissions <- unique(ma$uID)
    # first CBGs from all last (index) admissions
    index_dates <- s[last_admission_flag == 1 &ID_n == 1]
    index_dates <- index_dates %>% select(uID, unix_dateTime)
    index_dates <- index_dates %>% rename(index_first_unix = unix_dateTime)
    
    # merge with all CBG data
    # to get all prior CBGs for the uIDs of interest
    all_prior_cbgs <- merge(index_dates, data.orig, by.x = 'uID', by.y = 'uID')
    all_prior_cbgs <- all_prior_cbgs[unix_dateTime < index_first_unix]
    
    # process all prior admissions from the prior CBG dataset
    all_prior_cbgs[, 'prior_admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(uID)]
    all_prior_cbgs <- prior_process(all_prior_cbgs, days_n, minimum_n_cbgs) 
    
    all_prior_cbgs <- prior_admission_process(all_prior_cbgs)
    
    # add time from index
    all_prior_cbgs[, 'time_from_index' := index_first_unix - max(unix_dateTime), by=.(uID)]
    all_prior_cbgs[, 'last_prior_admission_flag' := ifelse(prior_admission_vec == max(prior_admission_vec), 1, 0), by=.(uID)]
    
    # only take the last admission
    last_prior_admission <- all_prior_cbgs[last_prior_admission_flag == 1]
    last_prior_admission[, 'single_row_n' := c(1:.N), by=.(uID)]
    # single row
    last_prior_admission <- last_prior_admission[single_row_n==1]
    
    # file for merge back to index admission dataset
    last_prior_admission <- last_prior_admission %>% select(uID, prior_overall_min, prior_overall_max, prior_duration, prior_gradient, prior_sd_cbg, prior_mean_cbg,  prior_cV, prior_N_tests, time_from_index)

# merge prior back to index dataset
s_merge <- merge(s, last_prior_admission, by.x = 'uID', by.y = 'uID', all.x = T)

## lastly
## ensure that only data from the last admission goes forward into training
s <- s_merge
s <- s[last_admission_flag == 1]

## per ID generate wide data
ids <-unique(s$ID)
for (i in c(1:length(ids))) {
  
  if (i %% 100 == 0) {print(i/length(ids))}
  
  sub <- s[ID == ids[i]]
  sub <- sub[day_n == 1]
  
  sub <- sub %>% select(ID, dateTime, day, Glu,
                        max_by_day, min_by_day, median_by_day, day_N, iqr_by_day, cV_by_day,
                        training_min, training_max, cV, gradient,
                        prior_overall_min, prior_overall_max, prior_duration, prior_gradient,
                        prior_sd_cbg, prior_mean_cbg, prior_cV, prior_N_tests, time_from_index,
                        label)
  
  if ((nrow(sub) > 0) & (as.numeric(min(diff(as.Date(sub$dateTime)))) <= 1)) { # second arguement needed to ensure that minimum time step is not greater than 1 day (thicken will no work if it is)
    
    #sub$day_numeric <- c(1:nrow(sub))
    
    if ((as.Date(sub$dateTime[nrow(sub)]) - as.Date(sub$dateTime[1]) < (days_n - 1))) {
      sub <- as.data.frame(sub)
      sub[nrow(sub)+1,] <- NA
      sub <- as.data.table(sub)
      
      sub$dateTime[nrow(sub)] <- (sub$dateTime[1] + ((days_n - 1) * 60*60*24))
      
    }
    
    cV <- sub$cV[1]
    gradient <- sub$gradient[1]
    training_min <- sub$training_min[1]
    training_max <- sub$training_max[1]
    id <- sub$ID[1]
    prior_overall_min <- sub$prior_overall_min[1]
    prior_overall_max <- sub$prior_overall_max[1]
    prior_duration <- sub$prior_duration[1]
    prior_gradient <- sub$prior_gradient[1]
    prior_sd_cbg <- sub$prior_sd_cbg[1]
    prior_mean_cbg <- sub$prior_mean_cbg[1]
    prior_cV <- sub$prior_cV[1]
    prior_N_tests <- sub$prior_N_tests[1]
    time_from_index <- sub$time_from_index[1]
    label <- sub$label[1]
    
    sub <- sub %>% select(dateTime, Glu, max_by_day, min_by_day, median_by_day, day_N, iqr_by_day, cV_by_day)
    sub <- sub %>% thicken("1 day") %>% select(-dateTime) %>% pad()
    
    # sub.locf <- na.locf(sub)
    sub.locf <- sub
    sub.locf$Glu <- NULL
    
    require(reshape2)
    melt.sub <- melt(sub.locf, id.vars='dateTime_day')
    m <- as.data.frame(t(melt.sub))
    m$cV <- cV
    m$gradient <- gradient
    m$training_min <- training_min
    m$training_max <- training_max
    m$id <- id
    m$prior_overall_min <- prior_overall_min
    m$prior_overall_max <- prior_overall_max
    m$prior_duration <- prior_duration
    m$prior_gradient <- prior_gradient
    m$prior_sd_cbg <- prior_sd_cbg
    m$prior_mean_cbg <- prior_mean_cbg
    m$prior_cV <- prior_cV
    m$prior_N_tests <- prior_N_tests
    m$time_from_index <- time_from_index
    m$label <- label
    m <- m[-1, ]
    colnames(m) <- c(m[1,1:(ncol(m) - 15)],'cV' ,'gradient' , 'training_min', 'training_max','id',
                     'prior_overall_min', 'prior_overall_max', 'prior_duration', 'prior_gradient',
                     'prior_sd_cbg', 'prior_mean_cbg', 'prior_cV', 'prior_N_tests', 'time_from_index',
                     'label')
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
cVday_names <- coln[coln == 'cV_by_day']

other_names <- coln[coln != 'max_by_day' & coln != 'min_by_day' & coln != 'iqr_by_day' &
                      coln != 'median_by_day' & coln != 'day_N'  & coln != 'cV_by_day']

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
for (c in c(1:length(cVday_names))) {
  cVday_names[c] <- paste0(cVday_names[c], '_', c)
}

coln_n <- c(max_names, min_names, med_names, N_names, iqr_names, cVday_names, other_names)
colnames(exp) <- coln_n

export <- exp

for (j in c(1:ncol(export))) {
  export[, j] <- as.numeric(export[, j])
}

write.table(export, file = paste0('Documents/data/CBGdata/abstract_exports/with_prior_export_admissionDuration_', days_n, '_days_hypothresh_NAs_included_', hypo_threshold, '_ratio_', ratio,'.csv'), sep = ',', row.names = F)

# export1 <- fread(paste0('Documents/data/CBGdata/abstract_exports/export_admissionDuration_', days_n, '_days_hypothresh_NAs_included_', hypo_threshold, '_ratio_', ratio, '.csv'))

##### build model - bayes opt
require(ParBayesianOptimization)

splits <- rsample::initial_split(export, prop = 0.7)
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
  select(!id, !label) %>%
  as.matrix()

# Get the target variable
y <- train_df %>%
  pull(label)

# Cross validation folds
folds <- list(fold1 = as.integer(seq(1, nrow(X), by = 5)),
              fold2 = as.integer(seq(2, nrow(X), by = 5)),
              fold3 = as.integer(seq(3, nrow(X), by = 5)),
              fold4 = as.integer(seq(4, nrow(X), by = 5)),
              fold5 = as.integer(seq(5, nrow(X), by = 5))
              )

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
                  nround = 600,
                  folds = folds,
                  prediction = TRUE,
                  early_stopping_rounds = 50,
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
               max_depth = c(1L, 30L),
               min_child_weight = c(1, 100),
               subsample = c(0.01, 1),
               lambda = c(1, 100),
               alpha = c(1, 100))
set.seed(42)
bayes_out <- bayesOpt(FUN = obj_func,
                      bounds = bounds,
                      initPoints = length(bounds) + 2,
                      iters.n = 60,
                      iters.k = 6,
                      plotProgress = T)

bayes_out$scoreSummary
data.frame(getBestPars(bayes_out))

best_params = data.frame(getBestPars(bayes_out))
## need to save out bayes score summary and best_params for later use / documentation
## also save out plot

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## leave one out

export.orig = export
sample_rows = 2000
export = export.orig[sample(nrow(export.orig), sample_rows), ]
export$prediction = NULL

export = data.table(export)
round_n = 100

use_bayes_opt = 1

k = nrow(export)

prediction <- rep(0, k)
ground_t <- rep(0, k)

n_vec <- rep(0, k)

for (kfold in c(1:nrow(export))) {
  
  print(kfold/k)
  
  train = export[-kfold,]
  test  = export[kfold,]
  
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
  if (use_bayes_opt == 1) {
    param <- list(max.depth = best_params$max_depth,
                  eta = best_params$eta,
                  min_child_weight = best_params$min_child_weight,
                  subsample = best_params$subsample,
                  lambda= best_params$lambda,
                  alpha = best_params$alpha
                  )
  } else {
    param <- list(max.depth = 2, eta = 0.1, nthread = 64, min_child_weight = 1)
  }
  
  model = xgb.train(param, data = xgb_train, watchlist=watchlist, nrounds = round_n, verbose = 0)
  
  # model = xgb.train(data = xgb_train, max.depth = maxd, watchlist=watchlist, nrounds = 100, nthread = 16)
  # plot(model$evaluation_log$test_rmse)
  
  n = which(model$evaluation_log$test_rmse == min(model$evaluation_log$test_rmse))
  print(n)
  n_vec[kfold] <- n[1]
  
  final = xgboost(data = xgb_train, nrounds = n[1], verbose = 0)
  
  pred_y = predict(final, xgb_test)
  
  prediction[kfold] = pred_y
  ground_t[kfold] = test_y$label
  # hist(pred_y, 100)
  
}

whole_data <-  data.matrix(export[, -c('id', 'label')])
whole_y    <- export[,'label']

whole_train = xgb.DMatrix(data = as.matrix(whole_data), label = whole_y$label)


whole_model <- xgb.train(param,
                         data = whole_train,
                         watchlist=watchlist,
                         nrounds = quantile(n_vec)[3],
                         verbose = 0)


    ## save out last trained model
    model_save_name = paste0(paste0('Documents/data/CBGdata/savedModels/wholeData_modelSAVE_', days_n,
                                    '_days_hypothresh_', hypo_threshold,
                                    '_ratio_', ratio,
                                    '_maxd_', best_params$max_depth,'.model'))
    model_dump_name = paste0(paste0('Documents/data/CBGdata/savedModels/wholeData_modelDUMP_', days_n,
                                    '_days_hypothresh_', hypo_threshold,
                                    '_ratio_', ratio,
                                    '_maxd_', best_params$max_depth,'.txt'))
    
    xgb.save(whole_model, model_save_name)
    xgb.dump(whole_model, model_dump_name, with_stats = FALSE)

print(quantile(n_vec))
hist(n_vec, 80)
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

hist(prediction, 60)

## generate validation accuracy metric
pred_acc <- ifelse(prediction >= c$threshold, 1, 0)
cm <- confusionMatrix(as.factor(pred_acc), as.factor(ground_t), positive = '1')
print(threshold)
print(cm)

print(cm$byClass[7])

importance_matrix = xgb.importance(colnames(xgb_train), model = model)
xgb.plot.importance(importance_matrix)
importance_matrix <- importance_matrix[order(-importance_matrix$Importance), ]
importance_matrix[1:20,]

## save out prediction data frame
export$prediction = prediction

write.table(export, file = paste0('Documents/data/CBGdata/abstract_exports/export_admissionDuration_', days_n, '_days_hypothresh_NAs_included_', hypo_threshold, '_ratio_', ratio,'_WITH_LOO_PRED.csv'), sep = ',', row.names = F)


plot(export$cV, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))

plot(export$min_by_day_13, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))
plot(export$min_by_day_6, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))
plot(export$min_by_day_1, export$prediction, cex = 0.4, col=rgb(0,0,0,0.4))

boxplot(export$prediction ~ cut(export$min_by_day_6, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$min_by_day_4, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$min_by_day_1, 60), varwidth=T)

boxplot(export$prediction ~ cut(export$time_from_index, 60), varwidth=T)

boxplot(export$prediction ~ cut(export$prior_overall_min, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$prior_overall_max, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$prior_duration, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$prior_gradient, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$prior_cV, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$prior_N_tests, 60), varwidth=T)

boxplot(export$prediction ~ cut(export$max_by_day_6, 20), varwidth=T)
boxplot(export$prediction ~ cut(export$max_by_day_3, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$max_by_day_1, 60), varwidth=T)

boxplot(export$prediction ~ cut(export$iqr_by_day_6, 20), varwidth=T)
boxplot(export$prediction ~ cut(export$iqr_by_day_3, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$iqr_by_day_1, 60), varwidth=T)

boxplot(export$prediction ~ export$day_N_6, varwidth=T)
boxplot(export$prediction ~ export$day_N_3, varwidth=T)
boxplot(export$prediction ~ export$day_N_1, varwidth=T)

boxplot(export$prediction ~ cut(export$gradient, 60), varwidth=T)
boxplot(export$prediction ~ cut(export$cV, 60), varwidth=T)


## visualise trees
require(DiagrammeR)
jpeg(paste0('Documents/data/CBGdata/treeplots/single_tree_plot_', days_n, '_days_hypothresh_NAs_included_', hypo_threshold, '_ratio_', ratio, '_maxd_', maxd,'.jpeg'))
xgb.plot.tree(model = model, trees = c(1), render = TRUE)
dev.off()

hist(outputTable_cat$pred, 100)


##################################
# standard metrics

base_predictions_last_cbg <- function(ID, Glu, label, hypothresh) {
  
  # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
  d <- data.table(ID, Glu, label)
  d[, 'last_below_thresh' := ifelse(Glu[length(Glu)] < hypothresh, 1, 0) ,by=.(ID)]
  d[, 'n' := c(1: .N), by=.(ID)]
  single_d <- d[n==1]
  # print(confusionMatrix(as.factor(single_d$last_below_thresh),
                        # as.factor(single_d$label),
                        # positive = '1',
                        # prevalence = (sum(single_d$label) / nrow(single_d))))
  
  cm = confusionMatrix(as.factor(single_d$last_below_thresh),
                       as.factor(single_d$label),
                       positive = '1',
                       prevalence = (sum(single_d$label) / nrow(single_d)))
  
  return(list(as.numeric(cm$overall[1]), # accuracy
              as.numeric(cm$byClass[11]), # balanced accuracy
              as.numeric(cm$byClass[7]),# f1
              as.numeric(cm$byClass[1]),# sens
              as.numeric(cm$byClass[2]))) # spec
  
}

base_predictions_prevHypo <- function(ID, Glu, label, hypothresh) {
  
  # ID = dat$id; Glu = dat$Glu; label = dat$label.y; hypothresh = 3
  d <- data.table(ID, Glu, label)
  d[, 'prior_below_3' := ifelse(min(Glu) < hypothresh, 1, 0) ,by=.(ID)]
  d[, 'n' := c(1: .N), by=.(ID)]
  single_d <- d[n==1]
  # print(confusionMatrix(as.factor(single_d$prior_below_3),
  #                       as.factor(single_d$label),
  #                       positive = '1',
  #                       prevalence = (sum(single_d$label) / nrow(single_d))))
  
  cm = confusionMatrix(as.factor(single_d$prior_below_3),
                       as.factor(single_d$label),
                       positive = '1',
                       prevalence = (sum(single_d$label) / nrow(single_d)))
  
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
  # print(confusionMatrix(as.factor(ifelse(single_d$cV > bestthresh, 1, 0)),
  #                       as.factor(single_d$label),
  #                       positive = '1',
  #                       prevalence = (sum(single_d$label) / nrow(single_d))))
  cm = confusionMatrix(as.factor(ifelse(single_d$cV > bestthresh, 1, 0)),
                       as.factor(single_d$label),
                       positive = '1',
                       prevalence = (sum(single_d$label) / nrow(single_d)))
  
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

bootst <- 200
bootst_n <- 1000

# last CBG
lasthV_acc <- rep(0, bootst)
lasthV_bal_acc <- rep(0, bootst)
lasthV_f1 <- rep(0, bootst)
lasthV_sens <- rep(0, bootst)
lasthV_spec <- rep(0, bootst)

# any hypo
predhV_acc <- rep(0, bootst)
predhV_bal_acc <- rep(0, bootst)
predhV_f1 <- rep(0, bootst)
predhV_sens <- rep(0, bootst)
predhV_spec <- rep(0, bootst)
# cV
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
  
  lasth <- base_predictions_last_cbg(dat$id, dat$Glu, dat$label.x, 3)
  predh <- base_predictions_prevHypo(dat$id, dat$Glu, dat$label.x, 3)
  predcv<- base_predictions_cv(dat$id, dat$Glu, dat$label.x)
  
  lasthV_acc[j] <- lasth[[1]]
  lasthV_bal_acc[j] <- lasth[[2]]
  lasthV_f1[j] <- lasth[[3]]
  lasthV_sens[j] <- lasth[[4]]
  lasthV_spec[j] <- lasth[[5]]
  
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


quantile(lasthV_acc)
quantile(lasthV_bal_acc)
quantile(lasthV_f1, na.rm = T)
quantile(lasthV_sens)
quantile(lasthV_spec)

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


