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
library(Metrics)


# 101355572

dataset <- 'huge'
days_n  <- 3
minimum_n_cbgs <- days_n + 2
# hypo_threshold <- 1.9

x <- fread(paste0('~/Documents/data/CBGdata/', dataset, '_unipoc_time_series_cohort_first_', days_n,'_days.csv'))
x <- x[order(x$uID, x$admission_vec)]

exclude <- c('dialysis', 'renal', 'abbott', 'Ketones', 'X Ray', 'X-ray', 'training', 'staff', 'returned', 'obsolete', 'old', 'Meterx', 'dead', 'BioChemistry', 'weldon', 'radiology', 'outpatient', 'rcas', 'Labour Ward')
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
x[, 'last_day_min' := min(Glu[flag_last_day == 1]), by=.(ID)]
x[, 'last_day_max' := max(Glu[flag_last_day == 1]), by=.(ID)]
# x[, 'label' := ifelse(min(Glu[flag_last_day == 1]) <= hypo_threshold, 1, 0), by=.(ID)]
# remove the last day from the training set
x <- x[flag_last_day == 0]

x[, 'N' := .N, by=.(ID)]
x <- x[N>days_n]

x[, 'n' := c(1 : .N), by=.(ID)]
# print(sum(x[n==1]$label))
print(nrow(x[n==1]))

# plot(x$dateTime, x$Glu, pch = 16, cex = 0.6, col = ifelse(x$label == 0,
#                                                           rgb(0,0,0,0.4, maxColorValue = 1),
#                                                           rgb(1,0,0,0.4, maxColorValue = 1)))

# subsample the full dataset to minimise time
idsx <- unique(x$ID)
sample_size <- length(idsx)
sample_ids <- idsx[sample(length(idsx), sample_size)]

s <- x[ID %in% sample_ids]

# produce fixed ratio of case to no case
# ratio = 2
# event_ids <- unique(x[label==1]$ID)
# no_event_ids <- unique(x[label==0]$ID)
# id_sample <- no_event_ids[sample(length(no_event_ids), round(length(event_ids) * ratio), 0)]
# 
# neg <- x[ID %in% id_sample]
# pos <- x[ID %in% event_ids]
# 
# x <- rbind(neg, pos)
# 
# print(sum(x[n==1]$label))
# print(nrow(x[n==1]))

# create value at 0000 first day (same as first value), and 23.59 (same as last recorded value) on the last day to allow minute fills of data

returnFirstDayTime <- function(sub) {
  first_day <- substr(sub$dateTime[1], 1, 10)
  first_val <- sub$Glu[1]
  
  first_point_time <- paste0(first_day, ' 0:00:01')
  first_point_time <- as.POSIXct(first_point_time, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  
  sub <- data.frame(sub$ID, sub$dateTime, sub$Glu, sub$loc)
  colnames(sub) <- c('ID', 'dateTime', 'Glu', 'loc')
  
  insert_top <- as.data.frame(matrix(nrow = 1, ncol = 4))
  colnames(insert_top) <- c('ID', 'dateTime', 'Glu', 'loc')
  insert_top$ID <- sub$ID[1]
  insert_top$dateTime <- first_point_time
  insert_top$Glu <- first_val
  insert_top$location <- sub$location[1]
  #insert_top$op <- sub$op[1]
  
  return(insert_top)
  
}

returnLastDayTime <- function(sub, days_n) {
  # last_day <- substr(sub$dateTime[nrow(sub)], 1, 10)
  last_day <- as.Date(substr(sub$dateTime[1], 1, 10)) + (days_n - 2)
  last_val <- sub$Glu[nrow(sub)]
  
  last_point_time <- paste0(last_day, ' 23:59:59')
  last_point_time <- as.POSIXct(last_point_time, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  
  sub <- data.frame(sub$ID, sub$dateTime, sub$Glu, sub$loc)
  colnames(sub) <- c('ID', 'dateTime', 'Glu', 'loc')
  
  insert_end <- as.data.frame(matrix(nrow = 1, ncol = 4))
  colnames(insert_end) <- c('ID', 'dateTime', 'Glu', 'loc')
  insert_end$ID <- sub$ID[1]
  insert_end$dateTime <- last_point_time
  insert_end$Glu <- last_val
  insert_end$location <- sub$location[nrow(sub)]
  #insert_end$op <- sub$op[nrow(sub)]
  
  return(insert_end)
  
}

densityMap <- function(v1, s, i, b, label) {
  
  #v1 = sub$ID[1]; s = out$Glu; i = 1000
  s <- diff(s)
  
  start_point <- 1
  offset <- i
  
  initial_point <- start_point + (2 * offset)
  end_padding <- length(s) - initial_point
  
  x <- s[start_point:end_padding]
  y <- s[(start_point + offset):(end_padding+offset)]
  z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
  
  # from paper 2D projection 
  u <- 1/3 * (x + y + z)
  v <- (1/sqrt(6)) * (x + y - (2 * z))
  w <- (1/sqrt(2)) * (x - y)
  
  # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  # 
  # v <- range01(v)
  # w <- range01(w)
  
  df <- data.frame(v, w)
  
  pixel_n = 200
  
  if (label == 1) {
    jpeg(paste0('~/Documents/data/plots/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    p <- ggplot(df, aes(x = v, y = w)) +
      stat_density2d(aes(fill = ..density..), geom = 'tile', n = pixel_n, contour = F) +
      scale_fill_viridis() + 
      p <- p + theme_void() + theme(legend.position = "none")
    print(p)
    dev.off()
  } else {
    jpeg(paste0('~/Documents/data/plots/nil.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    p <- ggplot(df, aes(x = v, y = w)) +
      stat_density2d(aes(fill = ..density..), geom = 'tile', n = pixel_n, contour = F) +
      scale_fill_viridis()
    p <- p + theme_void() + theme(legend.position = "none")
    print(p)
    dev.off()
  }
  
}

densityMap_diff <- function(v1, s, i, b, label) {
  
  #v1 = sub$ID[1]; s = out$Glu; i = 480
  s <- diff(s)
  
  start_point <- 1
  offset <- i
  
  initial_point <- start_point + (2 * offset)
  end_padding <- length(s) - initial_point
  
  x <- s[start_point:end_padding]
  y <- s[(start_point + offset):(end_padding+offset)]
  z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
  
  # library(rgl); plot3d(x,y,z)
  
  # from paper 2D projection 
  u <- 1/3 * (x + y + z)
  v <- (1/sqrt(6)) * (x + y - (2 * z))
  w <- (1/sqrt(2)) * (x - y)
  
  # plot(v, w, cex = 0.4, pch=16)
  
  # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  # 
  # v <- range01(v)
  # w <- range01(w)
  
  df <- data.frame(v, w)
  
  # pixel_n = 200
  
  pl_min = -0.1
  pl_max = 0.1
  
  if (label == 1) {
    jpeg(paste0('~/Documents/data/plots/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point() + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  } else {
    jpeg(paste0('~/Documents/data/plots/nil.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point() + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  }
  
}

densityMap_diff_plus <- function(v1, s, i, b, label, limit_min, limit_max) {
  
  #v1 = sub$ID[1]; s = out$Glu; i = 480
  s <- s
  
  jpeg(paste0('~/Documents/data/plots_simple/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
  plot(s)
  dev.off()
  
  
  sd <- c(0, diff(s))
  
  s <- s * sd
  
  start_point <- 1
  offset <- i
  
  initial_point <- start_point + (2 * offset)
  end_padding <- length(s) - initial_point
  
  x <- s[start_point:end_padding]
  y <- s[(start_point + offset):(end_padding+offset)]
  z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
  
  # library(rgl); plot3d(x,y,z)
  
  # from paper 2D projection 
  u <- 1/3 * (x + y + z)
  v <- (1/sqrt(6)) * (x + y - (2 * z))
  w <- (1/sqrt(2)) * (x - y)
  
  # plot(v, w, cex = 0.4, pch=16)
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  v <- range01(v)
  w <- range01(w)
  
  df <- data.frame(v, w)
  
  pixel_n = 200
  
  pl_min = limit_min
  pl_max = limit_max
  
  if (label == 1) {
    jpeg(paste0('~/Documents/data/plots/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point(alpha = 4/10) + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  } else {
    jpeg(paste0('~/Documents/data/plots/nil.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point(alpha = 4/10) + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  }
  
  
}

# idVec <- unique(x$ID)
# print(length(idVec))
# #format_out <- as.data.frame(matrix(0, nrow = length(idVec), ncol = ((days_n - 1) * (60*24)+1)))
# for (j in c(1:length(idVec))) {
# #for (j in c(1:100)) {
#   
#   if (j%%100 == 0) {print(j)}
#   
#   id <- idVec[j]
#   sub <- x[ID == id]
#   label <- sub$label[1]
#   
#   insert_top <- returnFirstDayTime(sub)
#   insert_end <- returnLastDayTime(sub, days_n)
#   
#   # truncate to maximum 30 day duration of data
#   max_day <- as.Date(insert_top$dateTime) + 31
#   sub <- sub[as.Date(sub$dateTime) <= max_day]
#   
#   sub <- data.frame(sub$ID, sub$dateTime, sub$Glu, sub$loc)
#   colnames(sub) <- c('ID', 'dateTime', 'Glu', 'loc')
#   
#   # deal with instances where all cbgs the same
#   if (prod(diff(sub$Glu)) == 0) {
#     sub$Glu[1] <- jitter(sub$Glu[1])
#   } 
#   
#   out <- rbind(insert_top, sub, insert_end)
#   
#   out <- out %>% thicken("1 min") %>% select(-dateTime) %>% pad()
#   
#   # out$Glu <- na.spline(out$Glu)
#   out$Glu <- imputeTS::na_locf(out$Glu)
#   
#   # print(nrow(out))
#   
#   #out <- fill(out, c(V1, Glu, location, op))
#   
#   # pass to density map function
#   # densityMap_diff(out$ID[1], jitter(out$Glu), 480, 100, label)
#   # densityMap_diff_plus(out$ID[1], jitter(out$Glu), 240, 100, label, 0, 1)
#   # densityMap_diff_plus(out$ID[1], out$Glu, 240, 100, label, 0, 1)
#   
#   report_col <- c(out$ID[1], label, out$Glu)
#   if (j == 1) {
#     export = as.data.frame(report_col)
#   } else {
#     export = cbind(export, report_col)
#   }
# }
# 
# final_frame = t(export)

## generate data to produce per day features for xgb model

dimGrad <- function(day, glu, ID) {
  
  #print(ID)
  
  # day = ss$unix_dateTime
  # glu = ss$Glu
  
  #dt <- data.table(day, glu)

  n_measurements <- length(day)
  grad_vec  <- rep(0, n_measurements)
  inter_vec <- rep(0, n_measurements)
  for (j in c(1:n_measurements)) {
    grad_vec[j] <- lm(glu[j:n_measurements] ~ day[j:n_measurements])$coefficients[2]
    inter_vec[j] <- lm(glu[j:n_measurements] ~ day[j:n_measurements])$coefficients[1]
  }
  
  grad_vec[is.na(grad_vec)] <- 0
  inter_vec[is.na(inter_vec)] <- 0
  # grad_vec <- grad_vec[1:length(grad_vec) - 1]
  #   g = as.data.frame(grad_vec)
  # inter_vec <- inter_vec[1:length(inter_vec) - 1]
  #   inter = t(as.data.frame(inter_vec))
  
  # g <- g[rep(seq_len(nrow(g)), each = ncol(g)), ]# dplyr package
  # inter <- inter[rep(seq_len(nrow(inter)), each = ncol(inter)), ]# dplyr package
  
  out <- list(grad_vec, inter_vec)
  return(out)
  
}

# s <- x

s[, 'max_by_day' := max(Glu), by=.(ID, day)]
s[, 'min_by_day' := min(Glu), by=.(ID, day)]
s[, 'median_by_day' := median(Glu), by=.(ID, day)]
s[, 'iqr_by_day' := quantile(Glu)[4]-quantile(Glu)[2], by=.(ID, day)]

#s[, c('gradient_vector', 'intercept_vector') := dimGrad(unix_dateTime, Glu, ID), by=.(ID)]

s[, 'day_n' := c(1:.N), by=.(ID, day)]
s[, 'day_N' := .N, by=.(ID, day)]
s[, 'gradient' := as.numeric(lm(Glu ~ dateTime)$coefficients[2]), by=.(ID)]
s[, 'sd_cbg' := sd(Glu) ,by=.(ID)]
s[, 'mean_cbg' := mean(Glu) ,by=.(ID)]
s[, 'cV' := sd_cbg/mean_cbg ,by=.(ID)]

s[, 'ID_n' := c(1:.N), by=.(ID)]
s[, 'ID_N' := .N, by=.(ID)]

## demog from chi
year <- ifelse(nchar(s$uID == 9), substr(s$uID, 5, 6), substr(6, 7))
currentYr <- as.numeric(substr(s$datetime, 9, 10)) + 100
s$age <- currentYr - as.numeric(year) # representation of age

sex_digit <- ifelse(nchar(s$uID == 9), substr(s$uID, 8, 8), substr(9, 9))
s$sex <- ifelse(as.numeric(sex_digit) %% 2 == 0, 1, 0) # female = 1

hist(s[ID_n==1]$cV, 1000)
hist(s[ID_n==1]$last_day_min, 1000)


s[, 'sum_diff_day' := as.numeric(sum(diff(day))), by=.(ID)]

## per ID generate wide data
ids <-unique(s$ID)
for (i in c(1:length(ids))) {
  
  if (i %% 100 == 0) {print(i/length(ids))}
  
  sub <- s[ID == ids[i]]
  sub <- sub[day_n == 1]
  
  sub <- sub %>% select(ID, age, sex, loc, dateTime, day, Glu, max_by_day, min_by_day, median_by_day, iqr_by_day, day_N, cV, gradient, last_day_min)
  
  if ((nrow(sub) > 2) & (as.numeric(min(diff(as.Date(sub$dateTime)))) <= 1)) { # second arguement needed to ensure that minimum time step is not greater than 1 day (thicken will no work if it is)
    
    #sub$day_numeric <- c(1:nrow(sub))
    
    if ((as.Date(sub$dateTime[nrow(sub)]) - as.Date(sub$dateTime[1]) < (days_n - 1))) {
      sub <- as.data.frame(sub)
      sub[nrow(sub)+1,] <- NA
      sub <- as.data.table(sub)
      
      sub$dateTime[nrow(sub)] <- (sub$dateTime[1] + ((days_n - 1) * 60*60*24))
      
    }
    
    age <- sub$age[1]
    sex <- sub$sex[1]
    cV <- sub$cV[1]
    gradient <- sub$gradient[1]
    id <- sub$ID[1]
    last_day_min <- sub$last_day_min[1]
    
    sub <- sub %>% select(dateTime, Glu, max_by_day, min_by_day, median_by_day, day_N, iqr_by_day, age, sex, loc)
    sub <- sub %>% thicken("1 day") %>% select(-dateTime) %>% pad()
    
    imputeFlag <- rep(0, nrow(sub))
    imputeFlag[is.na(sub$Glu)] <- 1
    
    sub$imputeFlag <- imputeFlag
    
    sub.locf = sub
    #sub.locf <- na.locf(sub)
    sub.locf$Glu <- NULL
    
    require(reshape2)
    melt.sub <- melt(sub.locf, id.vars='dateTime_day')
    m <- as.data.frame(t(melt.sub))
    m$age <- age
    m$sex <- sex
    m$cV <- cV
    m$gradient <- gradient
    m$id <- id
    m$last_day_min <- last_day_min
    m <- m[-1, ]
    colnames(m) <- c(m[1,1:(ncol(m) - 6)],'age','sex','cV' ,'gradient' ,'id', 'last_day_min')
    m <- m[-1, ]
    
    if (i == 1) {
      export <- m
    } else {
      export <- rbind(export, m)
    }
    
  }
  
  
}

for (j in c(1:ncol(export))) {
  export[, j] <- as.numeric(export[, j])
}

write.table(export, file = paste0('Documents/data/CBGdata/abstract_exports/export_REGRESSOR_admissionDuration_', days_n, '_days_sampleSize_', sample_size, '.csv'), sep = ',', row.names = F)

export1 <- fread(paste0('Documents/data/CBGdata/abstract_exports/export_REGRESSOR_admissionDuration_', days_n, '_days_sampleSize_', sample_size, '.csv'))

##### build model
##
export = data.table(export)

# take log of target
export$last_day_min <- log(export$last_day_min)

k = 4
flds <- createFolds(export$last_day_min, k, list = T, returnTrain = F)

rmse_vec <- rep(0, k)
n_vec <- rep(0, k)
threshold_frame <- as.data.frame(matrix(nrow = k, ncol = 3))

accuracy_v <- rep(0, k)
balanced_accuracy_v <- rep(0, k)
f1 <- rep(0, k)

maxd = 4

for (kfold in c(1:k)) {
  
  train = export[-flds[[kfold]], ]
  test  = export[flds[[kfold]], ]
  
  #define predictor and response variables in training set
  train_x = data.matrix(train[, -c('id', 'last_day_min')])
  train_y = train[,'last_day_min']
  
  #define predictor and response variables in testing set
  test_x = data.matrix(test[, -c('id', 'last_day_min')])
  test_y = test[,'last_day_min']
  
  #define final training and testing sets
  xgb_train = xgb.DMatrix(data = as.matrix(train_x), label = train_y$last_day_min)
  xgb_test = xgb.DMatrix(data = as.matrix(test_x), label = test_y$last_day_min)
  
  #define watchlist
  watchlist = list(train=xgb_train, test=xgb_test)
  
  #fit XGBoost model and display training and testing data at each round
  param <- list(max.depth = maxd, eta = 0.01, nthread = 64, objective = "reg:pseudohubererror")
  model = xgb.train(param, data = xgb_train, watchlist=watchlist, nrounds = 1000, verbose = 0)
  # plot(model$evaluation_log$test_rmse)
  
  n = which(model$evaluation_log$test_mphe[-1] == min(model$evaluation_log$test_mphe[-1]))
  print(n)
  n_vec[kfold] <- n
  
  final = xgboost(data = xgb_train, max.depth = maxd, nrounds = n, verbose = 0)
  
  pred_y = predict(final, xgb_test)
  # hist(pred_y, 100)
  
  outputTable <- data.table('pred' = pred_y, 'ref' = test_y$last_day_min)
  
  if (kfold == 1) {
    outputTable_cat <- outputTable
  } else {
    outputTable_cat <- rbind(outputTable_cat, outputTable)
  }
  
  rmse <- rmse(outputTable$ref, outputTable$pred)
  rmse_vec[kfold] = rmse
  
  print(summary(pred_y))
  print(summary(test_y$last_day_min))
  print(rmse)
  
  # ##
  # pROC_obj <- roc(test_y$label,pred_y,
  #                 # arguments for ci
  #                 ci=TRUE,
  #                 # arguments for plot
  #                 plot=FALSE, grid=TRUE,
  #                 print.auc=TRUE)
  # 
  # auc_vec[kfold] <- pROC_obj$auc
  # 
  # c = coords(pROC_obj, "best", "threshold")
  # threshold_frame[kfold, ] = c
  # 
  # print(pROC_obj$auc)
  
  if (kfold == 1) {
    plot(exp(test_y$last_day_min), exp(pred_y), xlim = c(1, 27), ylim = c(1, 27), cex = 0.4, pch = 16, col = rgb(0,0,0,0.4,maxColorValue = 1))
  } else{
    points(exp(test_y$last_day_min), exp(pred_y), cex = 0.4, pch = 16, col = rgb(0,0,0,0.4,maxColorValue = 1))
  }
  
  ## generate validation accuracy metric
  # pred_acc <- ifelse(pred_y >= c$threshold, 1, 0)
  # cm <- confusionMatrix(as.factor(pred_acc), as.factor(test_y$label))
  # 
  # accuracy_v[kfold] <- cm$overall[1]
  # balanced_accuracy_v[kfold] <- cm$byClass[11]
  # f1[kfold] <- cm$byClass[7]
  
  # if (kfold == k) {
  #   legend(0.2, 0.8, legend=paste0('kf auroc: ', round(auc_vec, 4)), cex=0.8)
  #   legend(0.2, 0.2, legend=paste0('t/sen/spec: ', round(apply(threshold_frame, 2, median), 2)), cex=0.8)
  # }
  
  
}

print(quantile(rmse_vec))
print(quantile(n_vec))

boxplot(exp(outputTable_cat$pred) ~ cut(exp(outputTable_cat$ref), 100), varwidth=T, ylim = c(0, 15))


print(quantile(threshold_frame[,2]))
print(quantile(threshold_frame[,3]))


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


