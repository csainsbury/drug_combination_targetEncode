# prep sample of CBG inpatient data for testing
## from aegisPOC

library(data.table)
library(tidyverse)
library(padr)

source('~/Documents/code/cbg_code/cbg_functions.R')

days_n  <- 4
minimum_n_cbgs <- days_n + 2
hypo_threshold <- 3

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

x1 <- fread('~/Documents/data/CBGdata/aegis_raw/aegis_010122-080822_qeuh.csv', header = F)
x2 <- fread('~/Documents/data/CBGdata/aegis_raw/aegis_010822_251122.csv', header = F)

x <- rbind(x1, x2)

subset <- x %>% select(V148, V150, V101, V104, V140)
colnames(subset) <- c('CHI', 'location', 'dateTime', 'Glu', 'op')
subset$CHI <- as.numeric(subset$CHI)
subset$Glu <- ifelse(subset$Glu == '<1.1', 1, subset$Glu)
subset$Glu <- ifelse(subset$Glu == '>27.8', 27.9, subset$Glu)
subset$dateTime <- as.POSIXct(subset$dateTime, format = "%d/%m/%Y %H:%M")

subset$numericGlu = as.numeric(subset$Glu)
subset$numericGlu[is.na(subset$numericGlu)] <- -99
subset <- subset[numericGlu > 0]

data.orig <- subset


# remove inappropriate locations
subset <- subset[location != 'RENALDIALYSIS']
# may add ITU to this

ids <- as.data.table(table(subset$CHI))
  ids$V1 <- as.numeric(ids$V1)
  ids$V1[is.na(ids$V1)] <- 0
  ids <- ids[V1>0]
  
ids <- ids[V1 != 1111111111 &
             V1 != 2222222222 &
             V1 != 3333333333 &
             V1 != 4444444444 &
             V1 != 5555555555 &
             V1 != 6666666666 &
             V1 != 7777777777 &
             V1 != 8888888888 &
             V1 != 9999999999 &
             V1 != 111111111 &
             V1 != 222222222 &
             V1 != 333333333 &
             V1 != 444444444 &
             V1 != 555555555 &
             V1 != 666666666 &
             V1 != 777777777 &
             V1 != 888888888 &
             V1 != 999999999
             ]

ids <- ids[nchar(V1) >8]
id_start <- 10000
ids$uID <- c(id_start:(id_start + (nrow(ids) - 1)))

m <- merge(ids, subset, by.x = 'V1', by.y = 'CHI')
m <- m[order(-m$N, m$dateTime), ]

# required output
# "uID", "admission_vector", "datetime", "Glu", "loc", "dateTime", "unix_dateTime", "n_admissions", "interval_seconds", "interval_days", "n", "N", "new_interval_seconds", "admission_vec"   

## for compatibility-> reshape here into unipoc format

# identify each admission using 7 day lockout period
# number of admissions (unix date time required):
lockout = 7
m$unix_dateTime <- returnUnixDateTime(m$dateTime)
m[, 'n_admissions' := admission_N(unix_dateTime, lockout), by=.(V1)]
# label each admission
m[, 'admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(V1)]

# add interval last-first measurement
m[, 'interval_seconds' := as.numeric(max(dateTime) - min(dateTime)), by=.(V1, admission_vec)]
m[, 'interval_days' := as.numeric(max(as.Date(dateTime)) - min(as.Date(dateTime))), by=.(V1, admission_vec)]
m <- m[order(m$V1, m$dateTime), ]
m[, 'n' := c(1:.N), by=.(V1, admission_vec)]

single_m <- m[n==1]

single_m <- single_m[order(-single_m$interval_days), ]
plot(single_m$interval_days)
abline(h = days_n)

# minimum number of admission duration days for study cohort
# will equal train + test period
cohort_selection_threshold <- days_n
last_id <- which(single_m$interval_days == cohort_selection_threshold)[length(which(single_m$interval_days == cohort_selection_threshold))]

cohort <- single_m[1:last_id, ]
print(dim(cohort))

cf <- merge(cohort, m,
                     by.x = c('V1', 'admission_vec'),
                     by.y = c('V1', 'admission_vec'))
cf <- cf %>% select(V1, N.y, uID.y, location.y, dateTime.y, Glu.y, op.y, unix_dateTime.y, n_admissions.y, interval_seconds.y, interval_days.y, n.y)
colnames(cf) <- c('V1', 'N', 'uID', 'location', 'dateTime', 'Glu', 'op', 'unix_dateTime', 'n_admissions', 'interval_seconds', 'interval_days', 'n')

# label admissions within cf
cf[, 'admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(V1)]

cf[, 'use_flag' := n_days_select(dateTime, cohort_selection_threshold, interval_days, V1, admission_vec), by =.(V1, admission_vec)]

use_cohort <- cf[use_flag == 1]
use_cohort$uID <- NULL # remove old uID col
use_cohort <- rename(use_cohort, uID = V1) # readd new uID - old CHI

write.table(use_cohort, file = paste0('~/Documents/data/CBGdata/aegis_processed/time_series_cohort_', cohort_selection_threshold, 'days.csv'), sep = ',', row.names = F)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### pipeline for running through model
x <- use_cohort
x$Glu <- as.numeric(x$Glu)

#library(dplyr)
# add unique ID by 2 cols - uID and admission vec
x <- x %>% mutate(ID = group_indices(x, .dots=c("uID", "admission_vec"))) 
# x[, 'N' := .N, by=.(uID)]

## here need to truncate all admissions to the right length!!
x$day <- as.Date(x$dateTime, tz = 'GMT')
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
# ratio = 4
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

s <- x

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

hist(s[ID_n==1]$cV, 60)

s[, 'sum_diff_day' := as.numeric(sum(diff(day))), by=.(ID)]

## per ID generate wide data
ids <-unique(s$ID)
for (i in c(1:length(ids))) {
  
  if (i %% 100 == 0) {print(i/length(ids))}
  
  sub <- s[ID == ids[i]]
  sub <- sub[day_n == 1]
  
  sub <- sub %>% select(ID, dateTime, day, Glu, max_by_day, min_by_day, median_by_day, day_N, iqr_by_day, cV_by_day, training_min, training_max, cV, gradient, label)
  
  if ((nrow(sub) > 0) & (as.numeric(min(diff(as.Date(sub$dateTime)))) <= 1)) { # second arguement needed to ensure that minimum time step is not greater than 1 day (thicken will not work if it is)
    
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
    label <- sub$label[1]
    
    sub <- sub %>% select(dateTime, Glu, max_by_day, min_by_day, median_by_day, day_N, iqr_by_day, cV_by_day)
    sub <- sub %>% thicken("1 day") %>% select(-dateTime) %>% pad()
    # hack to remove any duplicate days due to GMT/BST issues
    diff_days = c(diff(sub$dateTime_day), 1)
    sub = sub[diff_days == 1, ]
    
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
    m$label <- label
    m <- m[-1, ]
    colnames(m) <- c(m[1,1:(ncol(m) - 6)],'cV' ,'gradient' , 'training_min', 'training_max','id', 'label')
    m <- m[-1, ]
    
    if (i == 1) {
      export <- m
    } else {
      if (ncol(export) == ncol(m)) {
        export <- rbind(export, m)
      }
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

aegis_model_ready = data.table(export)

## load trained model here
require(xgboost)
require(pROC)
require(caret)

ratio = 4
maxd = 1
model_load_path = paste0(paste0('Documents/data/CBGdata/savedModels/wholeData_modelSAVE_', days_n,
                                '_days_hypothresh_', hypo_threshold,
                                '_ratio_', ratio,
                                '_maxd_', maxd,'.model'))
trained_model = xgb.load(model_load_path)

aegis_model_run <- data.matrix(aegis_model_ready[, -c('id', 'label')])

pred_y = predict(trained_model, aegis_model_run)

##
pROC_obj <- roc(aegis_model_ready$label,pred_y,
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

hist(pred_y, 100)

## generate validation accuracy metric
pred_acc <- ifelse(pred_y >= c$threshold, 1, 0)
cm <- confusionMatrix(as.factor(pred_acc), as.factor(aegis_model_ready$label), positive = '1')
print(threshold)
print(cm)

print(cm$byClass[7])

# 0.2350176

training_threshold = 0.235
pred_acc <- ifelse(pred_y >= training_threshold, 1, 0)
cm <- confusionMatrix(as.factor(pred_acc), as.factor(aegis_model_ready$label), positive = '1')
print(training_threshold)
print(cm)

review_preds <- data.table(cbind(export, pred_y))
review_preds[label==1]

## compare with any prior hypo
cmh <- confusionMatrix(as.factor(ifelse(aegis_model_ready$training_min < 4, 1, 0)), as.factor(aegis_model_ready$label), positive = '1')
print(cmh)
