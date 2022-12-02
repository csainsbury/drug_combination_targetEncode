# prep sample of CBG inpatient data for testing
## from aegisPOC

library(data.table)
library(tidyverse)
library(padr)

source('~/Documents/code/cbg_code/cbg_functions.R')

days_n  <- 21
minimum_n_cbgs <- days_n + 2
hypo_threshold <- 3
lockout = 7
ratio = -99

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
  x <- x[order(x$uID, x$unix_dateTime), ]
  
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

process_aegis <- function(x) {
  
  subset <- x %>% select(V148, V150, V101, V104, V140)
  colnames(subset) <- c('CHI', 'location', 'dateTime', 'Glu', 'op')
  subset$CHI <- as.numeric(subset$CHI)
  subset <- subset[CHI > 0]
  
  chi_remove = c(1111111111, 2222222222, 3333333333, 4444444444, 5555555555, 6666666666,
                 7777777777, 8888888888, 9999999999, 111111111, 222222222, 333333333,
                 444444444, 555555555, 666666666, 777777777, 888888888, 999999999, 3232323232,
                 2121212121)
  subset <- subset[!subset$CHI %in% chi_remove]
  
  subset$Glu <- ifelse(subset$Glu == '<1.1', 1, subset$Glu)
  subset$Glu <- ifelse(subset$Glu == '>27.8', 27.9, subset$Glu)
  subset$dateTime <- as.POSIXct(subset$dateTime, format = "%d/%m/%Y %H:%M")
  
  subset$numericGlu = as.numeric(subset$Glu)
  subset$numericGlu[is.na(subset$numericGlu)] <- -99
  subset <- subset[numericGlu > 0]
  
  subset$unix_dateTime <- returnUnixDateTime(subset$dateTime)
  
  # remove unwanted locations
  remove_vector <- c('XDUMP', 'RENAL DIALYSI', 'Renal Dialysis', 'RENALDIALYSIS', 'Renal Unit', 'Missing meters',
                     'Missing Metersx', 'RMissing Meters', '44 LABOUR WARD', 'Biochem (K)', 'Biochemistry',
                     'Biochemistry RA', 'Mother\\Baby Uni')
  subset <- subset[!subset$location %in% remove_vector]
  
  # rename CHI to uID
  subset <- subset %>% select(-Glu)
  subset <- subset %>% rename(uID = CHI,
                              Glu = numericGlu)
  
  return(subset)
  
}

x1 <- fread('~/Documents/data/CBGdata/aegis_raw/aegis_010122-080822_qeuh.csv', header = F)
x2 <- fread('~/Documents/data/CBGdata/aegis_raw/aegis_010822_251122.csv', header = F)

x <- rbind(x1, x2)

# process aegis data into correct format for analysis
x <- process_aegis(x)

data.orig <- x

x[, 'admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(uID)]

x <- initial_process_label(x, days_n, minimum_n_cbgs, hypo_threshold)

if (ratio >0 ) {
  x <- case_control_ratio(x, ratio)
}

s <- admission_process(x)

# get data from most recent admission, and all prior admission N
s[, 'qualifying_admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(uID)]

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

    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


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

aegis_model_ready = data.table(export)

## load trained model here
require(xgboost)
require(pROC)
require(caret)

ratio = 4
maxd = 1
model_load_path = paste0(paste0('Documents/data/CBGdata/savedModels/SAT_21_FULL_wholeData_modelSAVE_', days_n,
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
