## glucose steroid redo in R
library(data.table)
library(tidyverse)

source('~/Documents/code/cbg_code/cbg_functions.R')

x1 <- fread('~/Documents/data/glucose_steroid/queh_010122-010522.csv', header = T,fill = T)
x2 <- fread('~/Documents/data/glucose_steroid/queh_010521-010122.csv', header = T, fill = T)

x <- rbind(x1, x2)

x <- x[`Result Type` == "Patient"]

s <- data.table('uID' = x$`Patient ID`,
                'datetime' = x$`Test Date/Time`,
                'Glu' = x$GLU,
                'loc' = x$Location)

s$uID <- as.numeric(s$uID)
s$uID[is.na(s$uID)] <- 0
s <- s[uID > 0]
s <- s[nchar(s$uID) >= 9]

s <- s[uID != 1111111111 &
         uID != 2222222222 &
         uID != 3333333333 &
         uID != 4444444444 &
         uID != 5555555555 &
         uID != 6666666666 &
         uID != 7777777777 &
         uID != 8888888888 &
         uID != 9999999999 &
         uID != 111111111 &
         uID != 222222222 &
         uID != 333333333 &
         uID != 444444444 &
         uID != 555555555 &
         uID != 666666666 &
         uID != 777777777 &
         uID != 888888888 &
         uID != 999999999 &
         uID != 0000000001 &
         uID != 9999999993
]

s$Glu <- ifelse(s$Glu == '<1.1', 1, s$Glu)
s$Glu <- ifelse(s$Glu == '>27.8', 27.9, s$Glu)
s$Glu <- as.numeric(s$Glu)
s$Glu[is.na(s$Glu)] <- 0
s <- s[Glu > 0]

s$dateTime <- as.POSIXct(s$datetime, format = "%d/%m/%Y %I:%M%p", tz = 'GMT')

s <- unique(s)
#plot(s$dateTime, s$Glu, cex = 0.2, col = rgb(0,0,0,0.1, maxColorValue = 1))

# identify each admission using 7 day lockout period
# number of admissions (unix date time required):
m <- s
m <- m[order(m$uID, m$dateTime), ]

lockout = 7
m$unix_dateTime <- returnUnixDateTime(m$dateTime)
m[, 'n_admissions' := admission_N(unix_dateTime, lockout), by=.(uID)]
# label each admission
m[, 'admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(uID)]

# add interval last-first measurement
m[, 'interval_seconds' := as.numeric(max(unix_dateTime) - min(unix_dateTime)), by=.(uID, admission_vec)]
m[, 'interval_days' := as.numeric(max(as.Date(dateTime)) - min(as.Date(dateTime))), by=.(uID, admission_vec)]
m <- m[order(m$uID, m$dateTime), ]
m[, 'n' := c(1:.N), by=.(uID, admission_vec)]
m[, 'N' := .N, by=.(uID, admission_vec)]

single_m <- m[n==1]

single_m <- single_m[order(-single_m$interval_days), ]

## load steroid / dm drug data
dmst <- fread('~/Documents/data/glucose_steroid/dm_drugs5.csv')
dmst_ids <- as.numeric(unique(dmst$CHI))

dmst$earliest_date <- as.Date('1970-01-01')
# loop to sort Earliest dates:
for (j in c(1:nrow(dmst))) {
  date = dmst$Earliest[j]
  
  if (nchar(date) == 10) {
    dmst$earliest_date[j] <- as.Date(date, format = '%d/%m/%Y')
  } else {
    dmst$earliest_date[j] <- as.Date(date, format = '%d/%m/%y')
  }
  
}

# manage DOB
dmst$DOB_ <- as.Date(dmst$DOB, format = '%m/%d/%y')
# 2010 as cutoff : from : hist(dmst$DOB_, 100)
dmst$DOB_ <- ifelse(dmst$DOB_ > '2010-01-01', dmst$DOB_ - (100 * 365.25), dmst$DOB_)
dmst$DOB_ <- as.Date(dmst$DOB_, origin = '1970-01-01')

dmst$DOB <- dmst$DOB_
dmst$DOB_ <- NULL

write.table(dmst, file = '~/Documents/data/glucose_steroid/dmst.csv', sep = ',', row.names = F)
write.table(dmst, file = '~/Documents/data/glucose_steroid/WED/dmst.csv', sep = ',', row.names = F)

## in both
both = single_m[single_m$uID %in% dmst_ids]
both[, 'both_admission_n' := c(1:.N), by=.(uID)]
both[, 'both_admission_N' := .N, by=.(uID)]

## flag last admission per ID
both[, 'flag_last_admission' := ifelse(both_admission_n == max(both_admission_n), 1, 0), by=.(uID)]

## select last admissions only
require(chron)
both_last <- both[flag_last_admission==1]
both_last$time = times(substr(both_last$dateTime, 12,19))
both_last$day = as.Date(both_last$dateTime)
## set rule that if first CBG after 6pm, day of admission includes next day
both_last$day_advance <- ifelse(both_last$time > '18:00:00', 1, 0)
both_last$effective_day <- ifelse(both_last$day_advance == 1, (both_last$day + 1), both_last$day)
both_last$effective_day <- as.Date(both_last$effective_day, origin = '1970-01-01')

## for each last admission find the closest set of prescriptions
both_last$match = 0

for (k in c(1:nrow(both_last))) {
  
  uID = both_last$uID[k]
  effective_day <- both_last$effective_day[k]
  
  dmst_sub <- dmst[CHI == uID]
  if (
    length(which(dmst_sub$earliest_date == effective_day)) == 1 |
      length(which(dmst_sub$earliest_date == effective_day - 1)) == 1
      ) {
    print('yes')
    both_last$match[k] = 1
    } else {
      print('no')
      both_last$match[k] = 0
    }
}

dt = dmst
dt = dt %>% select(-c('V1', 'Name', 'Drug 1 End Date', 'Drug 2 End Date', 'Drug 3 End Date', 'Drug 4 End Date',
                        'Drug 5 End Date', 'Drug 6 End Date', 'Drug 7 End Date', 'Drug 8 End Date'))
dt = dt[, 1:19]

for (j in c(1:nrow(dt))) {
  
  sub = dt[j, ]
  new =  data.table('chi' = sub$CHI, 'drug' = sub$`Drug 1 Name`, date = sub$`Drug 1 Start Date`)
  new1 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 2 Name`, date = sub$`Drug 2 Start Date`)
  new2 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 3 Name`, date = sub$`Drug 3 Start Date`)
  new3 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 4 Name`, date = sub$`Drug 4 Start Date`)
  new4 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 5 Name`, date = sub$`Drug 5 Start Date`)
  new5 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 6 Name`, date = sub$`Drug 6 Start Date`)
  new6 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 7 Name`, date = sub$`Drug 7 Start Date`)
  new7 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 8 Name`, date = sub$`Drug 8 Start Date`)
  new8 = data.table('chi' = sub$CHI, 'drug' = sub$`Drug 9 Name`, date = sub$`Drug 9 Start Date`)
  
  new_out = rbind(new, new1, new2, new3, new4, new5, new6, new7, new8)
  
  if (j==1) {
    out = new_out
  } else {
    out = rbind(out, new_out)
  }
  
}

long = out
long$date = as.Date(long$date, format = '%d/%m/%Y')
long$drug <- ifelse(long$drug == '', 0, long$drug)
long$drug[is.na(long$drug)] <- 0

long = long[long$drug != '0']
long$drug <- as.factor(long$drug)

long$drug_new <- as.character(long$drug)
long$drug_new[grep("ABASAGLAR", long$drug, ignore.case = TRUE)] <- "glargine"
long$drug_new[grep("lantus", long$drug, ignore.case = TRUE)] <- "glargine"
long$drug_new[grep("toujeo", long$drug, ignore.case = TRUE)] <- "glargine"
long$drug_new[grep("apidra", long$drug, ignore.case = TRUE)] <- "apidra"
long$drug_new[grep("dulaglutide", long$drug, ignore.case = TRUE)] <- "dulaglutide"
long$drug_new[grep("tresiba", long$drug, ignore.case = TRUE)] <- "dulaglutide"
long$drug_new[grep("fiasp", long$drug, ignore.case = TRUE)] <- "novorapid"
long$drug_new[grep("novorapid", long$drug, ignore.case = TRUE)] <- "novorapid"
long$drug_new[grep("humalog", long$drug, ignore.case = TRUE)] <- "humalog"
long$drug_new[grep("Xultophy", long$drug, ignore.case = TRUE)] <- "ideglira"


table(long$drug_new)

long$drug = long$drug_new
long$drug_new <- NULL


# newd = one_hot(as.data.table(long[, 1:2]))
# 
# summary data per admission (for prior admission info)
m[, 'max_by_admission' := max(Glu), by=.(uID, admission_vec)]
m[, 'min_by_admission' := min(Glu), by=.(uID, admission_vec)]
m[, 'median_by_admission' := median(Glu), by=.(uID, admission_vec)]
m[, 'iqr_by_admission' := quantile(Glu)[4]-quantile(Glu)[2], by=.(uID, admission_vec)]

m[, 'gradient_by_admission' := as.numeric(lm(Glu ~ dateTime)$coefficients[2]), by=.(uID, admission_vec)]
m[, 'sd_cbg_by_admission' := sd(Glu) , by=.(uID, admission_vec)]
m[, 'mean_cbg_by_admission' := mean(Glu) , by=.(uID, admission_vec)]
m[, 'cV_by_admission' := sd_cbg_by_admission/mean_cbg_by_admission , by=.(uID, admission_vec)]

plot(m[n==1]$cV_by_admission, m[n==1]$iqr_by_admission, cex = 0.2, col=rgb(0,0,0,0.3))
# plot(m[n==1]$cV_by_admission, m[n==1]$median_by_admission, cex = 0.2, col=rgb(0,0,0,0.3))

# build X, y
bd <- both_last[match==1]

# index admission glucose params and y
indexAdmissionCBGs <- merge(bd, m, by.x=c('uID', 'admission_vec'), by.y = c('uID', 'admission_vec'))
indexAdmissionCBGs <- indexAdmissionCBGs %>% select(uID, admission_vec, day_advance, effective_day, dateTime.y, Glu.y)

indexAdmissionCBGs$day.y <- as.Date(indexAdmissionCBGs$dateTime.y)

indexAdmissionCBGs[, 'info_period' := ifelse(day.y == effective_day, 1, 0), by =.(uID, admission_vec)]
indexAdmissionCBGs[, 'info_period_extra' := ifelse(day_advance == 1 & day.y == (effective_day - 1), 1, 0), by =.(uID, admission_vec)]

indexAdmissionCBGs$combined_info_period <- ifelse(indexAdmissionCBGs$info_period == 1 | indexAdmissionCBGs$info_period_extra == 1, 1, 0)

# summary data per admission (for prior admission info)
indexAdmissionCBGs[, 'max_by_admission' := max(Glu.y), by=.(uID, admission_vec, combined_info_period)]
indexAdmissionCBGs[, 'min_by_admission' := min(Glu.y), by=.(uID, admission_vec, combined_info_period)]
indexAdmissionCBGs[, 'median_by_admission' := median(Glu.y), by=.(uID, admission_vec, combined_info_period)]
indexAdmissionCBGs[, 'iqr_by_admission' := quantile(Glu.y)[4]-quantile(Glu.y)[2], by=.(uID, admission_vec, combined_info_period)]

indexAdmissionCBGs[, 'gradient_by_admission' := as.numeric(lm(Glu.y ~ dateTime.y)$coefficients[2]), by=.(uID, admission_vec, combined_info_period)]
indexAdmissionCBGs[, 'sd_cbg_by_admission' := sd(Glu.y) , by=.(uID, admission_vec, combined_info_period)]
indexAdmissionCBGs[, 'mean_cbg_by_admission' := mean(Glu.y) , by=.(uID, admission_vec, combined_info_period)]
indexAdmissionCBGs[, 'cV_by_admission' := sd_cbg_by_admission/mean_cbg_by_admission , by=.(uID, admission_vec, combined_info_period)]

# create long data where all drugs prescribed during info period
long_ids <- unique(long$chi)
for (k in c(1:length(long_ids))) {
  
  indexSub <- indexAdmissionCBGs[uID == long_ids[k]]
  
  if (nrow(indexSub) > 0) {
    
    indexSub = indexSub[combined_info_period == 1]
    interest_days <- unique(indexSub$day.y)
    
    long_sub <- long[chi == indexSub$uID[1]]
    long_sub <- long_sub[date %in% interest_days]
    
    if (k==1) {
      long_out = long_sub
    } else {
      long_out = rbind(long_out, long_sub)
    }
    
  }

}
  
setDT(long_out[, 1:2])
multihot = dcast(setDT(melt(long_out[, 1:2],id.vars = c("chi")))[,ind:=1],chi~value,value.var = "ind",fill=0)
multihot$chi = as.numeric(multihot$chi)

# merge index_admission info with drug info and add label
indexAdmissionCBGs$label <- ifelse(indexAdmissionCBGs$combined_info_period == 0 & indexAdmissionCBGs$max_by_admission >= 20, 1, 0)
indexAdmissionCBGs[, 'label_max' := max(label), by=.(uID)]
  
# split out index info period data with label
cbg_data = indexAdmissionCBGs[combined_info_period == 1]
cbg_data[, 'n' := c(1:.N), by=.(uID)]

cbg_data_single = cbg_data[n==1]

cbg = merge(cbg_data_single, multihot, by.x = 'uID', by.y = 'chi')

## add features from prior admissions:
cbg$n_prior_admissions = cbg$admission_vec - 1

flag=1
for (r in  c(1:nrow(cbg))) {
  print(r)
  print(cbg$n_prior_admissions[r])
  
  if (cbg$n_prior_admissions[r] > 0) {
    required_admission_number = cbg$n_prior_admissions[r]
    sub = m[uID == cbg$uID[r] & admission_vec == required_admission_number]
    sub = sub %>% select(uID, max_by_admission, min_by_admission, median_by_admission, iqr_by_admission,
                         gradient_by_admission, sd_cbg_by_admission, mean_cbg_by_admission, cV_by_admission)
    out = sub[1,]
    colnames(out) = paste0('prior_', colnames(out))
    colnames(out)[1] <- c('uID')
    
    if (flag == 1) {
      prior_out = out
    } else {
      prior_out = rbind(prior_out, out)
    }
    flag = flag + 1
  }
}

final = merge(cbg, prior_out, by.x = 'uID', by.y = 'uID', all.x = T)

final_export <- final %>% select(-admission_vec, -day_advance, -effective_day, -dateTime.y, -Glu.y, -day.y, -info_period, -info_period_extra, -combined_info_period, -label)

write.table(final_export, file = '~/Documents/data/glucose_steroid/new_final.csv', sep = ',', row.names = F)


# final set


