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


