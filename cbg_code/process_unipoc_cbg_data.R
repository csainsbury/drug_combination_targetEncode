# process unipoc data for spar plotting

library(data.table)
library(tidyverse)

source('~/Documents/code/cbg_code/cbg_functions.R')

x1 <- fread('~/Documents/data/CBGdata/queh_010122-010522.csv', header = T,fill = T)
x2 <- fread('~/Documents/data/CBGdata/queh_010521-010122.csv', header = T, fill = T)

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

s$dateTime <- as.POSIXct(s$datetime, format = "%d/%m/%Y %I:%M%p")

s <- unique(s)
#plot(s$dateTime, s$Glu, cex = 0.2, col = rgb(0,0,0,0.1, maxColorValue = 1))

# identify each admission using 7 day lockout period
# number of admissions (unix date time required):
m <- s

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

# minimum number of admission duration days for study cohort
# will equal train + test period

cohort_selection_threshold <- 5

plot(single_m$interval_days, cex = 0); lines(single_m$interval_days, lwd = 3)
abline(h = cohort_selection_threshold, col = 'red')

last_id <- which(single_m$interval_days == cohort_selection_threshold)[length(which(single_m$interval_days == cohort_selection_threshold))]

cohort <- single_m[1:last_id, ]
print(dim(cohort))

cohort = data.table('uID' = cohort$uID,
                    'admission_vector' = cohort$admission_vec)

cf <- merge(cohort, m,
            by.x = c('uID', 'admission_vector'),
            by.y = c('uID', 'admission_vec'))

# recalculate admission duration intervals
cf[, 'new_interval_seconds' := max(unix_dateTime) - min(unix_dateTime), by =.(uID, admission_vector)]

# label admissions within cf
cf[, 'admission_vec' := admission_N_vector(unix_dateTime, lockout), by=.(uID)]

u_IDs <- unique(cf$uID)

options(warn = 1)  

export_concat <- cf[1,]
export_concat <- export_concat[-1,]
for (j in c(1:length(u_IDs))) {
  
  # if (j %% 1000 == 0) {print(j)}
  print(j)
  
  sub = cf[uID == u_IDs[j]]
  
  # find each admission
  admissionsList <- as.data.frame(table(sub$admission_vec))$Var1
  
  for (m in c(1:length(admissionsList))) {
    subsub <- sub[admission_vec == admissionsList[m]]
    
    duration <- as.numeric(max(as.Date(subsub$dateTime)) - min(as.Date(subsub$dateTime)))
    
    if (duration > cohort_selection_threshold) {
      
      dateList <- as.Date(subsub$dateTime)
        last_start_date <- max(dateList) - cohort_selection_threshold
          available_dates <- dateList[which(dateList <= last_start_date)]
          random_available_date <- available_dates[sample(length(available_dates), 1)]
          
          start_window <- random_available_date
          end_window <- random_available_date + cohort_selection_threshold
          
          export <- subsub[as.Date(dateTime) >= start_window &
                             as.Date(dateTime) <= end_window]
          
          export_concat <- rbind(export_concat, export)
          
    } else {
      export_concat <- rbind(export_concat, subsub)
    }
  }
  
}


use_cohort <- export_concat

write.table(use_cohort, file = paste0('~/Documents/data/CBGdata/unipoc_time_series_cohort_', cohort_selection_threshold, 'days.csv'), sep = ',', row.names = F)
