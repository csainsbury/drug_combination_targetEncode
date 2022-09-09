# prep sample of CBG inpatient data for testing
library(data.table)
library(tidyverse)

source('~/Documents/code/cbg_code/cbg_functions.R')

x <- fread('~/Documents/data/CBGdata/aegis_010122-080822_qeuh.csv', header = F)

subset <- x %>% select(V148, V150, V101, V104, V140)
colnames(subset) <- c('CHI', 'location', 'dateTime', 'Glu', 'op')
subset$CHI <- as.numeric(subset$CHI)
subset$Glu <- ifelse(subset$Glu == '<1.1', 1, subset$Glu)
subset$Glu <- ifelse(subset$Glu == '>27.8', 1, subset$Glu)
subset$dateTime <- as.POSIXct(subset$dateTime, format = "%d/%m/%Y %H:%M")

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
abline(h = 10)

# minimum number of admission duration days for study cohort
# will equal train + test period
cohort_selection_threshold <- 10
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

cf[, 'use_flag' := n_days_select(dateTime, cohort_selection_threshold, interval_days), by =.(V1, admission_vec)]

use_cohort <- cf[use_flag == 1]



  