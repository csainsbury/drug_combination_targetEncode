# prep sample of CBG inpatient data for testing

library(data.table)
library(tidyverse)

x <- fread('~/Documents/data/CBGdata/aegis_010122-080822_qeuh.csv', header = F)

subset <- x %>% select(V148, V150, V101, V104)
colnames(subset) <- c('CHI', 'location', 'dateTime', 'Glu')

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
