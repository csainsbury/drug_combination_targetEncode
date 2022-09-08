# prep sample of CBG inpatient data for testing

library(data.table)
library(tidyverse)

x <- fread('~/Documents/data/CBGdata/aegis_010122-080822_qeuh.csv', header = F)

subset <- x %>% select(V148, V150, V101, V104, V140)
colnames(subset) <- c('CHI', 'location', 'dateTime', 'Glu', 'op')
subset$CHI <- as.numeric(subset$CHI)
subset$Glu <- ifelse(subset$Glu == '<1.1', 1, subset$Glu)
subset$Glu <- ifelse(subset$Glu == '>27.8', 1, subset$Glu)
subset$dateTime <- as.POSIXct(subset$dateTime, format = "%d/%m/%Y %H:%M")

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
head(m)

  id_sort <- ids[order(-ids$N), ]
  id = id_sort$V1[4]
  plot(m[V1 == id]$dateTime, m[V1 == id]$Glu, col = as.factor(m[V1 == id]$location), pch = 16)
  lines(m[V1 == id]$dateTime, m[V1 == id]$Glu, lwd = 0.4)
  print(m[V1 == id]$location)
  
  