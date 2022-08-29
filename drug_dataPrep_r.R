## all drugs prediction

library(data.table)
library(tidyverse)

d1 <- fread('~/Documents/data/first_out.csv', stringsAsFactors = F, fill = TRUE, header = T, quote = "")
d1 <- d1 %>% dplyr::select(LinkId, PrescriptionDateTime, DrugName)
d1 <- unique(d1)

d2 <- fread('~/Documents/data/second_out.csv', stringsAsFactors = F, fill = TRUE, header = T, quote = "")
d2 <- d2 %>% dplyr::select(LinkId, PrescriptionDateTime, DrugName)
d2 <- unique(d2)

d3 <- fread('~/Documents/data/third_out.csv', stringsAsFactors = F, fill = TRUE, header = T, quote = "")
d3 <- d3 %>% dplyr::select(LinkId, PrescriptionDateTime, DrugName)
d3 <- unique(d3)

# for first up to 500
# for second and third up to 6000

t <- rbind(d1, d2, d3)

t$PrescriptionDateTime <- as.Date(t$PrescriptionDateTime)

earlyCut <- as.Date('2010-01-01')
lateCut <- as.Date('2018-01-01')
t <- t[PrescriptionDateTime > earlyCut & t$PrescriptionDateTime < lateCut]

  x <- as.data.frame(table(t$DrugName))
  x <- x[order(x$Freq, decreasing = T),]
  
  threshold_Rx_n <- 10000
  limit_row <- which(x$Freq <= threshold_Rx_n)[1]
  limit_vector <- x$Var1[1:limit_row]
  
limit_t <- t[t$DrugName %in% limit_vector,]

write.table(limit_t, file = paste0('~/Documents/data/limit_Rxn_threshold_', threshold_Rx_n, '.csv'), sep = ',', row.names = F)
