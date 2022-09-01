# prep sample of CBG inpatient data for testing

library(data.table)
library(tidyverse)

x <- fread('~/Documents/data/CBGdata/aegis_010122-080822_qeuh.csv', header = F)

subset <- x %>% select(V148, V150, V101, V104)
colnames(subset) <- c('CHI', 'location', 'dateTime', 'Glu')

