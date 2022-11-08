library(data.table)
library(tidyverse)
library(padr)
library(zoo)

path <- c('~/Documents/data/raw_libre/')
savepath <- c('~/Documents/data/processed_libre/')
file_list <- list.files(paste0(path))


returnFirstDayTime <- function(sub) {
  
  sub$Glu[is.na(sub$Glu)] <- 0
  sub <- sub[Glu>0]
  
  first_day <- substr(sub$dateTime[1], 1, 10)
  first_val <- sub$Glu[1]
  
  first_point_time <- paste0(first_day, ' 0:00')
  first_point_time <- as.POSIXct(first_point_time, format = "%Y-%m-%d %H:%M", tz = 'UTC')
  
  sub <- data.frame(sub$dateTime, sub$Glu)
  colnames(sub) <- c('dateTime', 'Glu')
  
  insert_top <- as.data.frame(matrix(nrow = 1, ncol = 2))
  colnames(insert_top) <- c('dateTime', 'Glu')
  insert_top$dateTime <- first_point_time
  insert_top$Glu <- first_val
  
  return(insert_top)
  
}

returnLastDayTime <- function(sub) {
  
  sub$Glu[is.na(sub$Glu)] <- 0
  sub <- sub[Glu>0]

  last_day <- substr(sub$dateTime[nrow(sub)], 1, 10)
  last_val <- sub$Glu[nrow(sub)]
  
  last_point_time <- paste0(last_day, ' 23:45')
  last_point_time <- as.POSIXct(last_point_time, format = "%Y-%m-%d %H:%M", tz = 'UTC')
  
  sub <- data.frame(sub$dateTime, sub$Glu)
  colnames(sub) <- c('dateTime', 'Glu')
  
  insert_end <- as.data.frame(matrix(nrow = 1, ncol = 2))
  colnames(insert_end) <- c('dateTime', 'Glu')
  insert_end$dateTime <- last_point_time
  insert_end$Glu <- last_val
  
  return(insert_end)
  
}

for (i in c(1:length(file_list))) {
  
  x <- fread(paste0(path, file_list[i]))
  
  x$dateTime <- as.POSIXct(x$`Device Timestamp`, format = "%d-%m-%Y %H:%M")
  x <- data.table('dateTime' = x$dateTime,
                  'Glu' = x$`Historic Glucose mmol/L`)
  
  x <- x[order(x$dateTime), ]
  
  insert_top <- returnFirstDayTime(x)
  insert_end <- returnLastDayTime(x)
  
  x <- rbind(insert_top, x, insert_end)
  x <- data.table(x)
  
  x$dateTime <- as.POSIXct(x$dateTime, format = "%d-%m-%Y %H:%M")
  
  x <- x[is.na(x$dateTime) == FALSE]
  
  x$Glu[is.na(x$Glu)] <- 0
  x <- x[Glu>0]
  
  x <- x %>% thicken("15 min") %>% select(-dateTime) %>% pad()
  x$Glu <- na.approx(x$Glu)
  
  x$Glu[is.na(x$Glu)] <- 0
  x[Glu>0]
  
  plot(x$dateTime_15_min, x$Glu, cex = ifelse(x$Glu<=3.9, 0.2, ifelse(x$Glu>20, 0.2, 0)), col = ifelse(x$Glu<=3.9, rgb(1,0,0,0.4), 'black')); lines(x$dateTime_15_min, x$Glu, col=rgb(0,0,0,0.2))
  
  name <- paste0('proc_', i)
  write.table(x, file = paste0(savepath, name, '.csv'), sep = ',', row.names = F)
  
}
