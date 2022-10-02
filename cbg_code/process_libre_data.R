library(data.table)

path <- c('~/Documents/data/raw_libre/')
file_list <- list.files(paste0(path))


returnFirstDayTime <- function(sub) {
  
  first_day <- substr(sub$dateTime[1], 1, 10)
  first_val <- sub$Glu[1]
  
  first_point_time <- paste0(first_day, ' 0:01')
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

  last_day <- substr(sub$dateTime[nrow(sub)], 1, 10)
  last_val <- sub$Glu[nrow(sub)]
  
  last_point_time <- paste0(last_day, ' 23:59')
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
  
  insert_top <- returnFirstDayTime(x)
  insert_end <- returnLastDayTime(x)
  
  x <- rbind(insert_top, x, insert_end)
  
}
