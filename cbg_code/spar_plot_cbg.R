## CBG spar plots
library(data.table)
library(tidyverse)
library(tidyquant)
library(timetk)
library(padr)
library(zoo)
library(viridis)

x <- fread('~/Documents/data/CBGdata/time_series_cohort_5days.csv')
x[, 'N' := .N, by=.(V1)]
x <- x[N>=6]

# split out last day and label
t <- x[V1 == 203446038]
x$day <- as.Date(x$dateTime)
x[, 'flag_last_day' := ifelse(day == max(day), 1, 0), by=.(V1)]
x[, 'label' := ifelse(min(Glu[flag_last_day == 1]) < 6, 1, 0), by=.(V1)]
# remove the last day from the training set
x <- x[flag_last_day == 0]

x[, 'n' := c(1 : .N), by=.(V1)]
print(sum(x[n==1]$label))
print(nrow(x[n==1]))

#
plot(x$Glu)
class(x$dateTime)
plot(x$dateTime, x$Glu, pch = 16, cex = 0.6, col = ifelse(x$label == 0,
                                                          rgb(0,0,0,0.4, maxColorValue = 1),
                                                          rgb(1,0,0,0.4, maxColorValue = 1)))

# create value at 0000 first day (same as first value), and 23.59 (same as last recorded value) on the last day to allow minute fills of data

returnFirstDayTime <- function(sub) {
  first_day <- substr(sub$dateTime[1], 1, 10)
  first_val <- sub$Glu[1]
  
  first_point_time <- paste0(first_day, ' 0:00:01')
  first_point_time <- as.POSIXct(first_point_time, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  
  sub <- data.frame(sub$V1, sub$dateTime, sub$Glu, sub$location, sub$op)
  colnames(sub) <- c('V1', 'dateTime', 'Glu', 'location', 'op')
  
  insert_top <- as.data.frame(matrix(nrow = 1, ncol = 5))
  colnames(insert_top) <- c('V1', 'dateTime', 'Glu', 'location', 'op')
  insert_top$V1 <- sub$V1[1]
  insert_top$dateTime <- first_point_time
  insert_top$Glu <- first_val
  insert_top$location <- sub$location[1]
  insert_top$op <- sub$op[1]
  
  return(insert_top)
  
}

returnLastDayTime <- function(sub) {
  last_day <- substr(sub$dateTime[nrow(sub)], 1, 10)
  last_val <- sub$Glu[nrow(sub)]
  
  last_point_time <- paste0(last_day, ' 23:59:59')
  last_point_time <- as.POSIXct(last_point_time, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  
  sub <- data.frame(sub$V1, sub$dateTime, sub$Glu, sub$location, sub$op)
  colnames(sub) <- c('V1', 'dateTime', 'Glu', 'location', 'op')
  
  insert_end <- as.data.frame(matrix(nrow = 1, ncol = 5))
  colnames(insert_end) <- c('V1', 'dateTime', 'Glu', 'location', 'op')
  insert_end$V1 <- sub$V1[1]
  insert_end$dateTime <- last_point_time
  insert_end$Glu <- last_val
  insert_end$location <- sub$location[nrow(sub)]
  insert_end$op <- sub$op[nrow(sub)]
  
  return(insert_end)
  
}

densityMap <- function(v1, s, i, b, label) {
  
  #v1 = sub$V1[1]; s = jitter(out$Glu); i = 1000
  s <- log(s)
  
  start_point <- 1
  offset <- i
    
    initial_point <- start_point + (2 * offset)
    end_padding <- length(s) - initial_point
    
    x <- s[start_point:end_padding]
    y <- s[(start_point + offset):(end_padding+offset)]
    z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
    
    # from paper 2D projection 
    u <- 1/3 * (x + y + z)
    v <- (1/sqrt(6)) * (x + y - (2 * z))
    w <- (1/sqrt(2)) * (x - y)
    
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    
    v <- range01(v)
    w <- range01(w)
    
    df <- data.frame(v, w)
    
    if (label == 1) {
      png(paste0('~/Documents/data/plots/pos/plot_', v1, '_', label, '.png'), width = 8, height = 8, units = 'in', res = 300)
      p <- ggplot(df, aes(x = v, y = w)) +
        stat_density2d(aes(fill = ..density..), geom = 'tile', contour = F) +
        scale_fill_viridis()
      p <- p + theme_void() + theme(legend.position = "none")
      print(p)
      dev.off()
    } else {
      png(paste0('~/Documents/data/plots/neg/plot_', v1, '_', label, '.png'), width = 8, height = 8, units = 'in', res = 300)
      p <- ggplot(df, aes(x = v, y = w)) +
        stat_density2d(aes(fill = ..density..), geom = 'tile', contour = F) +
        scale_fill_viridis()
      p <- p + theme_void() + theme(legend.position = "none")
      print(p)
      dev.off()
    }
    
  }

idVec <- unique(x$V1)
for (j in c(1:length(idVec))) {
#for (j in c(1:100)) {
  
  id <- idVec[j]
  sub <- x[V1 == id]
  label <- sub$label[1]
  
  insert_top <- returnFirstDayTime(sub)
  insert_end <- returnLastDayTime(sub)
  
  sub <- data.frame(sub$V1, sub$dateTime, sub$Glu, sub$location, sub$op)
  colnames(sub) <- c('V1', 'dateTime', 'Glu', 'location', 'op')
    
  out <- rbind(insert_top, sub, insert_end)
  out <- out %>% thicken("1 min") %>% select(-dateTime) %>% pad()
  
  out$Glu <- na.approx(out$Glu)
  
  #out <- fill(out, c(V1, Glu, location, op))
  
  # pass to density map function
  densityMap(out$V1[1], jitter(out$Glu), 480, 100, label)
  
  }
  
## try this as a classifier