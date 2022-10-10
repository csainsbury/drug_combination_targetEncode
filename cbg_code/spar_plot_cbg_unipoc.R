## CBG spar plots
library(data.table)
library(tidyverse)
library(tidyquant)
library(timetk)
library(padr)
library(zoo)
library(viridis)
library(imputeTS)

# 101355572

dataset <- 'huge'
days_n  <- 10
minimum_n_cbgs <- days_n + 2

x <- fread(paste0('~/Documents/data/CBGdata/', dataset, '_unipoc_time_series_cohort_first_', days_n,'_days.csv'))
x <- x[order(x$uID, x$admission_vec)]

exclude <- c('dialysis', 'renal', 'abbott', 'Ketones', 'X Ray', 'X-ray', 'training', 'staff', 'returned', 'obsolete', 'old', 'Meterx', 'dead', 'BioChemistry', 'weldon', 'radiology', 'outpatient', 'rcas')
x_ex <- x
for (e in c(1:length(exclude))) {
  x_ex <- x_ex[-grep(exclude[e], x_ex$loc, ignore.case = TRUE), ]
}
x <- x_ex

#library(dplyr)
# add unique ID by 2 cols - uID and admission vec
x <- x %>% mutate(ID = group_indices(x, .dots=c("uID", "admission_vec"))) 
# x[, 'N' := .N, by=.(uID)]

## here need to truncate all admissions to the right length!!
x$day <- as.Date(x$dateTime)
x[, 'correct_duration_days' := ifelse(day <= (min(day) + (days_n - 1)), 1, 0), by=.(ID)]
x <- x[correct_duration_days == 1]

# limit to those with minimum n CBGs
x[, 'N_truncated' := .N, by=.(ID)]
x <- x[N_truncated >= minimum_n_cbgs]

# split out last day and label
#t <- x[V1 == 203446038]
#x$day <- as.Date(x$dateTime)
x[, 'flag_last_day' := ifelse(day == max(day), 1, 0), by=.(ID)]
x[, 'label' := ifelse(min(Glu[flag_last_day == 1]) <= 3, 1, 0), by=.(ID)]
# remove the last day from the training set
x <- x[flag_last_day == 0]

x[, 'N' := .N, by=.(ID)]
x <- x[N>1]

x[, 'n' := c(1 : .N), by=.(ID)]
print(sum(x[n==1]$label))
print(nrow(x[n==1]))

# plot(x$dateTime, x$Glu, pch = 16, cex = 0.6, col = ifelse(x$label == 0,
#                                                           rgb(0,0,0,0.4, maxColorValue = 1),
#                                                           rgb(1,0,0,0.4, maxColorValue = 1)))

# produce fixed ratio of case to no case
ratio = 2
event_ids <- unique(x[label==1]$ID)
no_event_ids <- unique(x[label==0]$ID)
id_sample <- no_event_ids[sample(length(no_event_ids), round(length(event_ids) * ratio), 0)]

neg <- x[ID %in% id_sample]
pos <- x[ID %in% event_ids]

x <- rbind(neg, pos)

print(sum(x[n==1]$label))
print(nrow(x[n==1]))

# create value at 0000 first day (same as first value), and 23.59 (same as last recorded value) on the last day to allow minute fills of data

returnFirstDayTime <- function(sub) {
  first_day <- substr(sub$dateTime[1], 1, 10)
  first_val <- sub$Glu[1]
  
  first_point_time <- paste0(first_day, ' 0:00:01')
  first_point_time <- as.POSIXct(first_point_time, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  
  sub <- data.frame(sub$ID, sub$dateTime, sub$Glu, sub$loc)
  colnames(sub) <- c('ID', 'dateTime', 'Glu', 'loc')
  
  insert_top <- as.data.frame(matrix(nrow = 1, ncol = 4))
  colnames(insert_top) <- c('ID', 'dateTime', 'Glu', 'loc')
  insert_top$ID <- sub$ID[1]
  insert_top$dateTime <- first_point_time
  insert_top$Glu <- first_val
  insert_top$location <- sub$location[1]
  #insert_top$op <- sub$op[1]
  
  return(insert_top)
  
}

returnLastDayTime <- function(sub) {
  last_day <- substr(sub$dateTime[nrow(sub)], 1, 10)
  last_val <- sub$Glu[nrow(sub)]
  
  last_point_time <- paste0(last_day, ' 23:59:59')
  last_point_time <- as.POSIXct(last_point_time, format = "%Y-%m-%d %H:%M:%S", tz = 'UTC')
  
  sub <- data.frame(sub$ID, sub$dateTime, sub$Glu, sub$loc)
  colnames(sub) <- c('ID', 'dateTime', 'Glu', 'loc')
  
  insert_end <- as.data.frame(matrix(nrow = 1, ncol = 4))
  colnames(insert_end) <- c('ID', 'dateTime', 'Glu', 'loc')
  insert_end$ID <- sub$ID[1]
  insert_end$dateTime <- last_point_time
  insert_end$Glu <- last_val
  insert_end$location <- sub$location[nrow(sub)]
  #insert_end$op <- sub$op[nrow(sub)]
  
  return(insert_end)
  
}

densityMap <- function(v1, s, i, b, label) {
  
  #v1 = sub$ID[1]; s = out$Glu; i = 1000
  s <- diff(s)
  
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
    
    # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    # 
    # v <- range01(v)
    # w <- range01(w)
    
    df <- data.frame(v, w)
    
    pixel_n = 200
    
    if (label == 1) {
      jpeg(paste0('~/Documents/data/plots/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
      p <- ggplot(df, aes(x = v, y = w)) +
        stat_density2d(aes(fill = ..density..), geom = 'tile', n = pixel_n, contour = F) +
        scale_fill_viridis() + 
      p <- p + theme_void() + theme(legend.position = "none")
      print(p)
      dev.off()
    } else {
      jpeg(paste0('~/Documents/data/plots/nil.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
      p <- ggplot(df, aes(x = v, y = w)) +
        stat_density2d(aes(fill = ..density..), geom = 'tile', n = pixel_n, contour = F) +
        scale_fill_viridis()
      p <- p + theme_void() + theme(legend.position = "none")
      print(p)
      dev.off()
    }
    
}

densityMap_diff <- function(v1, s, i, b, label) {
  
  #v1 = sub$ID[1]; s = out$Glu; i = 480
  s <- diff(s)
  
  start_point <- 1
  offset <- i
  
  initial_point <- start_point + (2 * offset)
  end_padding <- length(s) - initial_point
  
  x <- s[start_point:end_padding]
  y <- s[(start_point + offset):(end_padding+offset)]
  z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
  
  # library(rgl); plot3d(x,y,z)
  
  # from paper 2D projection 
  u <- 1/3 * (x + y + z)
  v <- (1/sqrt(6)) * (x + y - (2 * z))
  w <- (1/sqrt(2)) * (x - y)
  
  # plot(v, w, cex = 0.4, pch=16)
  
  # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  # 
  # v <- range01(v)
  # w <- range01(w)
  
  df <- data.frame(v, w)
  
  # pixel_n = 200
  
  pl_min = -0.1
  pl_max = 0.1
  
  if (label == 1) {
    jpeg(paste0('~/Documents/data/plots/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point() + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  } else {
    jpeg(paste0('~/Documents/data/plots/nil.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point() + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  }
  
}

densityMap_diff_plus <- function(v1, s, i, b, label, limit_min, limit_max) {
  
  #v1 = sub$ID[1]; s = out$Glu; i = 480
  s <- s
  
      jpeg(paste0('~/Documents/data/plots_simple/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
          plot(s)
      dev.off()
  
  
  sd <- c(0, diff(s))
  
  s <- s * sd
  
  start_point <- 1
  offset <- i
  
  initial_point <- start_point + (2 * offset)
  end_padding <- length(s) - initial_point
  
  x <- s[start_point:end_padding]
  y <- s[(start_point + offset):(end_padding+offset)]
  z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
  
  # library(rgl); plot3d(x,y,z)
  
  # from paper 2D projection 
  u <- 1/3 * (x + y + z)
  v <- (1/sqrt(6)) * (x + y - (2 * z))
  w <- (1/sqrt(2)) * (x - y)
  
  # plot(v, w, cex = 0.4, pch=16)
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  v <- range01(v)
  w <- range01(w)
  
  df <- data.frame(v, w)
  
  pixel_n = 200
  
  pl_min = limit_min
  pl_max = limit_max
  
  if (label == 1) {
    jpeg(paste0('~/Documents/data/plots/event.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')

    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point(alpha = 4/10) + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)

    dev.off()
  } else {
    jpeg(paste0('~/Documents/data/plots/nil.', v1, '.', label, '.jpg'), width = 800, height = 800, units = 'px')

    m <- ggplot(df, aes(x = v, y = w)) +
      geom_point(alpha = 4/10) + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)

    dev.off()
  }

  
}

idVec <- unique(x$ID)
print(length(idVec))
for (j in c(1:length(idVec))) {
#for (j in c(1:100)) {
  
  if (j%%100 == 0) {print(j)}
  
  id <- idVec[j]
  sub <- x[ID == id]
  label <- sub$label[1]
  
  insert_top <- returnFirstDayTime(sub)
  insert_end <- returnLastDayTime(sub)
  
  # truncate to maximum 30 day duration of data
  max_day <- as.Date(insert_top$dateTime) + 31
  sub <- sub[as.Date(sub$dateTime) <= max_day]
  
  sub <- data.frame(sub$ID, sub$dateTime, sub$Glu, sub$loc)
  colnames(sub) <- c('ID', 'dateTime', 'Glu', 'loc')
  
  # deal with instances where all cbgs the same
  if (prod(diff(sub$Glu)) == 0) {
    sub$Glu[1] <- jitter(sub$Glu[1])
  } 
    
  out <- rbind(insert_top, sub, insert_end)
  
      # silence warnings for the padding script
      oldw <- getOption("warn")
      options(warn = -1)
      
        out <- out %>% thicken("1 min") %>% select(-dateTime) %>% pad()
      
      on.exit(options(warn = oldw))
  
  # out$Glu <- imputeTS::na_locf(out$Glu)
  out$Glu <- na.approx(out$Glu)
  
  #out <- fill(out, c(V1, Glu, location, op))
  
  # pass to density map function
  # densityMap_diff(out$ID[1], jitter(out$Glu), 480, 100, label)
  # densityMap_diff_plus(out$ID[1], jitter(out$Glu), 240, 100, label, 0, 1)
  densityMap_diff_plus(out$ID[1], jitter(out$Glu, amount = 10e-10), 240, 100, label, 0, 1)
  
  }
  
## try this as a classifier