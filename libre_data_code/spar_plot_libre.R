## CBG spar plots
library(data.table)
library(tidyverse)
library(tidyquant)
library(timetk)
library(padr)
library(zoo)
library(viridis)

# ek2_supervised_7days.csv 336 (7 days) points train. 48 points predict
x <- fread('~/Documents/data/libre_data/ek2_supervised_7days.csv')
t_col = 336

train <- x[, 1:t_col]
test <- x[, (t_col + 1):ncol(x)]

min_test <- apply(test,1,min)
label <- ifelse(min_test <3, 1, 0)

train$label <- label
train$id <- c(1:nrow(train))

# take every nth row (to prevent over similarity)
n = 48
train <- train[id%%n==0]

# down sample n
case_n = nrow(train[label ==  1]) # n cases
ratio = 2

cases    <- train[label == 1]
controls <- train[label == 0]

set.seed(42)
select_case <- cases[sample(nrow(cases), case_n), ]
select_control <- controls[sample(nrow(controls), case_n * ratio), ]

sc <- rbind(cases, select_control)

# create value at 0000 first day (same as first value), and 23.59 (same as last recorded value) on the last day to allow minute fills of data

densityMap_diff_plus <- function(v1, s, i, b, label, limit_min, limit_max) {
  
  # #v1 = sub$ID[1]; s = out$Glu; i = 480
  # s <- sc[1,]
  # id <- s$id
  # label <- s$label
  # s$label = NULL
  # s$id = NULL
  
  # jpeg(paste0('~/Documents/data/plots_simple/event_,', label, '.jpg'), width = 800, height = 800, units = 'px')
  # plot(as.numeric(s[1,]))
  # dev.off()
  
  s <- as.numeric(s)
  s <- s^2
  #s = (s - 1) / (27 - 1)
  sd <- c(0, diff(s))
  # sd = (sd - 0.1) / (27 - 0.1)
  # 
  #s <- s * sd
  # s = (s - min(s)) / (max(s) - min(s))
  # 
  start_point <- 1
  offset <- i
  
  initial_point <- start_point + (2 * offset)
  end_padding <- length(s) - initial_point
  
  x <- s[start_point:end_padding]
  y <- s[(start_point + offset):(end_padding+offset)]
  z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]
  
  #plot(s)
  
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
    jpeg(paste0('~/Documents/data/plots_by_class/hypo/', id, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x=v, y=w) ) +
      geom_hex(bins = 40) +
      scale_fill_continuous(type = "viridis") +
      theme_bw() + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  } else {
    jpeg(paste0('~/Documents/data/plots_by_class/nohypo/', id, '.', label, '.jpg'), width = 800, height = 800, units = 'px')
    
    m <- ggplot(df, aes(x=v, y=w) ) +
      geom_hex(bins = 40) +
      scale_fill_continuous(type = "viridis") +
      theme_bw() + theme_void()
    m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
    print(m)
    
    dev.off()
  }
  
  
}


idVec <- unique(sc$id)
for (j in c(1:length(idVec))) {
#for (j in c(1:100)) {
  
  if (j%%100 == 0) {print(j/nrow(sc))}
  
  id <- sc$id[j]
  sub <- sc[j, ]
  label <- sc$label[j]
  
  sub <- sub[, -c('id', 'label')]
  sub <- as.numeric(sub)
  
  #out <- fill(out, c(V1, Glu, location, op))
  
  # pass to density map function
  densityMap_diff_plus(id, sub, 24, 100, label, 0 ,1)
  
  }
  
## try this as a classifier


m <- ggplot(df, aes(x = v, y = w)) +
  geom_point(alpha = 4/10) + theme_void()
m <- m + geom_density_2d() + xlim(c(pl_min, pl_max)) + ylim(c(pl_min, pl_max))
print(m)

ggplot(df, aes(x=v, y=w) ) +
  geom_hex(bins = 60) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() + theme_void()
