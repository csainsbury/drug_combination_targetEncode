library(data.table)
library(tidyverse)

# train_val_test <- c(0.6,0.2,0.2)

path <- '~/Documents/data/plots/'

allfiles <- list.files(path)

class1 <- allfiles[substr(allfiles, 1, 3) == 'nil']
class2 <- allfiles[substr(allfiles, 1, 3) == 'eve']

class1L <- length(class1)
class2L <- length(class2)

r <- 1

class1 <- class1[sample(class1L, class2L * r)]



my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

move_files <- function(file_list, pathfrom, pathto) {
  
  for (i in c(1:length(file_list))) {
    filename = file_list[i]
    
    my.file.rename(from = paste0(pathfrom, filename),
                   to = paste0(pathto, filename))
  }
  
}


splitf <- function(data, train_val_test) {
  
  #data = class2
  
  d <- length(data)
  numbers <- floor(d * train_val_test)
  print(numbers)
  train_sample <- data[sample(d, numbers[1])]
      data <- data[-match(train_sample, data)]
  d <- length(data)
  val_sample <- data[sample(d, numbers[2])]
      data <- data[-match(val_sample, data)]
  d <- length(data)
  test_sample <- data[sample(d, numbers[2], replace = T)]
      data <- data[-match(val_sample, data)]
      
  move_files(train_sample, '~/Documents/data/plots/', '~/Documents/data/plots_output/train/')
  move_files(val_sample, '~/Documents/data/plots/', '~/Documents/data/plots_output/val/')
  move_files(test_sample, '~/Documents/data/plots/', '~/Documents/data/plots_output/test/')

}

splitf(class2, c(0.7, 0.25, 0.05))
splitf(class1, c(0.7, 0.25, 0.05))
