# download from google drive to downloads folder as .csv

library(data.table)
x <- fread('~/Downloads/BMI_weight_trajectory - current.csv', stringsAsFactors = F)
x$V1 <- as.Date(x$V1, format = '%d/%m/%Y')

plot(x$V1, x$weight)

x = x[70:nrow(x), ]
x$days <- c(1:nrow(x))

plot(x$days, x$weight)

abline(lm(x$weight ~ x$days))
print(lm(x$weight ~ x$days))
