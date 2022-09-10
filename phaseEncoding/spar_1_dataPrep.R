library(data.table)

returnUnixDateTime<-function(date) {
  returnVal<-as.numeric(as.POSIXct(date, format="%d-%m-%Y %H:%M", tz="GMT"))
  return(returnVal)
}

my_data <- fread('~/projects11/3d2d/MarkSmoker_glucose_29-6-2022.csv', skip = 2)
my_data <- fread('~/projects11/3d2d/EileenKesson_glucose_28-6-2022.csv', skip = 2)

my_data$unixTime <- returnUnixDateTime(my_data$`Device Timestamp`)

my_data <- my_data[order(my_data$unixTime), ]
my_data$`Historic Glucose mmol/L`[is.na(my_data$`Historic Glucose mmol/L`)] <- 0 
my_data <- my_data[`Historic Glucose mmol/L` > 0]

# organise into day cycles
m = my_data
x <- substr(m$`Device Timestamp`, 1, 10)
m$Date <- as.Date(x, format = c("%d-%m-%Y"))
m[, 'day_n' := seq(1, .N, 1), by=.(Date)]

# distribution of n tests/day
m[, 'n_measures' := max(day_n), by=.(Date)]

# 95 or 96 measures are complete days
m <- m[n_measures == 95 | n_measures == 96]

# add overall row number
m$diff_date <- c(1, diff(m$Date))
m$cumsum_diff <- cumsum(m$diff_date)

start_val = 1
i = 1

report_vector <- 1

increment = 1

for (j in c(2:nrow(m))) {
  
  if (diff(c(m$cumsum_diff[j], m$cumsum_diff[j-1])) == 0 |
      diff(c(m$cumsum_diff[j], m$cumsum_diff[j-1])) == 1) {
    
    increment = 0
    
  } else {
    
    increment = 1
    
  }
  
  report_vector <- c(report_vector, report_vector[j-1] + increment)
  
  
  
}




out <- data.frame(m$unixTime, m$`Historic Glucose mmol/L`)
colnames(out) <- c('datetime', 'glucose')

write.csv(out, '~/projects11/3d2d/export_2.csv', row.names = F)
