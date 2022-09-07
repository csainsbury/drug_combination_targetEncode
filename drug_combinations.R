library(data.table)
library(tidyverse)
library(caret)

x <- fread('~/Documents/data/x_out_3Mbins.csv')
x <- x[, 1:(ncol(x) - 1)]

x <- data.table(x)
x$PrescriptionDateTime <- as.Date(x$PrescriptionDateTime)

ids <- unique(x$LinkId)

for (i in c(1:length(ids))) {
  
  if(i%%1000 == 0) {print(paste0(i, ' of ', length(ids)))}
  
  sub <- x[LinkId == ids[i]]
  sub[, 'comb' := as.character(paste0(attr(table(DrugName), 'dimnames'))), by=.(PrescriptionDateTime)]
  
  t.first <- sub[match(unique(sub$PrescriptionDateTime), sub$PrescriptionDateTime),]
  t.first <- t.first %>% select(LinkId, PrescriptionDateTime, comb)
  t.first <- t.first[order(t.first$PrescriptionDateTime), ]
  
  if (i == 1) {
    export <- t.first
  } else {
    export <- rbind(export, t.first)
  }
  
}

## export represents a combination per month for all drugs prescribed for 15k IDs
write.table(export, file = paste0('~/Documents/data/export_drugCombinationsPerID_3Mbins.csv'), sep = ',', row.names = F)

export <- fread(paste0('~/Documents/data/export_drugCombinationsPerID_3Mbins.csv'))

# add a column (n_bins) that records the number of time bins that a drug combination is prescribed for
export[, 'n_bins' := .N, by=.(LinkId, comb)]
# number of time bins for which there is data per ID
export[, 'bins_of_data' := .N, by=.(LinkId)]

  # investigate how much data is available per ID
  e <- export
  e[, 'n' := c(1: .N), by =. (LinkId)]
  e <- e[n==1]
  h <- hist(e$bins_of_data, 100)

# explore most common combinations
x <- as.data.frame(table(export$comb))
x <- x[order(-x$Freq), ]

plot(x$Freq, xlim = c(0, 1000), ylim = c(0, 200))
abline(a = 100, b = 0)

## import outcome data
death <- fread('~/Documents/data/diagnosisDateDeathDate.txt', stringsAsFactors = F)
death$DeathDate <- as.Date(death$DeathDate)

last_date <- max(death$DeathDate, na.rm = T)
death <- death %>% select(LinkId, DeathDate)

study_group <- merge(export, death, by.x = 'LinkId', by.y = 'LinkId', all.x = T)

## define time periods
# train:   (last_date - (5y + 18 months)) to (last_date - 18 months)
# outcome: (last_date - 18 months) to last_date

# target encoding function
# from:: https://medium.com/@darnelbolaos/target-encoding-function-with-r-8a037b219fb7
target_enc <- function(data, enc_col, tar_col, k_col, kmin, kmax){
  require(tidyverse)
  col_name <- paste0("tar_enc_", enc_col)
  
  temp <- map_df(kmin:kmax, function(k) {
    
    xtrain <- data[data[k_col] != k, ]
    xvalid <- data[data[k_col] == k, ]
    
    feat <- xtrain %>% group_by_at(enc_col) %>%
      summarise_at(.vars = tar_col, mean)
    
    colnames(feat)[2] <- col_name
    temp_df <- xvalid %>% left_join(feat, by = enc_col) %>% 
      select(all_of(enc_col), all_of(col_name))
    
    return(temp_df)
    
  })
  
  temp_enc <- temp %>% group_by_at(enc_col) %>% 
    summarise_at(.vars = col_name, .funs = mean) %>% 
    ungroup()
  
  data %>% left_join(temp_enc, by = enc_col)
  
}

train_years = 5
outcome_years = 2

train_cut <- study_group[PrescriptionDateTime >= (last_date- ((train_years*365.25) + (outcome_years*365.25)))]

  # require alive at end train
  end_train_date <- (last_date - (outcome_years*365.25))
  train_cut$alive_flag <- ifelse(is.na(train_cut$DeathDate) == TRUE, 1, 0)
  train_cut$include <- ifelse(train_cut$alive_flag == 1 | train_cut$DeathDate >= end_train_date, 1, 0)
  train_cut <- train_cut[include == 1]
  # require an entry within the last (end_buffer_months) months
  end_buffer_months = 6
  buffer_cut <- end_train_date - (end_buffer_months * ((365.25/12)))
  train_cut[, 'flag_include' := ifelse(max(PrescriptionDateTime) > buffer_cut, 1, 0), by =.(LinkId)]
  train_cut <- train_cut[flag_include == 1]

  # add outcome flag
  train_cut$outcome <- ifelse(train_cut$DeathDate >= (last_date - (outcome_years*365.25)), 1, 0)
  train_cut$outcome[is.na(train_cut$outcome)] <- 0

  # replace all commas in the drug combination col
  train_cut$comb <- gsub(',', ';', train_cut$comb, ignore.case=TRUE)

  # ensure ID is numeric
  train_cut$LinkId <- as.numeric(train_cut$LinkId)
  
    # investigate how much data is available per ID
    e <- train_cut
    e[, 'n' := c(1: .N), by =. (LinkId)]
    e <- e[n==1]
    h <- hist(e$bins_of_data, 100)
  
  # simplify
  train_cut <- train_cut %>% select(LinkId, PrescriptionDateTime, comb, outcome)

## kfold split in here
  
target_encode_lookup <- function(dat, min_threshold = 10) {
  #reduce to single occurence of each drug combination per ID
  data_ids <- unique(dat$LinkId)
  for (j in c(1:length(data_ids))) {
    sub <- dat[LinkId == data_ids[j]]
    first_occurence <- sub[match(unique(sub$comb), sub$comb),]
    
    if (j == 1) {
      output <- first_occurence
    } else {
      output <- rbind(output, first_occurence)
    }
  }
  
  # target encoding
  data <- as_tibble(output)
  data <- data %>% mutate(comb = fct_lump_min(comb, min_threshold)) # 10 gives reasonable balance
  n_kfold <- 20 # target encoding kfold
  
  data <- data %>% mutate(k = sample(1:n_kfold, nrow(.), replace = TRUE))
  
  x <- target_enc(data = data, enc_col = "comb",
                  tar_col = "outcome", k_col = "k",
                  kmin = 1, kmax = n_kfold)
  
  hist(x$tar_enc_comb, 500, ylim = c(0, 1000))
  
  # from here can generate a lookup table of the target encoding value per combination
  lookup <- x[match(unique(x$comb), x$comb),]
  lookup <- lookup[order(-lookup$tar_enc_comb), ]
  
  lookup <- lookup %>% select(comb, tar_enc_comb)
  lookup <- as.data.frame(lookup)
  
  print(lookup$comb[1:10])
  print(lookup$tar_enc_comb[1:10])
  
  plot(lookup$tar_enc_comb)
  return(lookup)
}
  
  
      # kfold creation
      set.seed(100)
      k = 10
      
      # frame with IDs and outcome
      train_ids <- train_cut
      train_ids[, 'n' := c(1: .N), by=.(LinkId)]
      train_ids <- train_ids[n == 1]
      train_ids <- train_ids %>% select(LinkId, outcome)
      
      # define training control
      flds <- createFolds(train_ids$outcome, k, list = T, returnTrain = F)
      
      for (kfold in c(1:k)) {
        
        # define train/test by id
        dat  <- train_ids[-flds[[kfold]], ]
        test <- train_ids[flds[[kfold]], ]
        # expand to include all data
        dat <- merge(dat, train_cut, by.x = c('LinkId', 'outcome'), by.y = c('LinkId', 'outcome'))
        test <- merge(test, train_cut, by.x = c('LinkId', 'outcome'), by.y = c('LinkId', 'outcome'))
        
        lookup <- target_encode_lookup(dat, 60)
        
        # add target encoded values to training data
        dat <- merge(dat, lookup, by.x = 'comb', by.y = 'comb')
        
      }

# # do mean target encoding in R
# # summarise the rare categories
# n_rare_limit <- 5
# x <- as.data.frame(table(train_cut$comb))
# x <- x[order(-x$Freq), ]
# plot(x$Freq, xlim = c(0, 1000), ylim = c(0, 200))
# abline(a = n_rare_limit, b = 0)
# 
# print(which(x$Freq < n_rare_limit)[1])


    # # convert train_cut to tibble
    # data_train <- train_cut
    # 
    # #reduce to single occurence of each drug combination per ID
    # data_ids <- unique(data_train$LinkId)
    # for (j in c(1:length(data_ids))) {
    #   sub <- data_train[LinkId == data_ids[j]]
    #   first_occurence <- sub[match(unique(sub$comb), sub$comb),]
    #   
    #   if (j == 1) {
    #     output <- first_occurence
    #   } else {
    #     output <- rbind(output, first_occurence)
    #   }
    # }
    
        # # target encoding
        # data <- as_tibble(output)
        # 
        # data <- data %>% mutate(comb = fct_lump_min(comb, 10)) # 10 gives reasonable balance
        # 
        # n_kfold <- 20
        # 
        # data <- data %>% mutate(k = sample(1:n_kfold, nrow(.), replace = TRUE))
        # 
        # x <- target_enc(data = data, enc_col = "comb",
        #                 tar_col = "outcome", k_col = "k",
        #                 kmin = 1, kmax = n_kfold)
        # 
        # hist(x$tar_enc_comb, 1000)
        
        # # from here can generate a lookup table of the target encoding value per combination
        # lookup <- x[match(unique(x$comb), x$comb),]
        # lookup <- lookup[order(-lookup$tar_enc_comb), ]
        # 
        # head(lookup, 20)
        # tail(lookup, 20)

        
        
# export represents a combination per month for all drugs prescribed for 15k IDs
write.table(train_cut, file = paste0('~/Documents/data/analysis_data.csv'), sep = ',', row.names = F)


