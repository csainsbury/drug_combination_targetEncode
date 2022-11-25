require(berryFunctions)
require(harmonicmeanp)

## original functions

# functions
returnUnixDateTime<-function(date) {
  returnVal<-as.numeric(as.POSIXct(date, format="%d/%m/%Y %I:%M%p", tz="GMT"))
  return(returnVal)
}

# return location from xls file for unknown serials
locationReturn <- function(loc_label) {
  # loc_label <- x$Location
  location <- loc_label[1]
  if (is.na(match(location, locsOfInterest)) == TRUE) {output <- rep("no match", length(loc_label))}
  if (is.na(match(location, locsOfInterest)) == FALSE) {output <- rep(location, length(loc_label))}
  return(output)
}

# return number of admissions per ID
admission_N <- function(datetime, admissionThresholdDays) {
  datetime <- datetime[order(datetime)]
  difference <- diff(datetime)
  flag <- ifelse(difference > (60*60*24*admissionThresholdDays), 1, 0)
  admission_N <- sum(flag) + 1
  return(admission_N)
}

# return vector admission number per ID
admission_N_vector <- function(datetime, admissionThresholdDays=7) {
  datetime <- datetime[order(datetime)]
  difference <- diff(datetime)
  flag <- ifelse(difference > (60*60*24*admissionThresholdDays), 1, 0)
  flagVec <- c(1, flag)
  return(cumsum(flagVec))
}

# hypo episodes function
admission_hypoEpisodes <- function(datetime, glu, hypoThresh, hypoDurationThresholdMins, ID) {
  # datetime = dataDT[Patient.ID == 1508366284]$dateplustime1; hypoThresh = 4; glu = dataDT[Patient.ID == 1508366284]$GLUn; hypoDurationThresholdMins = 60
  hypoFrame <- data.frame(datetime, glu)
  # hypoFrame$episode_cc <- 0
  
  if (nrow(hypoFrame) > 1) {
    hypoFrame$ishypo <- ifelse(hypoFrame$glu < hypoThresh, 1, 0)
    hypoFrame$isNextCBG_hypo <- c(hypoFrame$ishypo[2:nrow(hypoFrame)], 0)
    hypoFrame$timeToPrior_mins <- c(0, (diff(hypoFrame$datetime) / (60)))
    hypoFrame$timeToNext_mins <- c((diff(hypoFrame$datetime) / (60)), 0)
    
    hypoFrame$endOfHypoEpisode <- ifelse((hypoFrame$ishypo == 1 & hypoFrame$isNextCBG_hypo == 0) |
                                           (hypoFrame$ishypo == 1 & hypoFrame$isNextCBG_hypo == 1 & hypoFrame$timeToNext_mins > hypoDurationThresholdMins), 1, 0)
    if (hypoFrame$ishypo[nrow(hypoFrame)] == 1) {hypoFrame$endOfHypoEpisode[nrow(hypoFrame)] = 1}
    hypoFrame$resolvedHypoEpisode <- ifelse(hypoFrame$ishypo == 1 & hypoFrame$isNextCBG_hypo == 0 & hypoFrame$timeToNext_mins < hypoDurationThresholdMins, 1, 0)
    
    n_hypo_episodes <- sum(hypoFrame$endOfHypoEpisode)
    n_hypo_episodes_resolved <- sum(hypoFrame$resolvedHypoEpisode)
    
    endOfHypoEpisode_output <- hypoFrame$endOfHypoEpisode
    resolvedHypoEpisode_output <- hypoFrame$resolvedHypoEpisode
    
  }
  
  if (nrow(hypoFrame) == 1 & hypoFrame$glu[1] < hypoThresh) {
    n_hypo_episodes = 1
    n_hypo_episodes_resolved = 0
    endOfHypoEpisode_output = 1
    resolvedHypoEpisode_output = 0
  }
  
  if (nrow(hypoFrame) == 1 & hypoFrame$glu[1] >= hypoThresh) {
    n_hypo_episodes = 0
    n_hypo_episodes_resolved = 0
    endOfHypoEpisode_output = 0
    resolvedHypoEpisode_output = 0
  } 
  
  hypoResolve_cc <- ifelse(endOfHypoEpisode_output == 1 & resolvedHypoEpisode_output == 1, 1, (ifelse(endOfHypoEpisode_output == 1 & resolvedHypoEpisode_output == 0, -1, 0)))
  
  output <- list(endOfHypoEpisode_output, resolvedHypoEpisode_output, hypoResolve_cc, n_hypo_episodes, n_hypo_episodes_resolved)
  
  return(output)
}

# TTR
admissionTTR <- function(datetime, glu) {
  # datetime = dataDT[Patient.ID == 0101366108]$dateplustime1; hypoThresh = 4; glu = dataDT[Patient.ID == 0101366108]$GLUn; hypoDurationThresholdMins = 60
  hypoFrame <- data.frame(datetime, glu)
  
  if (nrow(hypoFrame) > 1) {
    hypoFrame$timeToNext_mins <- c((diff(hypoFrame$datetime) / (60)), 0)
  }
  
  if (nrow(hypoFrame) == 1) {
    hypoFrame$timeToNext_mins <- 0
  }
  
  return(hypoFrame$timeToNext_mins)
}

# TTR - flagging for control chart
admissionTTR_controlChart <- function(datetime, glu, hypoThresh, TTR_standard) {
  # datetime = dataDT[Patient.ID == 0101366108]$dateplustime1; hypoThresh = 4; glu = dataDT[Patient.ID == 0101366108]$GLUn; hypoThresh = 4; TTR_standard = 60
  hypoFrame <- data.frame(datetime, glu)
  
  if (nrow(hypoFrame) > 1) {
    hypoFrame$timeToNext_mins <- c((diff(hypoFrame$datetime) / (60)), 0)
  }
  
  if (nrow(hypoFrame) == 1) {
    hypoFrame$timeToNext_mins <- 0
  }
  
  hypoFrame$flag_opportunity <- ifelse(hypoFrame$glu < hypoThresh, 1, 0)
  hypoFrame$cc_score <- 0
  hypoFrame$cc_score <- ifelse(hypoFrame$flag_opportunity == 1 & hypoFrame$timeToNext_mins < TTR_standard, 1, hypoFrame$cc_score)
  hypoFrame$cc_score <- ifelse(hypoFrame$flag_opportunity == 1 & hypoFrame$timeToNext_mins >= TTR_standard, -1, hypoFrame$cc_score)
  
  output <- list(hypoFrame$timeToNext_mins, hypoFrame$flag_opportunity, hypoFrame$cc_score)
  
  return(output)
  
}

priorMeasurementTime <- function(datetime) {
  # datetime = dataDT[serial  == "KAAT339D0230"]$dateplustime1
  # datetime = dataDT[locBySerial  == "47"]$dateplustime1
  
  dateTable <- data.table(datetime)
  dateTable$reportingOrder <- c(1: length(datetime))
  
  dateTable <- dateTable[order(dateTable$datetime), ]
  
  
  if (nrow(dateTable) > 1) {
    dateTable$timeFromPriorMeasurement <- c(0, (diff(dateTable$datetime) / (60)))
  }
  
  if (length(datetime) == 1) {
    timeFromPriorMeasurement <- 0
  }
  
  dateTable <- dateTable[order(dateTable$reportingOrder), ]
  output <- dateTable$timeFromPriorMeasurement
  
  return(output)
  
}

funnelPlotPrep <- function(inputFrame) {
  # inputFrame = dataDT[dateplustime1 < returnUnixDateTime("01/03/2019 12:00AM")]
  inputFrame[, c("admissionN_vector") := admission_N_vector(dateplustime1, 5), by=.(Patient.ID)]
  inputFrame[, c("CBGn_by_admission") := seq(1, .N, 1), by =. (Patient.ID, admissionN_vector)]
  inputFrame <- inputFrame[CBGn_by_admission == 1]
  inputFrame[, c("admissionN_perLocation") := seq(1, .N, 1), by =. (locBySerial)]
  inputFrame[, c("hypoEps_perLocation") := sum(hypoEps), by =. (locBySerial)]
  inputFrame[, c("resolvedHypoEps_perLocation") := sum(resolvedHypoEps), by =. (locBySerial)]
  
  return(inputFrame)
}

controlChartFunction <- function(inputDT, metric, type, plotname, xlim_lower, xlim_upper, opacity){
  
  # inputDT = dataDT[locBySerial == "4A"]; metric = "resolve"; type = "unit"; plotname = "test"; xlim_lower = -200; xlim_upper = 400; opacity = 0.4
  
  if(metric == "TTR") {
    TTR_sub <- inputDT[TTR_opportunity == 1]
  }
  
  if(metric == "resolve") {
    TTR_sub <- inputDT[endOfHypoEpisode == 1]
  }
  
  if (type == "unit") {units <- as.data.frame(table(TTR_sub$locBySerial))}
  if (type == "operator") {units <- as.data.frame(table(TTR_sub$Operator.ID))}
  if (type == "studyFlag") {units <- as.data.frame(table(TTR_sub$unitStudyFlag))}
  
  pdf(file = paste0("~/projects/sanofi_inpatient/plots/", metric, "_", plotname, ".pdf"), width=16, height=9)
  
  #plotting loop
  for (j in seq(1, nrow(units), 1)) {
    
    if (type == "unit") {TTR_unit_sub <- TTR_sub[locBySerial == units$Var1[j]]}
    if (type == "operator") {TTR_unit_sub <- TTR_sub[Operator.ID == units$Var1[j]]}
    if (type == "studyFlag") {TTR_unit_sub <- TTR_sub[unitStudyFlag == units$Var1[j]]}
    
    TTR_unit_sub <- TTR_unit_sub[order(TTR_unit_sub$dateplustime1), ]
    
    if(metric == "TTR") {
      TTR_unit_sub$cc <- cumsum(TTR_unit_sub$TTR_cc_score)
    }
    
    if(metric == "resolve") {
      TTR_unit_sub$cc <- cumsum(TTR_unit_sub$hypoResolve_cc)
    }
    
    print(paste("unit: ", units$Var1[j], " total: ", TTR_unit_sub$cc[nrow(TTR_unit_sub)], sep = ""))
    
    nhypos = TTR_unit_sub$isHypo_byLocation[1]
    lineW = sqrt(nhypos) * 10e-2
    
    if (j == 1) {
      plot(TTR_unit_sub$dateplustime1,
           TTR_unit_sub$cc,
           xlim = c(min(dataDT$dateplustime1), max(dataDT$dateplustime1)),
           ylim = c(xlim_lower, xlim_upper),
           cex = 0,
           lwd = lineW,
           main = paste0("control chart for ", type," - success in repeating CBG within ", TTR_assessment_time, " mins post CBG < ", hypo_threshold))
      lines(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc)
      text(max(TTR_unit_sub$dateplustime1) + 1000000, TTR_unit_sub$cc[nrow(TTR_unit_sub)], labels = units$Var1[j], cex = 0.6)
    }
    
    if (j > 1) {
      points(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, cex = 0)
      lines(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, lwd = lineW, col = ifelse(TTR_unit_sub$locBySerial == "5B" | TTR_unit_sub$locBySerial == "5A", rgb(1, 0, 0, opacity, maxColorValue = 1), ifelse(TTR_unit_sub$locBySerial == "ARU2", rgb(0, 1, 0, opacity, maxColorValue = 1), rgb(0, 0, 0, opacity, maxColorValue = 1))))
      text(max(TTR_unit_sub$dateplustime1) + 1000000, TTR_unit_sub$cc[nrow(TTR_unit_sub)], labels = units$Var1[j], cex = 0.6)
      abline(v = returnUnixDateTime("01/03/2019 12:00AM"), col = "red", lwd = 2)
      abline(v = min(TTR_unit_sub$dateplustime1), col = rgb(0, 0, 0, 0.2))
      
    }
    
  }
  
  dev.off()
  
}

normalisedControlGradient <- function(input_1, input_2, timeBins, metric, type, plotname, opacity){
  
  # input_1 = dataDT[locBySerial == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")]; input_2 = dataDT[locBySerial == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")]; timeBins = 1; metric = "TTR"; type = "unit"; plotname = "aa_test_dddd"; xlim_lower = -200; xlim_upper = 400; opacity = 0.4
  
  ## setup
  
  if(metric == "TTR") {
    TTR_sub_1 <- input_1[TTR_opportunity == 1]
    TTR_sub_2 <- input_2[TTR_opportunity == 1]
  }
  
  if(metric == "resolve") {
    TTR_sub_1 <- input_1[endOfHypoEpisode == 1]
    TTR_sub_2 <- input_2[endOfHypoEpisode == 1]
  }
  
  if (type == "unit") {
    units_1 <- as.data.frame(table(TTR_sub_1$locBySerial))
    units_2 <- as.data.frame(table(TTR_sub_2$locBySerial))
  }
  if (type == "operator") {
    units_1 <- as.data.frame(table(TTR_sub_1$Operator.ID))
    units_2 <- as.data.frame(table(TTR_sub_2$Operator.ID))
  }
  if (type == "studyFlag") {
    units_1 <- as.data.frame(table(TTR_sub_1$unitStudyFlag))
    units_2 <- as.data.frame(table(TTR_sub_2$unitStudyFlag))
  }
  
  ### normalise gradients
  
  ## max / min time value in either input
  maxTime = max(max(input_1$dateplustime1), max(input_2$dateplustime1))
  minTime = min(min(input_1$dateplustime1), min(input_2$dateplustime1))
  
  ## time bins
  binSeconds = (maxTime - minTime) / timeBins
  
  ## time loop
  for (ii in seq(0, timeBins - 1, 1)) {
    
    pdf(file = paste0("~/projects/sanofi_inpatient/plots/normalised_", metric, "-", plotname, "_",ii, ".pdf"), width=16, height=9)
    
    print(ii)
    
    timeSub_1 <- TTR_sub_1[dateplustime1 >= (minTime + (binSeconds * (ii))) &
                             dateplustime1 < (minTime + (binSeconds * (ii+1)))]
    timeSub_2 <- TTR_sub_2[dateplustime1 >= (minTime + (binSeconds * (ii))) &
                             dateplustime1 < (minTime + (binSeconds * (ii+1)))]
    
    timeSub_1 <- timeSub_1[order(timeSub_1$dateplustime1), ]
    timeSub_2 <- timeSub_2[order(timeSub_2$dateplustime1), ]
    
    # sandboxing
    if(metric == "TTR") {
      gradient_1 = cumsum(timeSub_1$TTR_cc_score)[nrow(timeSub_1)] / nrow(timeSub_1)
      timeSub_1$gradient <- gradient_1
      
      gradient_2 = cumsum(timeSub_2$TTR_cc_score)[nrow(timeSub_2)] / nrow(timeSub_2)
      timeSub_2$gradient <- gradient_2
    }
    
    if(metric == "resolve") {
      gradient_1 = cumsum(timeSub_1$hypoResolve_cc)[nrow(timeSub_1)] / nrow(timeSub_1)
      timeSub_1$gradient <- gradient_1
      
      gradient_2 = cumsum(timeSub_2$hypoResolve_cc)[nrow(timeSub_2)] / nrow(timeSub_2)
      timeSub_2$gradient <- gradient_2
    }
    
    if (nrow(timeSub_1) > nrow(timeSub_2)) {
      smallSub <- timeSub_2
      largeSub <- timeSub_1
    }
    if (nrow(timeSub_2) > nrow(timeSub_1)) {
      smallSub <- timeSub_1
      largeSub <- timeSub_2
    }
    
    ## proportion of +1 vs total for each and chi square
    if(metric == "TTR") {
      small_pos = subset(as.data.frame(table(smallSub$TTR_cc_score)), Var1 == 1)$Freq
      large_pos = subset(as.data.frame(table(largeSub$TTR_cc_score)), Var1 == 1)$Freq
    }
    if(metric == "resolve") {
      small_pos = subset(as.data.frame(table(smallSub$hypoResolve_cc)), Var1 == 1)$Freq
      large_pos = subset(as.data.frame(table(largeSub$hypoResolve_cc)), Var1 == 1)$Freq
    }
    # prop test
    prop = prop.test(c(small_pos, large_pos), c(nrow(smallSub), nrow(largeSub)))
    
    mainText <- paste0(metric, " | time bin ", ii, ". smaller unit (black), n=", nrow(smallSub)," : ", smallSub$locBySerial[1], " - gradient: ", round(smallSub$gradient[1], 2), " || larger unit (red / blue), n=", nrow(largeSub)," : ", largeSub$locBySerial[1], " - gradient: ", round(largeSub$gradient[1], 2)," || normalised p = ", round(prop$p.value, 4))
    
    if(metric == "TTR") {
      plot(cumsum(smallSub$TTR_cc_score), main = mainText, xlim = (c(0, max(nrow(smallSub), nrow(largeSub)))), ylim = c(min(cumsum(smallSub$TTR_cc_score), cumsum(largeSub$TTR_cc_score)), max(cumsum(smallSub$TTR_cc_score), cumsum(largeSub$TTR_cc_score))))
      points(cumsum(largeSub$TTR_cc_score), col = "red")
    }
    
    if(metric == "resolve") {
      plot(cumsum(smallSub$hypoResolve_cc), main = mainText, xlim = (c(0, max(nrow(smallSub), nrow(largeSub)))), ylim = c(min(cumsum(smallSub$hypoResolve_cc), cumsum(largeSub$hypoResolve_cc)), max(cumsum(smallSub$hypoResolve_cc), cumsum(largeSub$hypoResolve_cc))))
      points(cumsum(largeSub$hypoResolve_cc), col = "red")
    }
    
    averageM = as.data.frame(matrix(nrow = nrow(largeSub), ncol = nrow(smallSub)))
    average_median = rep(0, nrow(smallSub))
    pval_vector = rep(0, nrow(largeSub))
    
    for (kk in seq(1, nrow(largeSub), 1)) {
      print(kk)
      if(metric == "TTR") {
        x = sample(largeSub$TTR_cc_score, nrow(smallSub))
      }
      if(metric == "resolve") {
        x = sample(largeSub$hypoResolve_cc, nrow(smallSub))
      }
      averageM[kk, ] = cumsum(x)
      # points(cumsum(x), col = rgb(0, 1, 0, opacity, maxColorValue = 1), pch = 16, cex = 0.2)
      prop_k <- prop.test(c(small_pos, sum(x==1)), c(nrow(smallSub), length(x)))
      pval_vector[kk] = prop_k$p.value
    }
    
    for (am in seq(1, length(average_median), 1)) {
      average_median[am] <- quantile(averageM[, am])[3]
    }
    
    points(average_median, col = rgb(0, 0, 1, 1, maxColorValue = 1), cex = 0)
    lines(average_median, col = rgb(0, 0, 1, 1, maxColorValue = 1), lwd = 4)
    
    print(prop$p.value)
    print(paste0("median pval, ", nrow(smallSub)," samples of larger set"))
    print(quantile(pval_vector)[3])
    
    print(paste0("harmonic mean multiple tests pval, ", nrow(smallSub)," samples of larger set"))
    print(p.hmp(pval_vector, L = length(pval_vector)))
    
    legendText = paste0('pval: props uncorrected: ', prop$p.value, "\npval: median ", nrow(smallSub), " sample: ", quantile(pval_vector)[3], "\npval: harmonic mean ", nrow(smallSub), " samples: ", p.hmp(pval_vector, L = length(pval_vector)))
    legend('topright', legend = legendText, bty = 'n')
    
    legend2Text = paste0('n input 1 / pre: ', nrow(TTR_sub_1), '\nn input 2 / post: ', nrow(TTR_sub_2))
    legend('bottomleft', legend = legend2Text, bty = 'n')
    
    dev.off()
    
    # if (timeBins == 1) {
    #   
    #   i1 <- unique(input_1$Operator.ID)
    #   i2 <- unique(input_2$Operator.ID)
    #   
    #   imatch = match(i2, i1)
    #   imatch[is.na(imatch)] <- 0
    #   imatch_ones <- ifelse(imatch > 0, 1, 0)
    #   
    #   total2 <- length(imatch_ones)
    #   retained <- sum(imatch_ones)
    #   new = total2 - retained
    #   
    #   churnF = new / total2
    #   
    #   if (metric == "TTR") {
    #     grad_1 <- cumsum(input_1$TTR_cc_score)[nrow(input_1)] / nrow(input_1[TTR_opportunity == 1])
    #     grad_2 <- cumsum(input_2$TTR_cc_score)[nrow(input_2)] / nrow(input_2[TTR_opportunity == 1])
    #   }
    #   
    #   return(list(grad_1, grad_2, churnF))
    #   
    # }
    
  }
  
  return(list(TTR_sub_1, TTR_sub_2))
  
}

## new functions

# select random 10 period from each admission
n_days_select <- function(dateTime, min_threshold, interval_seconds, uID, admission_vec) {
  
  print(uID)
  print(admission_vec)
  
  interval_days = interval_seconds / (60*60*24)
  
  window_flag <- rep(0, length(dateTime))
  
  if(interval_days[1] > min_threshold) {
    days <- as.Date(dateTime)
    days <- days[order(days)]
    last_start_date <- max(days) - min_threshold
    last_date_index <- which(days < last_start_date)[length(which(days < last_start_date))]
    random_start <- sample(last_date_index, 1)
    start_window <- days[random_start]
    end_window <- start_window + min_threshold
    window_flag <- ifelse((days >= start_window) & (days <= end_window), 1, 0)
  } else {
    window_flag <- rep(1, length(dateTime))
  }
  
  return(window_flag)
  
}
