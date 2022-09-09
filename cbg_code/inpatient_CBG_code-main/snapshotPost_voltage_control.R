# intervention wards are:-
  # 5A,5C,7A, 8B,9A,10A,11B,
  # Langlands 54, 56
interventionVector <- c("5A" ,"5C" ,"7A" , "8B" ,"9A" ,"10A" ,"11B" , "54", "56")

# Control wards are:-
  # 5B,5D,7B,8C,9B,10B,11C
  # Langlands 55,57
#
controlVector <- c("5B" ,"5D" ,"7B" ,"8C" ,"9B" ,"10B" ,"11C" , "55" ,"57")
# controlVector <- c("obsGyn")

TTR_assessment_time = 20
hypo_threshold = 3
resolve_window = 15 # (mins)
# folder <- "CBG_data_combined"
folder <- "controlData"


library(berryFunctions)
library(harmonicmeanp)

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
admission_N_vector <- function(datetime, admissionThresholdDays) {
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

# load meter serial data
meterSerial <- read.csv("~/projects/sanofi_inpatient/nonCBGdata/meterSerial.csv", header = FALSE, row.names = NULL, stringsAsFactors = FALSE)
colnames(meterSerial) <- c("ward", "serial")

# load data files and concatenate
file_list <- list.files(path=paste0("~/projects/sanofi_inpatient/", folder, "/"))

for (i in seq(1, length(file_list), 1)) {
  print(i)
  name <- paste("~/projects/sanofi_inpatient/", folder, "/", file_list[i], sep = "")
  dataFile <- read.csv(name, header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
  if (i == 1) {data <- dataFile}
  if (i > 1) {data <- rbind(data, dataFile)}
}

library(data.table)
dataDT <- data.table(data)
dataDT <- unique(dataDT)
# dataDT <- dataDT[substr(Patient.ID, 2, 2) != "/"]
# dataDT <- dataDT[substr(Patient.ID, 1, 2) != "00"]

# dupList <- c("22", "33", "44", "55", "66", "77", "88", "99")
# for (j in c(1:length(dupList))) {
#   print(j)
#   dataDT <- dataDT[substr(Patient.ID, 3, 4) != dupList[j]]
# }

if (folder != "controlData") {
dataDT$Patient.ID <- as.numeric(dataDT$Patient.ID)
dataDT <- dataDT[!is.na(dataDT$Patient.ID)]
}

dataDT$dateplustime1 <- returnUnixDateTime(dataDT$Test.Date.Time)
dataDT <- dataDT[order(Patient.ID, dateplustime1), ]

# cut all data from pre 1st June 2018
# dataDT <- dataDT[dateplustime1 > returnUnixDateTime("01/06/2018 12:00AM")]

# remove duplicates
# generate unique identifier comprising product of ID and datetimestamp
dataDT$uniqueID = dataDT$Patient.ID * dataDT$dateplustime1
if (folder == "controlData") {
  dataDT$sequentialID <- c(1:nrow(dataDT))
  dataDT$uniqueID <- dataDT$sequentialID * dataDT$dateplustime1
}
dataDT = dataDT[order(dataDT$uniqueID), ]
dataDT[, c("uniqueID_n") := seq(1, .N, 1) , by=.(uniqueID)]
print(dim(dataDT))
dataDT = dataDT[uniqueID_n == 1]
print(dim(dataDT))

dataDT[, c("CBGn") := seq(1, .N, 1) , by=.(Patient.ID)]
dataDT[, c("maxCBGn_byID") := max(CBGn) , by=.(Patient.ID)]
dataDT$serial <- paste(substr(dataDT$Instrument.Serial.No., 1, 7), substr(dataDT$Instrument.Serial.No., 9, 13), sep = "")
dataDT$locBySerial <- meterSerial$ward[match(dataDT$serial, meterSerial$serial)]
dataDT <- dataDT[GLU != ""]
dataDT$GLUn <- dataDT$GLU
dataDT$GLUn <- ifelse(dataDT$GLUn == ">27.8", 27.9, dataDT$GLUn)
dataDT$GLUn <- ifelse(dataDT$GLUn == "<1.1", 1.0, dataDT$GLUn)
dataDT$GLUn <- as.numeric(dataDT$GLUn)
dataDT[, c("medianGlu") := quantile(GLUn)[3] , by=.(Patient.ID)]
dataDT[, c("admissionN") := admission_N(dateplustime1, 5), by =. (Patient.ID)]

dataDT[, c("admissionDuration") := (max(dateplustime1) - min(dateplustime1)) / (60*60*24) , by=.(Patient.ID)]
dataDT$isHypo <- ifelse(dataDT$GLUn < hypo_threshold, 1, 0)
dataDT[, c("nHypoCBGs") := sum(isHypo) , by=.(Patient.ID)]
dataDT[, c("anyHypoCBGs_perID") := ifelse(sum(isHypo)>0, 1, 0) , by=.(Patient.ID)]
dataDT[, c("totalNCBGs") := max(CBGn) , by=.(Patient.ID)]
dataDT[, c("CBGn_by_locBySerial") := seq(1, .N, 1) , by=.(locBySerial)]
dataDT[, c("maxCBGn_by_locBySerial") := max(CBGn_by_locBySerial), by=.(locBySerial)]
dataDT[, c("isHypo_byLocation") := sum(isHypo), by=.(locBySerial)]

# replace NAs in location.by.serial by ward idents from location data (taken from meterSerial file)
# na_location_set <- dataDT[is.na(locBySerial)]
# dataDT_noNAs <- dataDT[!is.na(locBySerial)]
# 
# locsOfInterest <- as.character(as.data.frame(table(meterSerial$ward))$Var1)
# 
# na_location_set[, c("locBySerial") := locationReturn(Location) , by=.(Instrument.Serial.No.)]
# 
# coerced_location_set <- na_location_set[locBySerial != "no match"]
# 
# dataDT <- rbind(dataDT_noNAs, coerced_location_set)
# dataDT <- dataDT[order(Patient.ID, dateplustime1), ]
# 
#     # flag intervention (2) vs control (1) vs neither (0)
#     dataDT$unitStudyFlag <- 0
#     for (m in seq(1, length(interventionVector), 1)) {
#       print(m)
#       interestUnit <- interventionVector[m]
#       dataDT$unitStudyFlag <- ifelse(dataDT$locBySerial == interestUnit, 2, dataDT$unitStudyFlag)
#     }
#     for (n in seq(1, length(controlVector), 1)) {
#       print(n)
#       interestUnit <- controlVector[n]
#       dataDT$unitStudyFlag <- ifelse(dataDT$locBySerial == interestUnit, 1, dataDT$unitStudyFlag)
#       }

# # add n admssions per ID col
# dataDT[, c("nAdmssions_per_ID") := admission_N(dateplustime1, 5), by=.(Patient.ID)]
# # add vector of admission number per ID
# dataDT[, c("admissionN_vector") := admission_N_vector(dateplustime1, 5), by=.(Patient.ID)]
# # numberCBGs by admission
# dataDT[, c("CBGn_by_admission") := seq(1, .N, 1), by =. (Patient.ID, admissionN_vector)]
# dataDT[, c("max_CBGn_by_admission") := max(CBGn_by_admission), by =. (Patient.ID, admissionN_vector)]
# dataDT[, c("CBG_datetime_to_discharge") := max(dateplustime1) - dateplustime1, by =. (Patient.ID, admissionN_vector)]
# dataDT[, c("admissionDuration") := (max(dateplustime1) - min(dateplustime1)) / (60*60*24), by =. (Patient.ID, admissionN_vector)]
# dataDT[, c("nHypos_by_admission") := sum(isHypo), by =. (Patient.ID, admissionN_vector)]
# 
# dataDT[, c("TTR") := admissionTTR(dateplustime1, GLUn), by =. (Patient.ID, admissionN_vector)]
# dataDT[, c("endOfHypoEpisode", "resolvedHypoEpisode", "hypoResolve_cc", "hypoEps", "resolvedHypoEps") := admission_hypoEpisodes(dateplustime1, GLUn, hypo_threshold, resolve_window, Patient.ID), by =. (Patient.ID, admissionN_vector)]
# 
# dataDT[, c("TTR_time_to_next_mins", "TTR_opportunity", "TTR_cc_score") := admissionTTR_controlChart(dateplustime1, GLUn, hypo_threshold, TTR_assessment_time), by =. (Patient.ID, admissionN_vector)]
# 
# dataDT[, c("meter_prior_measurement") := priorMeasurementTime(dateplustime1), by =. (serial)]
# dataDT[, c("location_prior_measurement") := priorMeasurementTime(dateplustime1), by =. (locBySerial)]

# 
# # plot cc for TTR - whole system
# TTR_sub <- dataDT[TTR_opportunity == 1]
# TTR_sub <- TTR_sub[order(TTR_sub$dateplustime1), ]
# TTR_sub$cc <- cumsum(TTR_sub$TTR_cc_score)
# plot(TTR_sub$dateplustime1 / (60*60*24), TTR_sub$cc, main = paste0("control chart for QEUH - success in repeating CBG within ", TTR_assessment_time, " mins post CBG < ", hypo_threshold))
# abline(v = returnUnixDateTime("01/03/2019 12:00AM") / (60*60*24))
# 
# 
#     unitVec <- c("1", "2")
#     for (u in c(1:2)) {
#       
#           ## pre/post - TTR
#           interestUnit = unitVec[u]
#           preN = nrow(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & TTR_opportunity == 1])
#           print("pre n = "); print(nrow(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & TTR_opportunity == 1]))
#           x = normalisedControlGradient(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")],
#                                         dataDT[unitStudyFlag == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")],
#                                         1, "TTR", "studyFlag", paste0("prePost_preN_", preN,"_", TTR_assessment_time, ifelse(interestUnit == "1", "control", ifelse(interestUnit == "2", "intervention", "other"))), 0.02)
#           
#           ## pre/post - resolve
#           preN = nrow(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & endOfHypoEpisode  == 1])
#           print("pre n = "); print(preN)
#           normalisedControlGradient(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")],
#                                     dataDT[unitStudyFlag == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")],
#                                     1, "resolve", "studyFlag", paste0("prePost_preN_", preN,"_", TTR_assessment_time, ifelse(interestUnit == "1", "_control", ifelse(interestUnit == "2", "_intervention", "_other"))), 0.02)
#       
#     }
#     
#     
#     ## check proportion of events pre/post
#     firstHalf <- dataDT[dateplustime1 > returnUnixDateTime("01/06/2018 12:00AM") & dateplustime1 <=returnUnixDateTime("01/03/2019 12:00AM")]
#     secondHalf <- dataDT[dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM") & dateplustime1 <=returnUnixDateTime("01/12/2019 12:00AM")]
#     
#     firstHalf_intervention <- dataDT[dateplustime1 > returnUnixDateTime("01/06/2018 12:00AM") & dateplustime1 <=returnUnixDateTime("01/03/2019 12:00AM") & locBySerial %in% interventionVector]
#     firstHalf_control <- dataDT[dateplustime1 > returnUnixDateTime("01/06/2018 12:00AM") & dateplustime1 <=returnUnixDateTime("01/03/2019 12:00AM") & locBySerial %in% controlVector]
#     secondHalf_intervention <- dataDT[dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM") & dateplustime1 <=returnUnixDateTime("01/12/2019 12:00AM") & locBySerial %in% interventionVector]
#     secondHalf_control <- dataDT[dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM") & dateplustime1 <=returnUnixDateTime("01/12/2019 12:00AM") & locBySerial %in% controlVector]
#     
#     prop.test(c(nrow(firstHalf_intervention[GLUn < hypoT]), nrow(firstHalf_control[GLUn < hypoT])),
#               c(nrow(firstHalf_intervention), nrow(firstHalf_control)))
#     prop.test(c(nrow(secondHalf_intervention[GLUn < hypoT]), nrow(secondHalf_control[GLUn < hypoT])),
#               c(nrow(secondHalf_intervention), nrow(secondHalf_control)))
#     
#     prop.test(c(nrow(firstHalf_intervention[GLUn < hypoT]), nrow(secondHalf_intervention[GLUn < hypoT])),
#               c(nrow(firstHalf_intervention), nrow(secondHalf_intervention)))
#     prop.test(c(nrow(firstHalf_control[GLUn < hypoT]), nrow(secondHalf_control[GLUn < hypoT])),
#               c(nrow(firstHalf_control), nrow(secondHalf_control)))
    
    
  
    ## voltage and temperature relationships
    ###
    tempSetDT_low <- dataDT[GLUn >=1 &GLUn < 5]
    tempSetDT_high <- dataDT[GLUn > 10]
    
    
    tempSetDT <- tempSetDT_low
    tempSetDT <- tempSetDT_high
    tempSetDT$timeOnly <- substr(tempSetDT$Test.Date.Time, 12, 18)
    
    tempSetDT$fractionalDay <- ((as.numeric(as.POSIXct(paste("01/01/2000", tempSetDT$timeOnly), format="%d/%m/%Y %I:%M%p", tz="GMT"))) - 
            as.numeric(as.POSIXct("01/01/2000 12:00AM",  format="%d/%m/%Y %I:%M%p", tz="GMT"))) / (60*60*24)
    tempSetDT$roundedFractionalDay <- round(tempSetDT$fractionalDay, 4)
    
    # distribution of differing voltages across day (uncorrected for n of tests)
    hist(tempSetDT[BatteryVoltage > 3.4]$fractionalDay, breaks = 40)
    hist(tempSetDT[BatteryVoltage > 3.38 & BatteryVoltage < 3.4]$fractionalDay, breaks= 40)
    
    # match to randomly selected equal sized subset of measures at lower voltages
    highV <- tempSetDT[BatteryVoltage > 3.4]
    # lowerV <- tempSetDT[BatteryVoltage > 3.38 & BatteryVoltage <= 3.4]
    lowerV <- tempSetDT[BatteryVoltage <= 3.4]
    
    sampledLowerV <- lowerV[sample(nrow(lowerV), nrow(highV))]
    
    hist(highV$fractionalDay)
    hist(sampledLowerV$fractionalDay)
    
    wilcox.test(highV$fractionalDay, sampledLowerV$fractionalDay)
    
    # match each time with voltage > 3.4 with a time with voltage below this threshold
    
    pVal_vector <- rep(0, 10)
    matchedCBG_frame <- as.data.frame(matrix(nrow = nrow(highV), ncol = length(pVal_vector)))
    timeFromPrior_frame <- as.data.frame(matrix(nrow = nrow(highV), ncol = length(pVal_vector)))
    n_matchingPool <- as.data.frame(matrix(nrow = nrow(highV), ncol = length(pVal_vector)))
    
    findMatch <- function(glu_1_case, glu_2_pool, fractionalTime_1_case, fractionalTime_2_pool, timeWindow) {
    # findMatch <- function(glu_1_case, glu_2_pool, fractionalTime_1_case, fractionalTime_2_pool, priorTime_1_case, priorTime2_pool, timeWindow) {
      # glu_1_case = highV$GLUn[64]; glu_2_pool = lowerV$GLUn; fractionalTime_1_case = highV$roundedFractionalDay[64]; fractionalTime_2_pool = lowerV$roundedFractionalDay; timeWindow = 1
      
      poolDT <- data.table(glu_2_pool, fractionalTime_2_pool)
      matchPool <- poolDT[fractionalTime_1_case == fractionalTime_2_pool]
      output = matchPool[sample(nrow(matchPool), 1)]
      
      return(output$glu_2_pool)
      
    }
    
    for (tt in c(1:length(pVal_vector))) {
      
      print(tt)
      
      highV[, c("matched_glu") := findMatch(GLUn, lowerV$GLUn, roundedFractionalDay, lowerV$roundedFractionalDay, 2), by=.(uniqueID)]
      
      matchedCBG_frame[, tt] <- highV$matched_glu
      # timeFromPrior_frame[, tt] <- highV$matched_prior
      
      highV$matched_glu <- NULL # ; highV$matched_prior <- NULL
    }
    
    
    # # for (t in c(1:length(pVal_vector))) {
    #   print(t)
    #     lowerV$pool <- 1
    #     lowerV$roundedFractionalDay <- round(lowerV$fractionalDay, 3)
    #     highV$roundedFractionalDay <- round(highV$fractionalDay, 3)
    #     matchedFrame <- as.data.frame(matrix(nrow = nrow(highV), ncol = ncol(lowerV)))
    #     for (r in seq(1, nrow(highV), 1)) {
    #       if(r%%100 == 0) {print(r)}
    #       interestTime <- highV$roundedFractionalDay[r]
    #       interestPriorUse <- highV$meter_prior_measurement[r]
    #         timeWindow = interestPriorUse * 0.1
    #       matchPool <- lowerV[roundedFractionalDay == interestTime &
    #       #                      pool == 1 &
    #                             meter_prior_measurement >= (interestPriorUse - (timeWindow/2)) &
    #                             meter_prior_measurement <= (interestPriorUse + (timeWindow/2))]
    #       selectedMatch <- matchPool[sample(nrow(matchPool), 1)]
    #       # lowerV[dateplustime1 == selectedMatch$dateplustime1[1]]$pool <- 0
    #       matchedFrame[r, ] <- selectedMatch[1, ]
    #       n_matchingPool[r, t] <- nrow(matchPool)
    #     }
    #     
    #     colnames(matchedFrame) <- colnames(lowerV)
    #     matchedFrame <- data.table(matchedFrame)
    #     
    #     hist(highV$GLUn)
    #     hist(matchedFrame$GLUn)
    #     
    #     distTest <- wilcox.test(highV$GLUn, matchedFrame$GLUn)
    #     pVal_vector[t] <- distTest$p.value
    #     print(distTest$p.value)
    #     
    #     matchedCBG_frame[, t] <- matchedFrame$GLUn
    #     timeFromPrior_frame[, t] <- matchedFrame$meter_prior_measurement
    #     
    #     
    # }
    
    matchedCBG_frame$highV_values <- highV$GLUn
   # timeFromPrior_frame$highV_values <- highV$location_prior_measurement
   # write.csv(matchedCBG_frame, file = paste0("~./projects/sanofi_inpatient/nonCBGdata/matchedCBG_withOrginalhighVVals_", length(pVal_vector), "_matches.csv"))
    matchedCBG_DT <- data.table(matchedCBG_frame)
    
    # # harmonic mean of pvalues:
    # mean_p = p.hmp(pVal_vector, L = length(pVal_vector))
    # print(mean_p)
    # boxplot(matchedCBG_frame, legend = text(paste0('topright', 'harmonic mean p: ', mean_p, '/nn samples: ', length(pVal_vector))))
    
    # plot of distributions of target and matched frames
    for (jj in seq(1, ncol(matchedCBG_frame), 1)) {
      h = hist(matchedCBG_frame[, jj], plot = F, breaks = seq(0, 40, 0.1))
      end_col = ncol(matchedCBG_frame)
      if (jj == 1) {
        plot(h$counts, cex = 0, xlim = c(17, 40),  ylim = c(0, 1000))
        lines(h$counts, col = jj)
      }
      if (jj > 1 & jj < ncol(matchedCBG_frame)) {
        points(h$counts, cex = 0)
        lines(h$counts, col = jj)
      }
      if (jj == ncol(matchedCBG_frame)) {
        points(h$counts, cex = 0)
        lines(h$counts, col = "black", lwd = 4)
      }
    }
    
    ## matrix of wilcox tests
    wilcoxFrame <- as.data.frame(matrix(nrow = ncol(matchedCBG_frame), ncol = ncol(matchedCBG_frame)))
   # wilcoxFrame_time <- as.data.frame(matrix(nrow = ncol(matchedCBG_frame), ncol = ncol(matchedCBG_frame)))
    wilcoxFrame_hypoValues <- as.data.frame(matrix(nrow = ncol(matchedCBG_frame), ncol = ncol(matchedCBG_frame)))
    
    for (w in c(1:nrow(wilcoxFrame))) {
      print(w)
      for (v in c(1:ncol(wilcoxFrame))) {
        test <- wilcox.test(matchedCBG_frame[, w], matchedCBG_frame[, v])
        wilcoxFrame[w, v] <- test$p.value
        
       # test_time <- wilcox.test(timeFromPrior_frame[, w], timeFromPrior_frame[, v])
       # wilcoxFrame_time[w, v] <- test_time$p.value
        
        # hypo1 <- matchedCBG_frame[, w]
        # hypo1  = hypo1[hypo1 <4]
        # print(quantile(hypo1))
        # hypo2 <- matchedCBG_frame[, v]
        # hypo2  = hypo2[hypo2 <4]
        # hypoTest <- wilcox.test(hypo1, hypo2)
        # wilcoxFrame_hypoValues[w, v] <- hypoTest$p.value
      }
    }
    
    
harmonicP_vector <- rep(1, ncol(wilcoxFrame_hypoValues) - 1) 
medianP <- rep(1, ncol(wilcoxFrame_hypoValues) - 1)

for (w in seq(1, nrow(wilcoxFrame_hypoValues), 1)) {
  testVector <- wilcoxFrame_hypoValues[w,][-w]
  harmonicP_vector[w] <- p.hmp(testVector, L = length(testVector))
  medianP[w] <- as.numeric(quantile(testVector)[3])
}    

    write.csv(wilcoxFrame, file = paste0("~./projects/sanofi_inpatient/nonCBGdata/pval_multipleWilcoxTests_", length(pVal_vector), "_matches.csv"))
    write.csv(wilcoxFrame_hypoValues, file = paste0("~./projects/sanofi_inpatient/nonCBGdata/pval_multipleWilcoxTests_", length(pVal_vector), "_matches_CBGbelow4.csv"))
    
    
    x = boxplot(tempSetDT$GLUn ~ cut(tempSetDT$BatteryVoltage, breaks = seq(3.3, 3.42, 0.02)), varwidth = T)
    x = boxplot(tempSetDT$Temperature ~ cut(tempSetDT$BatteryVoltage, breaks = seq(3.3, 3.42, 0.02)), varwidth = T)
    
    
  
    
    x_lower = dataDT[BatteryVoltage > 3.38 & BatteryVoltage <= 3.4]
    x_high = dataDT[BatteryVoltage > 3.4 & BatteryVoltage <= 3.42]
    wilcox.test(log(x_lower$GLUn), log(x_high$GLUn))
    
    ## 
    x = boxplot(log(dataDT$GLUn) ~ cut(dataDT$Temperature, breaks = 100), varwidth = T)
    
    x_lower = dataDT[BatteryVoltage > 3.38 & BatteryVoltage <= 3.4]
    x_high = dataDT[BatteryVoltage > 3.4 & BatteryVoltage <= 3.42]
    wilcox.test(log(x_lower$GLUn), log(x_high$GLUn))
    
    x = boxplot(log(dataDT$BatteryVoltage) ~ cut(dataDT$Temperature, breaks = 100), varwidth = T)
    
    
    
    
    
    
    
    
    

#### generate control charts and save out outputs
controlChartFunction(dataDT, "TTR", "studyFlag", "studyFlag_TTR_45_all", -0, 1000, 0.4)
controlChartFunction(dataDT, "resolve", "studyFlag", "studyFlag_resolveAll_4_60_all", -1200, 200, 0.4)

controlChartFunction(dataDT, "TTR", "unit", "all_units_TTR_15", -1000, 100, 0.4)
controlChartFunction(dataDT, "resolve", "unit", "all_units_resolve_15", -800, 400, 0.4)

  ## tests for DTS abstract
    DTS_sub <- dataDT
    DTS_sub$unitStudyFlag <- 0
    DTS_sub$unitStudyFlag <- ifelse(DTS_sub$locBySerial == "11A", 1, 0)
    
    # plot control charts
    controlChartFunction(DTS_sub, "TTR", "studyFlag", "11A_TTR_11A_vs_all", -10596, 100, 0.4)
    controlChartFunction(DTS_sub, "resolve", "studyFlag", "11A_resolve_11A_vs_all", -10596, 100, 0.4)
    # units against each other
    normalisedControlGradient(DTS_sub[unitStudyFlag == 1], DTS_sub[unitStudyFlag == 0], 1, "TTR", "studyFlag", "AA__1", 0.2)
    



controlChartFunction(dataDT[locBySerial == "5B"], "resolve", "operator", "test5B", -100, 100, 0.4)
controlChartFunction(dataDT[locBySerial == "5B"], "TTR", "unit", "test5B_unit", -100, 100, 0.4)

controlChartFunction(dataDT[locBySerial == "11A"], "TTR", "operator", "11A_op", -100, 100, 0.4)
controlChartFunction(dataDT[locBySerial == "11D"], "TTR", "operator", "11D_op", -100, 100, 0.4)

controlChartFunction(dataDT[locBySerial == "11D"], "TTR", "operator", "11D_op", -75, 20, 0.4)
controlChartFunction(dataDT[locBySerial == "11D"], "TTR", "unit", "11D_unit", -1000, 100, 0.4)



controlChartFunction(dataDT[locBySerial == "11D"], "TTR", "operator", "test11D_oper2", -30, 42, 0.2)
controlChartFunction(dataDT[locBySerial == "11A"], "TTR", "operator", "test11A_oper2", -30, 42, 0.2)
controlChartFunction(dataDT[locBySerial == "4A"], "TTR", "operator", "test4A_oper2", -30, 42, 0.2)


controlChartFunction(dataDT[locBySerial == "obsGyn"], "TTR", "unit", "testobsGyn", -200, 0)
controlChartFunction(dataDT[locBySerial == "11D"], "TTR", "operator", "test11D_op", -30, 42)


# individual units against each other
normalisedControlGradient(dataDT[locBySerial == "5A"], dataDT[locBySerial == "5B"], 2, "TTR", "unit", "AAA_1", 0.2)
normalisedControlGradient(dataDT[locBySerial == "5A"], dataDT[locBySerial == "5B"], 2, "resolve", "unit", "AAA_1", 0.2)


normalisedControlGradient(dataDT[unitStudyFlag == 2], dataDT[unitStudyFlag == 1], 2, "TTR", "studyFlag", "test_TTR_30_intervention_vs_control", 0.2)

## pre/post
interestUnit = "5B"
normalisedControlGradient(dataDT[locBySerial == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")],
                          dataDT[locBySerial == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")],
                          1, "resolve", "unit", paste0("prePost_", interestUnit, "_TTRthresh_", TTR_assessment_time), 0.2)

    unitVec <- c("1", "2")
    for (u in c(1:2)) {
      
        ## pre/post - TTR
        interestUnit = unitVec[u]
        preN = nrow(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & TTR_opportunity == 1])
        print("pre n = "); print(nrow(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & TTR_opportunity == 1]))
        x = normalisedControlGradient(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")],
                                  dataDT[unitStudyFlag == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")],
                                  1, "TTR", "studyFlag", paste0("prePost_preN_", preN,"_", TTR_assessment_time, ifelse(interestUnit == "1", "control", ifelse(interestUnit == "2", "intervention", "other"))), 0.02)
        
        ## pre/post - resolve
        preN = nrow(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & endOfHypoEpisode  == 1])
        print("pre n = "); print(preN)
        normalisedControlGradient(dataDT[unitStudyFlag == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")],
                                  dataDT[unitStudyFlag == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")],
                                  1, "resolve", "studyFlag", paste0("prePost_preN_", preN,"_", TTR_assessment_time, ifelse(interestUnit == "1", "_control", ifelse(interestUnit == "2", "_intervention", "_other"))), 0.02)
    
    }


# all locs
locTable <- as.data.frame(table(dataDT$locBySerial))
locTable$grad1 <- 0
locTable$grad2 <- 0
locTable$churn <- 0

for (l in seq(1, nrow(locTable), 1)) {
  
  interestUnit = as.character(locTable$Var1[l])
  
  preN = nrow(dataDT[locBySerial == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & TTR_opportunity == 1])
  
  if (preN > 80) {
  print("pre n = "); print(nrow(dataDT[locBySerial == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM") & TTR_opportunity == 1]))
  x = normalisedControlGradient(dataDT[locBySerial == interestUnit & dateplustime1 <= returnUnixDateTime("01/03/2019 12:00AM")],
                                dataDT[locBySerial == interestUnit & dateplustime1 > returnUnixDateTime("01/03/2019 12:00AM")],
                                1, "TTR", "studyFlag", paste0("prePost_preN_", preN,"_", TTR_assessment_time, ifelse(interestUnit == "1", "control", ifelse(interestUnit == "2", "intervention", "other"))), 0.02)
  
  locTable$grad1[l] <- x[[1]]
  locTable$grad2[l] <- x[[2]]
  locTable$churn[l] <- x[[3]]
  
  }

}

locT <- data.table(locTable)
locT <- locT[grad1 != 0]

plot((locT$grad2 - locT$grad1), locT$churn)
cor.test((locT$grad2 - locT$grad1), locT$churn)


# plot for each unit / operator
#type = "unit"
#
type = "operator"
#type = "studyFlag"
TTR_sub <- dataDT[TTR_opportunity == 1]
if (type == "unit") {units <- as.data.frame(table(TTR_sub$locBySerial))}
if (type == "operator") {units <- as.data.frame(table(TTR_sub$Operator.ID))}
if (type == "studyFlag") {units <- as.data.frame(table(TTR_sub$unitStudyFlag))}


for (j in seq(1, nrow(units), 1)) {
  if (type == "unit") {TTR_unit_sub <- TTR_sub[locBySerial == units$Var1[j]]}
  if (type == "operator") {TTR_unit_sub <- TTR_sub[Operator.ID == units$Var1[j]]}
  if (type == "studyFlag") {TTR_unit_sub <- TTR_sub[unitStudyFlag == units$Var1[j]]}
  
  TTR_unit_sub <- TTR_unit_sub[order(TTR_unit_sub$dateplustime1), ]
  TTR_unit_sub$cc <- cumsum(TTR_unit_sub$TTR_cc_score)

  print(paste("unit: ", units$Var1[j], " total: ", TTR_unit_sub$cc[nrow(TTR_unit_sub)], sep = ""))
  
  
 # insertionRow <- c(as.character(units$Var1[j]), TTR_unit_sub$cc[nrow(TTR_unit_sub)], TTR_unit_sub$isHypo_byLocation[1])
  
#  rankingFrame = rbind(rankingFrame, insertionRow)
  
  if (type == "unit") {
    xlim_lower = -100
    xlim_upper = 650
  }
  if (type == "operator") {
    xlim_lower = -50
    xlim_upper = 50
  }
  if (type == "studyFlag") {
    xlim_lower = -100
    xlim_upper = 1000
  }
  
  nhypos = TTR_unit_sub$isHypo_byLocation[1]
  lineW = sqrt(nhypos) * 10e-2
  
  if (j == 1) {
    plot(TTR_unit_sub$dateplustime1,
         TTR_unit_sub$cc,
         xlim = c(min(dataDT$dateplustime1), max(dataDT$dateplustime1)),
         ylim = c(xlim_lower, xlim_upper),
         cex = 0,
         lwd = lineW,
         main = paste0("control chart for ", type," - success in repeating CBG within ", TTR_assessment_time, " mins post CBG < 4"))
    lines(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc)
    text(max(TTR_unit_sub$dateplustime1) + 1000000, TTR_unit_sub$cc[nrow(TTR_unit_sub)], labels = units$Var1[j], cex = 0.6)
  }
  
  if (j > 1) {
    points(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, cex = 0)
    lines(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, lwd = lineW, col = ifelse(TTR_unit_sub$locBySerial == "5B" | TTR_unit_sub$locBySerial == "5A", rgb(1, 0, 0, 0.5, maxColorValue = 1), ifelse(TTR_unit_sub$locBySerial == "ARU2", rgb(0, 1, 0, 0.5, maxColorValue = 1), rgb(0, 0, 0, 0.5, maxColorValue = 1))))
    text(max(TTR_unit_sub$dateplustime1) + 1000000, TTR_unit_sub$cc[nrow(TTR_unit_sub)], labels = units$Var1[j], cex = 0.6)
    abline(v = returnUnixDateTime("01/03/2019 12:00AM"))
    
  }
  
}

# gradient plots n bins - needs work. gradients not making sense
for (j in seq(1, nrow(units), 1)) {
  if (type == "unit") {TTR_unit_sub <- TTR_sub[locBySerial == units$Var1[j]]}
  if (type == "operator") {TTR_unit_sub <- TTR_sub[Operator.ID == units$Var1[j]]}
  if (type == "studyFlag") {TTR_unit_sub <- TTR_sub[unitStudyFlag == units$Var1[j]]}
  
  TTR_unit_sub <- TTR_unit_sub[order(TTR_unit_sub$dateplustime1), ]
  TTR_unit_sub$cc <- cumsum(TTR_unit_sub$TTR_cc_score)
  
  print(paste("unit: ", units$Var1[j], " total: ", TTR_unit_sub$cc[nrow(TTR_unit_sub)], sep = ""))
  
  n_bins = 2
  stepSize = (max(dataDT$dateplustime1) - min(dataDT$dateplustime1)) / n_bins
  checkPoints = seq(min(dataDT$dateplustime1), max(dataDT$dateplustime1), stepSize)
  #checkPoints = (checkPoints - min(dataDT$dateplustime1)) / (60*60*24)
  
  gradientFrame = as.data.frame(matrix(nrow = 1, ncol = n_bins))
  
  for (k in seq(1, n_bins, 1)) {
    print("step: "); print(k)
    if(k <= n_bins) {
      measurementSub <- TTR_unit_sub[dateplustime1 > checkPoints[k] & dateplustime1 <= checkPoints[k+1]]
    }
    
    
    subGradient <- (measurementSub$cc[nrow(measurementSub)] - measurementSub$cc[1]) / ((checkPoints[k+1] - checkPoints[k]) / (60*60*24))
    
    gradientFrame[k] = subGradient
    
  }
  
  col_val = gradientFrame[1, ncol(gradientFrame)] - gradientFrame[1, 1]
  print(col_val)
  

  if (type == "unit") {
    xlim_lower = -0.6
    xlim_upper = 0.6
  }
  if (type == "operator") {
    xlim_lower = -0.2
    xlim_upper = 0.2
  }
  if (type == "studyFlag") {
    xlim_lower = -5
    xlim_upper = 0
  }
  
  nhypos = TTR_unit_sub$isHypo_byLocation[1]
  lineW = sqrt(nhypos) * 10e-2
  
  if (j == 1) {
    plot(c(1: n_bins),
         gradientFrame[1, ],
         ylim = c(xlim_lower, xlim_upper),
         lwd = lineW,
         main = paste0("gradient chart for ", type," - success in repeating CBG within ", TTR_assessment_time, " mins post CBG < 4"))
    lines(c(1: n_bins), gradientFrame[1, ], col = ifelse(col_val < 0, rgb(1, 0, 0, 0.4, maxColorValue = 1), rgb(0, 1, 0, 0.4, maxColorValue = 1)))
    # text(max(TTR_unit_sub$dateplustime1) + 1000000, TTR_unit_sub$cc[nrow(TTR_unit_sub)], labels = units$Var1[j], cex = 0.6)
  }
  
  if (j > 1) {
    points(c(1: n_bins), gradientFrame[1, ])
    lines(c(1: n_bins), gradientFrame[1, ], lwd = lineW, col = ifelse(col_val < 0, rgb(1, 0, 0, 0.4, maxColorValue = 1), rgb(0, 1, 0, 0.4, maxColorValue = 1)))
    text(n_bins + 0.02, gradientFrame[n_bins], labels = units$Var1[j], cex = 0.5)
    text( (1-0.02), gradientFrame[1], labels = units$Var1[j], cex = 0.5)
    abline(v = returnUnixDateTime("01/03/2019 12:00AM"))
    
  }
  
}


# set up per admission per location dataset
peradmission_set <- funnelPlotPrep(dataDT)
peradmission_set_period1 <- funnelPlotPrep(dataDT[dateplustime1 < returnUnixDateTime("01/03/2019 12:00AM")])
peradmission_set_period2 <- funnelPlotPrep(dataDT[dateplustime1 >= returnUnixDateTime("01/03/2019 12:00AM")])


# funnel plot hypo eps vs resolved eps per location
funnelPlot(peradmission_set[admissionN_perLocation == 1]$resolvedHypoEps_perLocation, peradmission_set[admissionN_perLocation == 1]$hypoEps_perLocation, labels = peradmission_set[admissionN_perLocation == 1]$locBySerial)

funnelPlot(peradmission_set[admissionN_perLocation == 1 & locBySerial != "IAU"]$resolvedHypoEps_perLocation, peradmission_set[admissionN_perLocation == 1 & locBySerial != "IAU"]$hypoEps_perLocation, labels = peradmission_set[admissionN_perLocation == 1 & locBySerial != "IAU"]$locBySerial, main = "funnel plot for proportions: resolved hypo episodes / all hypo episodes vs n all hypo episodes")

  # pre intervention
  # funnel plot hypo eps vs resolved eps per location
 x= funnelPlot(peradmission_set_period1[admissionN_perLocation == 1 & locBySerial != "IAU" & hypoEps_perLocation > 0]$resolvedHypoEps_perLocation,
             peradmission_set_period1[admissionN_perLocation == 1 & locBySerial != "IAU" & hypoEps_perLocation > 0]$hypoEps_perLocation,
             labels = peradmission_set_period1[admissionN_perLocation == 1 & locBySerial != "IAU" & hypoEps_perLocation > 0]$locBySerial,
             main = "funnel plot for proportions: resolved hypo episodes / all hypo episodes vs n all hypo episodes",
             ylim = c(0, 100),
             xlim = c(1, 700))
  
  # post intervention
  # funnel plot hypo eps vs resolved eps per location
  funnelPlot(peradmission_set_period2[admissionN_perLocation == 1 & locBySerial != "IAU" & hypoEps_perLocation > 0]$resolvedHypoEps_perLocation,
             peradmission_set_period2[admissionN_perLocation == 1 & locBySerial != "IAU" & hypoEps_perLocation > 0]$hypoEps_perLocation,
             labels = peradmission_set_period2[admissionN_perLocation == 1 & locBySerial != "IAU" & hypoEps_perLocation > 0]$locBySerial,
             main = "funnel plot for proportions: resolved hypo episodes / all hypo episodes vs n all hypo episodes",
             ylim = c(0, 100),
             xlim = c(1, 700))

## per operator
peradmission_set[, c("admissionN_perOpID") := seq(1, .N, 1), by =. (Operator.ID)]
peradmission_set[, c("hypoEps_perOpID") := sum(hypoEps), by =. (Operator.ID)]
peradmission_set[, c("resolvedHypoEps_perOpID") := sum(resolvedHypoEps), by =. (Operator.ID)]

funnelPlot(peradmission_set[admissionN_perOpID == 1 & hypoEps_perOpID > 0]$resolvedHypoEps_perOpID, peradmission_set[admissionN_perOpID == 1 & hypoEps_perOpID > 0]$hypoEps_perOpID)

## per ID
funnelPlot(peradmission_set[hypoEps > 0]$resolvedHypoEps, peradmission_set[hypoEps>0]$hypoEps, labels = peradmission_set[hypoEps > 0]$locBySerial)


# numbers
plot(table(dataDT[Department == "QEUH Wards"]$Location), las = 3, main = "number CBG measures jul-sep18 by indicated ward", ylab = "n")
plot(table(dataDT[Department == "QEUH Wards"]$locBySerial), las = 3, main = "number CBG measures jul-sep18 by ward via serial", ylab = "n")

plot(table(dataDT[CBGn == 1 & Department == "QEUH Wards"]$Location), las = 3, main = "unique IDs jul-sep18 by indicated ward", ylab = "n")
plot(table(dataDT[CBGn == 1 & Department == "QEUH Wards"]$locBySerial), las = 3, main = "unique IDs jul-sep18 by ward via serial", ylab = "n")

# average CBG by ward via serial
boxplot(dataDT[Department == "QEUH Wards"]$GLUn ~ dataDT[Department == "QEUH Wards"]$locBySerial, varwidth = T, las = 3, main = "boxplot measured glucose values per ward via serial")
boxplot(dataDT[CBGn == 1 & Department == "QEUH Wards"]$medianGlu ~ dataDT[CBGn == 1 & Department == "QEUH Wards"]$locBySerial, varwidth = T, las = 3, main = "boxplot median glucose values per ID by ward via serial")

# average admission duration by ward
boxplot(dataDT[admissionDuration > 0 & CBGn == 1 & Department == "QEUH Wards"]$admissionDuration ~ dataDT[admissionDuration > 0 & CBGn == 1 & Department == "QEUH Wards"]$locBySerial, varwidth = T, las = 3, ylim = c(0, 30),  main = "boxplot admission duration (days) per ID by ward via serial")

# n CBGs measured per day
boxplot((dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$totalNCBGs / dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$admissionDuration)
        ~ dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$locBySerial, varwidth = T, las = 3, ylim = c(0, 10), main = "boxplot number CBG measurements / day per ID by ward via serial")

# hypos per day
boxplot((dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$nHypoCBGs / dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$admissionDuration)
        ~ dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$locBySerial, varwidth = T, las = 3, ylim = c(0, 1), main = "boxplot measured hypo CBGs / day per ID by ward via serial")

# hypos per day vs total n tests

# hypos per ID
boxplot((dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$anyHypoCBGs_perID / dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$admissionDuration)
        ~ dataDT[admissionDuration > 1 & CBGn == 1 & Department == "QEUH Wards"]$locBySerial, varwidth = T, las = 3, ylim = c(0, 1), main = "boxplot measured hypo CBG per ID by ward via serial")

# rough funnel plot
loc <- dataDT[CBGn_by_locBySerial == 1]
plot(loc$maxCBGn_by_locBySerial, (loc$isHypo_byLocation / loc$maxCBGn_by_locBySerial))


hist(tempSetDT$BatteryVoltage, breaks = seq(3.29, 3.43, 0.001))

