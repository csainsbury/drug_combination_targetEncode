# intervention wards are:-
  # 5A,5C,7A, 8B,9A,10A,11B,
  # Langlands 54, 56
interventionVector <- c("5A" ,"5C" ,"7A" , "8B" ,"9A" ,"10A" ,"11B" , "54", "56")

# Control wards are:-
  # 5B,5D,7B,8C,9B,10B,11C
  # Langlands 55,57
controlVector <- c("5B" ,"5D" ,"7B" ,"8C" ,"9B" ,"10B" ,"11C" , "55" ,"57")

library(berryFunctions)

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
  
    }
  
  if (nrow(hypoFrame) == 1 & hypoFrame$glu[1] < hypoThresh) {
    n_hypo_episodes = 1
    n_hypo_episodes_resolved = 0
  }
  
  if (nrow(hypoFrame) == 1 & hypoFrame$glu[1] >= hypoThresh) {
    n_hypo_episodes = 0
    n_hypo_episodes_resolved = 0
  } 
  
  output <- list(n_hypo_episodes, n_hypo_episodes_resolved)
  
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
  
  output <- list(hypoFrame$flag_opportunity, hypoFrame$cc_score)
  
  return(output)
  
  }


# load meter serial data
meterSerial <- read.csv("~/projects/sanofi_inpatient/nonCBGdata/meterSerial.csv", header = FALSE, row.names = NULL, stringsAsFactors = FALSE)
colnames(meterSerial) <- c("ward", "serial")

# load data files and concatenate
file_list <- list.files(path="~/projects/sanofi_inpatient/CBGdata/")

for (i in seq(1, length(file_list), 1)) {
  print(i)
  name <- paste("~/projects/sanofi_inpatient/CBGdata/", file_list[i], sep = "")
  dataFile <- read.csv(name, header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
  if (i == 1) {data <- dataFile}
  if (i > 1) {data <- rbind(data, dataFile)}
}

library(data.table)
dataDT <- data.table(data)
dataDT <- unique(dataDT)
dataDT <- dataDT[substr(Patient.ID, 2, 2) != "/"]
dataDT <- dataDT[substr(Patient.ID, 1, 2) != "00"]

dupList <- c("22", "33", "44", "55", "66", "77", "88", "99")
for (j in c(1:length(dupList))) {
  print(j)
  dataDT <- dataDT[substr(Patient.ID, 3, 4) != dupList[j]]
}

dataDT$Patient.ID <- as.numeric(dataDT$Patient.ID)
dataDT <- dataDT[!is.na(dataDT$Patient.ID)]

dataDT$dateplustime1 <- returnUnixDateTime(dataDT$Test.Date.Time)
dataDT <- dataDT[order(Patient.ID, dateplustime1), ]

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
dataDT$isHypo <- ifelse(dataDT$GLUn < 4, 1, 0)
dataDT[, c("nHypoCBGs") := sum(isHypo) , by=.(Patient.ID)]
dataDT[, c("anyHypoCBGs_perID") := ifelse(sum(isHypo)>0, 1, 0) , by=.(Patient.ID)]
dataDT[, c("totalNCBGs") := max(CBGn) , by=.(Patient.ID)]
dataDT[, c("CBGn_by_locBySerial") := seq(1, .N, 1) , by=.(locBySerial)]
dataDT[, c("maxCBGn_by_locBySerial") := max(CBGn_by_locBySerial), by=.(locBySerial)]
dataDT[, c("isHypo_byLocation") := sum(isHypo), by=.(locBySerial)]

# flag intervention (2) vs control (1) vs neither (0)
dataDT$unitStudyFlag <- 0
for (m in seq(1, length(interventionVector), 1)) {
  print(m)
  interestUnit <- interventionVector[m]
  dataDT$unitStudyFlag <- ifelse(dataDT$locBySerial == interestUnit, 2, dataDT$unitStudyFlag)
}
for (n in seq(1, length(controlVector), 1)) {
  print(n)
  interestUnit <- controlVector[n]
  dataDT$unitStudyFlag <- ifelse(dataDT$locBySerial == interestUnit, 1, dataDT$unitStudyFlag)
}

# replace NAs in location.by.serial by ward idents from location data (taken from meterSerial file)
na_location_set <- dataDT[is.na(locBySerial)]
dataDT_noNAs <- dataDT[!is.na(locBySerial)]

locsOfInterest <- as.character(as.data.frame(table(meterSerial$ward))$Var1)

na_location_set[, c("locBySerial") := locationReturn(Location) , by=.(Instrument.Serial.No.)]

coerced_location_set <- na_location_set[locBySerial != "no match"]

dataDT <- rbind(dataDT_noNAs, coerced_location_set)
dataDT <- dataDT[order(Patient.ID, dateplustime1), ]

# add n admssions per ID col
dataDT[, c("nAdmssions_per_ID") := admission_N(dateplustime1, 5), by=.(Patient.ID)]
# add vector of admission number per ID
dataDT[, c("admissionN_vector") := admission_N_vector(dateplustime1, 5), by=.(Patient.ID)]
# numberCBGs by admission
dataDT[, c("CBGn_by_admission") := seq(1, .N, 1), by =. (Patient.ID, admissionN_vector)]
dataDT[, c("max_CBGn_by_admission") := max(CBGn_by_admission), by =. (Patient.ID, admissionN_vector)]
dataDT[, c("CBG_datetime_to_discharge") := max(dateplustime1) - dateplustime1, by =. (Patient.ID, admissionN_vector)]
dataDT[, c("admissionDuration") := (max(dateplustime1) - min(dateplustime1)) / (60*60*24), by =. (Patient.ID, admissionN_vector)]
dataDT[, c("nHypos_by_admission") := sum(isHypo), by =. (Patient.ID, admissionN_vector)]

dataDT[, c("TTR") := admissionTTR(dateplustime1, GLUn), by =. (Patient.ID, admissionN_vector)]
dataDT[, c("hypoEps", "resolvedHypoEps") := admission_hypoEpisodes(dateplustime1, GLUn, 4, 60, Patient.ID), by =. (Patient.ID, admissionN_vector)]

dataDT[, c("TTR_opportunity", "TTR_cc_score") := admissionTTR_controlChart(dateplustime1, GLUn, 4, 30), by =. (Patient.ID, admissionN_vector)]


# plot cc for TTR - whole system
TTR_sub <- dataDT[TTR_opportunity == 1]
TTR_sub <- TTR_sub[order(TTR_sub$dateplustime1), ]
TTR_sub$cc <- cumsum(TTR_sub$TTR_cc_score)
plot(TTR_sub$dateplustime1 / (60*60*24), TTR_sub$cc, main = "control chart for QEUH - success in repeating CBG within 15 mins post CBG < 4")

# plot for each unit / operator
# type = "unit"
# type = "operator"
# 
type = "studyFlag"
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
  
  if (type == "unit") {
    xlim_lower = -50
    xlim_upper = 600
  }
  if (type == "operator") {
    xlim_lower = -100
    xlim_upper = 100
  }
  if (type == "studyFlag") {
    xlim_lower = -500
    xlim_upper = 0
  }
  
  if (j == 1) {
    plot(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, xlim = c(min(dataDT$dateplustime1), max(dataDT$dateplustime1)), ylim = c(xlim_lower, xlim_upper), cex = 0, main = "control chart for individual clinical operators - success in repeating CBG within 15 mins post CBG < 4")
    lines(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc)
  }
  
  if (j > 1) {
    points(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, cex = 0)
    lines(TTR_unit_sub$dateplustime1, TTR_unit_sub$cc, col = ifelse(TTR_unit_sub$locBySerial == "5B" | TTR_unit_sub$locBySerial == "5A" | TTR_unit_sub$locBySerial == "ARU2", "red", "black"))
  }
  
}

# set up per admission per location dataset
peradmission_set <- dataDT[CBGn_by_admission == 1]
peradmission_set[, c("admissionN_perLocation") := seq(1, .N, 1), by =. (locBySerial)]
peradmission_set[, c("hypoEps_perLocation") := sum(hypoEps), by =. (locBySerial)]
peradmission_set[, c("resolvedHypoEps_perLocation") := sum(resolvedHypoEps), by =. (locBySerial)]

# funnel plot hypo eps vs resolved eps per location
funnelPlot(peradmission_set[admissionN_perLocation == 1]$resolvedHypoEps_perLocation, peradmission_set[admissionN_perLocation == 1]$hypoEps_perLocation, labels = peradmission_set[admissionN_perLocation == 1]$locBySerial)

funnelPlot(peradmission_set[admissionN_perLocation == 1 & locBySerial != "IAU"]$resolvedHypoEps_perLocation, peradmission_set[admissionN_perLocation == 1 & locBySerial != "IAU"]$hypoEps_perLocation, labels = peradmission_set[admissionN_perLocation == 1 & locBySerial != "IAU"]$locBySerial, main = "funnel plot for proportions: resolved hypo episodes / all hypo episodes vs n all hypo episodes")

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

