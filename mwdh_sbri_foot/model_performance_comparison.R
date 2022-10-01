
logit_progression <-c(0.6230937, 0.5969812, 0.6293591, 0.6200913, 0.6510151, 0.6018030, 0.6160818, 0.6056275,
                      0.6418712, 0.6499395)
xgb_progression <-c(0.6333029, 0.5887949, 0.6254172, 0.6156177, 0.6601937, 0.6113001, 0.6175233, 0.6143655, 0.6263067, 0.6534490)
xgb_extended_progression <- c(0.7288392, 0.7386518, 0.7146532, 0.7209416, 0.7208915, 0.7296956, 0.7252198, 0.7317549, 0.7475230, 0.7556619)
logit_extended_progression <- c(0.7126809, 0.7310789, 0.7111960, 0.7166578, 0.6937824, 0.7215373, 0.7282594, 0.7118849, 0.7522872, 0.7402830)


logit_extended_amputation <- c(0.8779524, 0.8947624, 0.8954126, 0.8638813)
xgb_extended_amputation <- c(0.8865769, 0.9097616, 0.8609386, 0.9201181)
logit_amputation <- c(0.8190521, 0.7622871, 0.7130072, 0.8319407)
xgb_amputation <- c(0.8635091, 0.7605276, 0.7088735, 0.8170678)

boxplot(logit_extended_amputation, xgb_extended_amputation)
boxplot(logit_extended_progression, xgb_extended_progression)

boxplot(logit_progression, xgb_progression)


ID <- c(1,1,1,1,2,2,2,2,3,3,3)
sex <- c(rep("F", 8), rep("M", 3))
screen_time_since_diagnosis <- c(48, 400, 760, 900, 20, 360, 800, 1400, 0, 400, 600)
outcome <- c(0,0,0,1,0,0,0,0,0,0,1)

diabetes <- as.data.frame(cbind(ID, sex, screen_time_since_diagnosis, outcome), 
                          stringsAsFactors = T)

model <- glm(outcome~ sex + screen_time_since_diagnosis,
             data = diabetes, 
             family = binomial(link = "cloglog"))
