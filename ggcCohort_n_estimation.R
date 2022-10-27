
xr <- read.csv('~/Documents/data/anonSD/demog_all2.txt')
x <- data.table(xr)

print(paste0('uniqueN_all: ', uniqueN(x$LinkId)))

x$DeathDate <- as.Date(x$DeathDate)
x$deathFlag <- ifelse(x$DeathDate < as.Date('2012-01-01'), 1, 0)
x$deathFlag[is.na(x$deathFlag)] <- 1
print(paste0('not dead / dead post 2012: ', uniqueN(x[deathFlag == 1]$LinkId)))
x <- x[deathFlag == 1]

x <- x[DiabetesMellitusType_Mapped == 'Type 2 Diabetes Mellitus']
print(paste0('uniqueN_DMT2: ', uniqueN(x$LinkId)))

# number of new diagnoses 2012-2017
new <- x[DateOfDiagnosisDiabetes_Date > as.Date('2012-01-01')]
paste0('number new diagnoses: ', uniqueN(new$LinkId))

##
eta 0.1 nrounds 200
0%        25%        50%        75%       100% 
0.01462878 0.07484655 0.10467344 0.13382010 0.19199714 
> 
  
  
eta 0.1, nrounds 400
0%         25%         50%         75%        100% 
0.003846773 0.069201936 0.102707370 0.139346426 0.214102457 

eta 0.01, nrounds 400
0%         25%         50%         75%        100% 
0.004715304 0.069946870 0.097903834 0.131786022 0.214793262 

eta 0.1, nrounds 800
0%        25%        50%        75%       100% 
0.01736599 0.07258174 0.10270737 0.13998417 0.21665343 
> 