##Read background risk from FMF table
##background <- read.csv2(file="backgroundrisk-FMF.csv",skip=1)
##Function to compute background risk based on age.
##background_risk <- approxfun( background$Age, background$W12)

##Read background risk from Snijders et al. (1999) table. This is a global variable.
data(snijders_etal1999-t21)
background <- read.csv2(file="data/snijders_etal1999-t21.csv")
devtools::use_data(background, internal = TRUE)