##Read background risk from FMF table
##background <- read.csv2(file="backgroundrisk-FMF.csv",skip=1)
##Function to compute background risk based on age.
##background_risk <- approxfun( background$Age, background$W12)

##Read background risk from Snijders et al. (1999) table. This is a global variable.
background <- read.csv2(file="inst/extdata/snijders_etal1999_t21.csv")
devtools::use_data(background, internal = TRUE, overwrite=TRUE)
snijders_etal1999_t21 <- background
save(file="data/snijders_etal1999_t21.RData", list=c("snijders_etal1999_t21"))
