#######################################################
## Functions of the trisomy21risk package
##
## Author: Michael Höhle <http://www.math.su.se/~hoehle>
## Date:   June 2018
## License: GNU-GPL v3.0
#######################################################

#' Manual least squares fit to the Merz et al. (2016) data using a 10th order polynomial 
#'
#' @param crl crown rump length (in mm)
#' @return The expected nuchal translucency (nt) from a 10th order polynomial fitted to the data of Merz et al.
#' @references 
#' Merz, E., C. Thode, B. Eiben, and S. Wellek. 2016. “Prenatal Risk Calculation (PRC) 3.0: An Extended Doe-Based First-Trimester Screening Algorithm Allowing for Early Blood Sampling.” Ultrasound International Open 2: E19–E26. doi:10.1055/s-0035-1569403.
#' 
#' @export
nt_ls <- function(crl) {
  (86.8495636412728*1) + (-6.45603285194204*crl) + (0.178676863518831*crl^2) + (-0.0019492362713131*crl^3) + (1.25983144232195e-07*crl^5) + (-6.84398153615137e-12*crl^7) + (1.53020537087493e-18*crl^10)
}

#' Manual least squares fit to the Merz et al. (2016) data using a 10th order polynomial 
#'
#' @param crlbe crown rump length (in mm) at the time of the serum sampling
#' @return The expected value of logbetahCG (from a 10th order polynomial fitted to the data of Merz et al.)
#' @references 
#' Merz, E., C. Thode, B. Eiben, and S. Wellek. 2016. “Prenatal Risk Calculation (PRC) 3.0: An Extended Doe-Based First-Trimester Screening Algorithm Allowing for Early Blood Sampling.” Ultrasound International Open 2: E19–E26. doi:10.1055/s-0035-1569403.
#' 
#' @export
logbetahCG <- function(crlbe) {
  5.05336 - 0.03595*crlbe +  2.18174e-4*crlbe^2 - 8.68749e-15*crlbe^7 - 2.36687e-21*crlbe^10
}

#' Manual least squares fit to the Merz et al. (2016) data using a 10th order polynomial 
#'
#' @param crlbe crown rump length (in mm) at the time of the serum sampling
#' @return The expected value of logPAPPA (from a 10th order polynomial fitted to the data of Merz et al.)
#' @references 
#' Merz, E., C. Thode, B. Eiben, and S. Wellek. 2016. “Prenatal Risk Calculation (PRC) 3.0: An Extended Doe-Based First-Trimester Screening Algorithm Allowing for Early Blood Sampling.” Ultrasound International Open 2: E19–E26. doi:10.1055/s-0035-1569403.
#' @export
logPAPPA <- function(crlbe) {
  -1.95000 + 0.07313*crlbe - 3.63693e-4* crlbe^2 + 1.39983e-20*crlbe^10
}

#' Compute the DoE (Degrees of Extremeness) of an observation
#' @param y The observed value 
#' @param y_pred This is the model based expected value
#' @param ref_lower Lower value of the reference band, i.e. the 5% quantile value 
#' @param ref_upper Upper value of the reference band, i.e. the 95% quantile value
#' @references 
#' Merz, E., C. Thode, B. Eiben, and S. Wellek. 2016. “Prenatal Risk Calculation (PRC) 3.0: An Extended Doe-Based First-Trimester Screening Algorithm Allowing for Early Blood Sampling.” Ultrasound International Open 2: E19–E26. doi:10.1055/s-0035-1569403.
#'
#' @export
doe <- function(y, y_pred, ref_lower, ref_upper) {
  ifelse(y>y_pred,
    (y - y_pred) / (ref_upper - y_pred),
    -( y_pred - y) / (y_pred - ref_lower)
  )
}

#' Inverse function of the DoE (Degrees of Extremeness), i.e. compute observed value given DoE
#'
#' @param y The observed value 
#' @param y_pred This is the model based expected value
#' @param ref_lower Lower value of the reference band, i.e. the 5% quantile value 
#' @param ref_upper Upper value of the reference band, i.e. the 95% quantile value
#' @seealso [trisomy21risk::doe]
#' @export
doe_inv <- function(doe, y_pred, ref_lower, ref_upper) {
  ifelse(doe >= 0,
         doe * (ref_upper - y_pred) + y_pred,
         doe * (y_pred - ref_lower) + y_pred
  )
}

#' Translate from age to CRL using the Wikipedia formula https://en.wikipedia.org/wiki/Crown-rump_length,
#' where CRL is measured in mm and gestational age in days.
#'
#'              Gestational age = (CRL x 1.04)^0.5 x 8.05 + 23.7
#'
#' Values on FMF page https://fetalmedicine.org/research/pregnancyDating with JavaScript https://fetalmedicine.org/assets/7a3311bf/js/base.min.js?v=51860
#' The page attributes https://www.ncbi.nlm.nih.gov/pubmed/1182090
#'
#' gaFromCrl:function(t,e,n){
#'   return e=e||45,n=n||84,t=this.parseFloat(t),null===t||t<e||t>n?null:23.53+8.052*Math.sqrt(1.037*t)
#' }
#'
#' We take the formula from http://journals.sagepub.com/doi/pdf/10.1179/174313409X448543
#' Gestational age = (CRL x 1.037)^0.5 x 8.052 + 23.73
#'
#' @param gestage Gestational age (in days)
#' @return Crown-rump-length (in mm)
#' @references http://journals.sagepub.com/doi/pdf/10.1179/174313409X448543
#' @export
gestage2crl <- function(gestage) {
  ##Wikipedia
  ##((gestage - 23.7)/8.05)^2 / 1.04
  ##FMF
  ##((gestage - 23.53)/8.052)^2 / 1.037
  ##Loughna et al. (2009)
  ((gestage - 23.73)/8.052)^2 / 1.037
}


#' Reverse function of gestage2crl 
#' 
#' @param gestage Gestational age (in days)
#' @return Crown-rump-length (in mm)
#' @seealso [trisomy21risk::gestage2crl]
#' @export
crl2gestage <- function(crl) {
  if (crl<30 | crl > 84) {
    warning("Conversion formula is only applicable for CRL values between 30 and 84 mm.")
    #return(NA)
  }
  ##Wikipedia
  ##(crl * 1.04)^0.5 * 8.05 + 23.7
  ##FMF
  ##(crl * 1.037)^0.5 * 8.052 + 23.53
  ##Loughna et al. (2009)
  (crl * 1.037)^0.5 * 8.052 + 23.73
  
}

#' Likelihood ratio between two univariate normal distributions
#'
#' @param x measured value
#' @param eu eucaryte (i.e. value in non-diseased population)
#' @param aneu anormal eucaryte (i.e. value in diseased population)
#' @return dnorm(x, mean=eu[1], sd=eu[2]) / dnorm(x, mean=aneu[1], sd=aneu[2])
#' @export
lr <- function(x, eu, aneu) {
  eu[2] / aneu[2] * exp(-0.5*((x-aneu[1])/aneu[2])^2 + 0.5*((x-eu[1])/eu[2])^2)
}

#' Compute posterior distribution for an event using the formula:
#' 
#'     posterior odds = likelihood ratio * prior odds
#'
#' @param lr_joint likelihood ratio
#' @param pi prior probability
#' @return The posterior odds as \eqn{1/(1+1/lr_joint*(1-pi)/pi)}
#' @export
ppost <- function(lr_joint, pi=1/100) {
  1/(1+1/lr_joint*(1-pi)/pi)
}


## Bandwidths, see Merz et al. (2016) paper
nt_qBands <-      c(-0.57159, 0.65869)
pappa_qBands <-   c(-0.93073, 0.87383)
betahCG_qBands <- c(-0.94259, 0.99854)

#' Compute background risk for T21 as a function of mother's age and gestational age.
#'
#' @param age_mother Age of the mother (in years)
#' @param gestational age (in days) at the time measurements were taken
#' @return Computes prior risk for Trisomy 21 based on the table in Snijders et al. (1999)
#' @export
background_risk <- function(age_mother, gestage) {
  ##Rounded week of gestation.
  gestweek_rounded <- max(10,min(14,floor(gestage / 7)))
  gestweek_exact <- gestage / 7
  
  ##Sanity checks
  if (age_mother < 20 | age_mother > 45) { warning("Age of mother out of bounds."); return(NA) }
  if (gestweek_rounded < 10 | gestweek_rounded > 14) { warning("Gestation age out of bounds."); return(NA) }
  
  ##If week 10,12,14 then extract background risk directly from snijders table
  if (gestweek_rounded %in% c(10,12,14)) {
    return(approxfun(background$Age, background[,paste0("W",gestweek_rounded)])(age_mother))
  }
  
  ##Otherwise interpolate between the two surrounding weeks
  if (gestweek_rounded %in% c(11,13)) {
    w_minus1 <- paste0("W",gestweek_rounded - 1)
    w_plus1 <- paste0("W",gestweek_rounded + 1)
    
    r_minus1  <- approxfun(background$Age, background[,w_minus1])(age_mother)
    r_plus1 <- approxfun(background$Age, background[,w_plus1])(age_mother)
    r <- approxfun(x=gestweek_rounded + c(-1,1), y=c(r_minus1,r_plus1))(gestweek_exact)
    return(r)
  }
  ##Otherwise
  return(NA)
}

#' Compute weight corrected logPAPPA as done in Merz et al. (2016)
#' 
#' @param pappa Concentration of PAPP-A in maternal serum (in miU/ml)
#' @param weight Weight of the mother (in kg)
#' @return log(PAPP-A) concentration corrected for the weight of the mother
#' @export
logPAPPA_corr   <- function(pappa, weight) {
  log(pappa) + 0.01654 * (weight - 68.341242)
}

#' Compute weight corrected logbetahCG as done in Merz et al. (2016)
#' 
#' @param betahCG Concentration of beta-hCG in maternal serum (in miU/ml)
#' @param weight Weight of the mother (in kg)
#' @return log(beta-hCG) concentration corrected for the weight of the mother
#' @export
logbetahCG_corr <- function(betahCG, weight) {
  log(betahCG) + 0.00972 * (weight - 68.341242)
}

######################################################################
#' Function to compute the trisomy21 risk based on 3 biomarkers (NT, PAPP-A, betahCG)
#' and visibility of the nasal bone on the ultrasound. Several covariates (age of mother, weight of mother, CRL at measurement)
#' are taken into account. 
#'  
#' @param age Age of the mother (in years)
#' @param weight Weight of the mother (in kg)
#' @param crlbe Crown-rump-length (in mm) at the time of the PAPP-A and betahCG measurements
#' @param crl Crown-rump-length (in mm) at the time of the NT measurement
#' @param pappa Measured value of PAPP-A in maternal serum (in miU/ml)
#' @param betahCG Measured value of betahCG in maternal serum (in miU/ml)
#' @param nasalbone Visibility of a nasal bone on the ultrasound (Boolean)
#' @param background Background probability for T21
#' @return Posterior risk of T21.
#' @references  
#' Merz, E., C. Thode, B. Eiben, and S. Wellek. 2016. “Prenatal Risk Calculation (PRC) 3.0: An Extended Doe-Based First-Trimester Screening Algorithm Allowing for Early Blood Sampling.” Ultrasound International Open 2: E19–E26. doi:10.1055/s-0035-1569403.
#' 
#' @export
trisomy21 <- function(age, weight, crlbe, crl, nt, pappa, betahCG, nasalbone, background) {
  ##Translate background number (i.e. the odds given as 1:x) to a prior probability
  prior <- 1/(background+1)
  
  ##Compute weight corrected logPAPPA and logbetahCG
  logPAPPA_corr   <- logPAPPA_corr(pappa, weight)
  logbetahCG_corr <- logbetahCG_corr(betahCG, weight) 
 
  ##Compute predicted values and reference range
  y <- c(nt=nt, logPAPPA=logPAPPA_corr, pappa=exp(logPAPPA_corr),logbetahCG=logbetahCG_corr, betahCG=exp(logbetahCG_corr))
  y_pred <- c(nt=nt_ls(crl), logPAPPA=logPAPPA(crlbe),logbetahCG=logbetahCG(crlbe))
  y_lower <- y_pred + c(nt_qBands[1],pappa_qBands[1],betahCG_qBands[1])
  y_upper <- y_pred + c(nt_qBands[2],pappa_qBands[2],betahCG_qBands[2])
  
  ##Compute degree of extremeness for the 3 markers
  doe_nt <- doe( nt, y_pred["nt"], y_lower["nt"], y_upper["nt"])
  doe_pappa <- doe( logPAPPA_corr, y_pred["logPAPPA"], y_lower["logPAPPA"], y_upper["logPAPPA"])
  doe_betahCG <- doe( logbetahCG_corr, y_pred["logbetahCG"], y_lower["logbetahCG"], y_upper["logbetahCG"])

  ##Store doEs as a list
  doe_list <- list(nt=doe_nt, pappa=doe_pappa,betahCG=doe_betahCG)

  ##Compute likelihood ratios - values from Merz et al. (2016)
  lr_nt <- lr( (doe_nt + 2.5)^(1/7), eu=c(1.13263, 0.04320),  aneu=c(1.23730, 0.09251))
  lr_pappa <- lr( doe_pappa,         eu=c(0.01629, 0.61358),  aneu=c(-0.85335, 0.73057))
  lr_betahCG <- lr( doe_betahCG,     eu=c(-0.01405, 0.61453), aneu=c(0.78042, 0.62854))
  
  #########
  ##Likelihood ratio of the event that no nasal bone is observed
  #########
  noNosalBoneObs <- (nasalbone == "No")
  ##Likelihood ratio for this according to Nicolaides (2004)
  lr_NoNasalBoneObs <- dbinom(noNosalBoneObs, size=1, prob=0.65) / dbinom(noNosalBoneObs, size=1, prob=0.02)
  lr_nasalbone <- ifelse(nasalbone == "Unknown", 1, lr_NoNasalBoneObs)
  
  ##Compute joint likelihood ratio
  lr_joint <- lr_nt * lr_pappa * lr_betahCG * lr_nasalbone
  
  ##Store LRs as a list
  lr_list <- list(nt=lr_nt, pappa=lr_pappa, betahCG=lr_betahCG, nasalbone=lr_nasalbone, joint=lr_joint)
  
  ##Compute posterior from likelihood ratio and prior
  posterior_probability <- ppost(lr_joint, pi=prior)
  
  ##Risk is given in odds form as 1:y, with y=
  tri21_onetox <- round((1-posterior_probability)/posterior_probability)

  ##Return everything as a long list.  
  return(list(PAPPA_corr=exp(logPAPPA_corr),
              betahCG_corr=exp(logbetahCG_corr),
              y=y, y_pred=y_pred, y_lower=y_lower, y_upper=y_upper,
              posterior=posterior_probability, 
              onetox=tri21_onetox, doe_list=doe_list,
              lr_list=lr_list))
}


#' Helper function to plot a likelihood ratio between the eu and aneu population
#' 
#' @param what Which biomarker to plot
#' @param what_title The title of the plot
#' @param what_xlab The xlab to use
#' @param eu Vector of length two containing the eu mean and standard deviation
#' @param aneu Vector of length two containing the aneu mean and standard deviation
#' @param res21 Names list containing results of the DoE, likelihood ratio, etc. for each of the three biomarkers (+ nasal bone information)
#' @param f Transformation function for the DoE
#' @return A `ggplot` object
#' @seealso [trisomy21risk::trisomy21]
#' @export
plotLR <- function(what, what_title, what_xlab, eu, aneu, res21, f=identity, f_inv=identity, diff_f_inv=function(x) rep(1,length(x)), f_obs=f, pop_quantiles=c(0.001,0.99)) {
  
  ##Range of doE values to look at
  range <- range(res21$doe_list[[what]],
                 f(qnorm(pop_quantiles, mean=eu[1], sd=eu[2])),
                 f(qnorm(pop_quantiles, mean=aneu[1], sd=aneu[2])))
 
  ##Compute densities on a grid
  densities <- data.frame(x=seq(range[1],range[2],length=1000)) %>% 
    mutate(Normal =   dnorm(f_inv(x), mean=eu[1], sd=eu[2])     * abs(diff_f_inv(x)),
           T21 = dnorm(f_inv(x), mean=aneu[1], sd=aneu[2]) * abs(diff_f_inv(x)))
  
  ##Convert to long format
  densities_long <- densities %>% tidyr::gather(key="Population", value="pdf",-x)
  
  ##90 tolerance interval in normal population
  q90 <- f(qnorm(c(0.1,0.9), mean=eu[1], sd=eu[2]))
  poly <- densities_long %>% filter(x >= q90[1] & x <= q90[2] & Population == "Normal")
  polygon <- data.frame( x=c(poly$x, rev(poly$x)),pdf=c(poly$pdf,rep(0,nrow(poly))),Population="Normal")
  
  ##Make the plot
  ggplot(densities_long, aes(x=x, y=pdf,color=Population)) + 
    geom_polygon(data=polygon, aes(x=x, y=pdf, color=Population),fill="aquamarine3") +
    geom_line(lwd=1.2) +
    geom_vline(xintercept=res21$y[[what]], pdf, lty=2, col="black") + 
    xlab(what_xlab) + ylab("Density") + 
    ggtitle(paste0(what_title," - LR = ",sprintf("%.3f",res21$lr_list[[what]]))) + 
    theme(legend.position="bottom") + 
    #theme(legend.position="bottom", text = element_text(size=16))
    scale_color_manual(values=c("Normal"="aquamarine3","T21"="darkgoldenrod3")) #darkorange3
  
}
