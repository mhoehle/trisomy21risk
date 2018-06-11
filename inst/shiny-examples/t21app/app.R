#######################################################
##
## This is a Shiny web application to calculate the risk
## of chromosone abnormalities, in particular trisomy 21,
## after first trimester scanning using a combined
## test.
##
## Author: Michael Höhle <http://www.math.su.se/~hoehle>
## Date:   June 2018
## License: GNU-GPL v3.0
#######################################################

library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(gridExtra)

##Load special trisomy 21 functions as part of a github package
library(trisomy21risk)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Risk of trisomy 21 after First-Trimester Screening"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         numericInput("age_mother","Age of Mother [years]:", value=33, min=20, max=45),
         numericInput("weight_mother","Weight of Mother [kg]:", value=65.1, min=20, max=120,step=0.1),
         numericInput("background","Background risk as 1:x", value=349, min=1, max=2000,step=1),
         numericInput("nt","Nuchal Translucency (NT) [mm]:", value=2.0, min=0, max=6, step=0.1),
         numericInput("crl","Crown Rump Length (CRL) at NT [mm]:", value=round(gestage2crl(86),digits=0), min=45, max=85,step=1),
         numericInput("crlbe","Crown Rump Length at blood sampling [mm]:", value=round(gestage2crl(58),digits=0), min=16, max=85),
         numericInput("PAPPA","PAPP-A [miU/ml]:", value=0.526, min=0, max=5,step=0.1),
         numericInput("freebetahCG","βhCG [miU/ml]:", value=124.7, min=0, max=200,step=1),
         selectInput("nasalbone","Nasal Bone visible on Ultrasound:", choices=c("Yes","No","Unknown"),selected="Unknown",multiple=FALSE)
      ),
      
      # Show a plot of the generated distribution
      mainPanel = mainPanel(
        tabsetPanel(
          tabPanel("Risk Assessment",  
                   htmlOutput("dataInfo"),
                   hr(),
                   htmlOutput("posteriorInfo")),
          tabPanel("Likelihood Ratios", 
                   p(),
                   plotlyOutput("plotLR_nt"),
                   p(), p(), p(),
                   plotlyOutput("plotLR_pappa"),
                   p(), p(), p(),
                   plotlyOutput("plotLR_betahCG"),
                   h3("Description"),
                   br("The plots above illustrate the computation of the likelihood ratio (LR) for the NT, PAPP-A and the βhCG measurements, respectively. For each of the three measurements a likelihood ratio is computed by taking the ratio between the density of the trisomy 21 fetus population (T21) vs. the density of the normal fetus population (Normal) evaluated at the observed measurement value. Likelihood ratios > 1 are indicative of trisomy 21. For better visualization an equi-tailed 90% tolerance interval for the normal population is shown as shaded green area, but what matters for the diagnostic procedure is how the T21 population differs from the Normal population. For this purpose the vertical dashed line shows the observed value at which the two densities are evaluated to get the LR. Note that in the case of PAPP-A and βhCG the weight-corrected value is used. Each plot can be hovered, panned and zoomed using ",a(href="https://help.plot.ly/zoom-pan-hover-controls/","plotly."))
                   
                   
          ),
          tabPanel("Information",
                   h3("Information"),
                   "This app is based on the procedure described in the paper",
                   a(href="https://www.thieme-connect.com/products/ejournals/abstract/10.1055/s-0035-1569403",em("Prenatal Risk Calculation (PRC) 3.0: An Extended DoE-Based First-Trimester Screening Algorithm Allowing For Early Blood Sampling")),
                   "by Merz et al. (2016) and which is used by the PRC 3.0 software of the ",a(href="http://www.fmf-deutschland.info","Fetal Medicine Foundation Germany.")," Numerical values will deviate slightly from the paper due to problems with the NT polynomial formula in the paper. Instead, a manual least squares fit of a similar polynomial as in the paper has been performed to the NT curve presented in the paper's Fig. 3.",
                   p(),
                   "The statistical details of the Bayesian diagnostic test procedure of the form ",p(),
                   p("posterior odds = likelihood ratio × prior odds",style="margin-left: 120px"),
                   "can be found in the paper and the blog entry ",a(href="","First Trimester Combined Testing for Trisomy 21"),".", 
                   "The R code for the calculations as well as the Shiny App are available under an open-source GPL v3.0 license as part of the ", code("trisomy21risk"), "package available from ",a(href="","github."),
                   p(),
                   "The background risk by maternal age and week of gestation is taken from Table 4 of the paper ",
                   a(href="http://onlinelibrary.wiley.com/doi/10.1046/j.1469-0705.1999.13030167.x/abstract", em("Maternal age- and gestation-specific risk for trisomy 21")),
                   "by Snijders et al. (1999), which is behind a paywall. The crown rump length is translated to gestational age using the equation ",a(href="http://journals.sagepub.com/doi/pdf/10.1179/174313409X448543", em("Gestational age = sqrt(CRL x 1.037) x 8.052 + 23.73.")),"If this gestational age is in week 10, 12 or 14 the values are taken directly from the table. For weeks 11 and 13 a linear interpolation between the risk of the neighbouring weeks is done. The background risks obtained by this procedure appear to differ slightly from the background risk calculation in the Merz et al. (2016) paper which states to use the same tables. You also have the opportunity to enter your individual brackground risk. However, whenever you change the maternal age the background risk is automatically re-calculated.", 
                   p(),
                   "The nasal bone information is taken from p.21 in ",a(href="http://www.fetalmedicine.com/fmf/FMF-English.pdf","Nicolaides (2004)")," where it says that the nasal bone is ",em("not visible by ultrasound in about 60-70% of fetuses with trisomy 21 and in about 2% of chromosomally normal fetuses."), "If the no information about the visibility of the nasal bone is available a LR of 1 is used.",
                   p(),
                   h3("Further Information"),
                   "Further information about first trimester screening can be obtained from the ",
                   a(href="https://fetalmedicine.org/", "Fetal Medicine Foundation"),
                   "in particular their ",
                   a(href="http://www.fetalmedicine.com/fmf/FMF-English.pdf","book"),
                   "on the 11-13 weeks scan (available in several languages).",
                   "The FMF also offers a ",
                   a(href="https://fetalmedicine.org/research/assess/trisomies","trisomies risk assessment calculator."),
                   p(),
                   h3("Disclaimer"),
                   "Aim of the present app is to visualize the information provided by the three biomarkers for educational and transparency purposes. The app provides comparable, but in no way identical, values to the PRC 3.0 software. Please note: Use of the app is no replacement for the actual supervised calculations performed by such software nor for the subsequent consultation with a gynaecologist. Absolutely no warranty is given for the results."
          )
        ))
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
  observe({
    # This will change the value of input$inText, based on x
    if (!is.na(input$crl) & !is.na(input$age_mother)) {
      updateNumericInput(session, "background", value = round(background_risk(age_mother=input$age_mother, gestage=crl2gestage(input$crl)),digits=0))
    }
  })
  
   output$dataInfo <- renderUI({
     ##Calculate gestation age
     derived_gestage <- crl2gestage(input$crl) #round(crl2gestage(input$crl),digits=0)
   
     
     HTML(paste0("<h4>Input Data</h4>",
                 "Age of Mother: ", input$age_mother, " years <br>",
                 "Weight of Mother: ", input$weight_mother, " kg <br>",
                 "Current gestational age (calculated from CRL at NT): ",round(derived_gestage,digits=0), " days (",floor(derived_gestage/7),"<sup>+",round(derived_gestage - floor(derived_gestage/7)*7,digits=0),"</sup> weeks)<br>",
                 "Background risk: 1:", input$background, "<br>",
                 "<p><p>",
                 "Nuchal Translucency: ",sprintf("%.1f",input$nt), " mm<br>",
                 "Crown Rump Length at NT: ",input$crl, " mm<br>",
                 "<p>",
                 "PAPP-A: ",input$PAPPA, " miU/ml<br>",
                 "Weight corrected PAPPA-A: ",sprintf("%.3f",exp(logPAPPA_corr(input$PAPPA, input$weight_mother))),  " miU/ml<br>",
                 "βhCG: ",input$freebetahCG, " miU/ml<br>",
                 "Weight corrected βhCG: ",sprintf("%.1f",exp(logPAPPA_corr(input$freebetahCG, input$weight_mother))),  " miU/ml<br>",
                 "Crown Rump Length at blood sampling: ",input$crlbe, " mm<br>",
                 "<p>",
                 "Nasal bone visible on ultrasound: ", input$nasalbone))
   })
   
   output$posteriorInfo <- renderUI({
     res21 <- trisomy21(input$age_mother, input$weight_mother, crlbe=input$crlbe, crl=input$crl, nt=input$nt, pappa=input$PAPPA, betahCG=input$freebetahCG, nasalbone=input$nasalbone, background=input$background) 
     
     ##print(res21)
        
     HTML(paste0("<h4>Risk Assessment</h4>",
                "Prior Odds of T21:", "1:",input$background,"<br>",
                "LR contribution NT: ",sprintf("%.3f",res21$lr_list$nt), " (DoE ",sprintf("%.2f",res21$doe_list$nt),")<br>",
                "LR contribution weight corrected PAPP-A: ",sprintf("%.3f",res21$lr_list$pappa), " (DoE ",sprintf("%.2f",res21$doe_list$papp),")<br>",
                "LR contribution weight corrected betahCG: ",sprintf("%.3f",res21$lr_list$betahCG), " (DoE ",sprintf("%.2f",res21$doe_list$betahCG),")<br>",
                "LR contribution nasal bone: ",sprintf("%.3f",res21$lr_list$nasalbone), "<p><p>",
                "LR joint = (LR NT) × (LR PAPP-A) × (LR betahCG) × (LR Nasal Bone): ",sprintf("%.3f",res21$lr_list$joint), "<p>",
                #"1/(LR joint): ",sprintf("%.3f",1/res21$lr_list$joint), "<p><p>",
                "Posterior Odds of T21 = (LR joint) × (Prior Odds of T21): <b>1:", res21$onetox,"</b>"))
   })
   
   #output$plotLR_nt <- renderPlot({
   output$plotLR_nt <- renderPlotly({
     ##Compute the risk
     res21 <- trisomy21(input$age_mother, input$weight_mother, crlbe=input$crlbe, crl=input$crl, nt=input$nt, pappa=input$PAPPA, betahCG=input$freebetahCG, nasalbone=input$nasalbone, background=input$background) 
     
     ##NT values from Table 3 of Merz et al. (2016)
     eu <- c(mean=1.13263, sd=0.0432)
     aneu <- c(mean=1.23730, sd=0.09251)
     
     ##Extract nt values
     y_pred <- res21$y_pred["nt"]
     y_ref <- c(res21$y_lower["nt"], res21$y_upper["nt"])
     
     ##Transform functions of the doe of NT, see Table 3 of Merz et al. (2016)
     # f <- function(x) x^7 - 2.5
     # f_inv <- function(y) (y+2.5)^(1/7)
     # diff_f_inv <- function(y) 1/(7 * (y+2.5)^(6/7))
     # 
     
     ##The transformation functions. Identity if not transformed.
     f <- function(x) {
       doe_inv(x^7-2.5, y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2])
     }
     
     f_inv <- function(y) {
       (doe(y, y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2]) + 2.5)^(1/7)
     }
     
     diff_f_inv <- function(y) {
       ref_lower <- y_ref[1] ; ref_upper <- y_ref[2]
       r <- ifelse(y>y_pred,
                   1 / (ref_upper - y_pred),
                   1 / (y_pred - ref_lower)
       )
       1/(7 * (doe(y, y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2])+2.5)^(6/7)) * r
     }
     
     p <- plotLR("nt", "NT", "NT [mm]", eu=eu, aneu=aneu, res21=res21, f=f, f_inv=f_inv, diff_f_inv=diff_f_inv, f_obs=identity)
     plotly::ggplotly(p)
   })
   
   output$plotLR_pappa <- renderPlotly({
     ##Compute the risk
     res21 <- trisomy21(input$age_mother, input$weight_mother, crlbe=input$crlbe, crl=input$crl, nt=input$nt, pappa=input$PAPPA, betahCG=input$freebetahCG, nasalbone=input$nasalbone, background=input$background) 
    
     ##PAPP-A alues from Table 3 of Merz et al. (2016)
     eu <- c(mean=0.01629, sd=0.61358)
     aneu <- c(mean=-0.85335, sd=0.73057)

     y_pred <- res21$y_pred["logPAPPA"]
     y_ref <- c(res21$y_lower["logPAPPA"], res21$y_upper["logPAPPA"])
     
     ##The transformation functions. Identity if not transformed.
     f <- function(x) {
       exp(doe_inv(x, y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2]))
     }
     
     f_inv <- function(x) {
       doe(log(x), y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2])
     }
     
     diff_f_inv <- function(x) {
       y <- log(x)
       ref_lower <- y_ref[1] ; ref_upper <- y_ref[2]
       ifelse(y>y_pred,
              1 / (ref_upper - y_pred) * abs(1/x),
              1 / (y_pred - ref_lower) * abs(1/x)
       )
     }
     
     p <- plotLR("pappa", "PAPP-A", "PAPP-A weight corrected [miU/ml]", eu=eu, aneu=aneu, res21=res21, f=f, f_inv=f_inv, diff_f_inv=diff_f_inv)
     plotly::ggplotly(p)
   })
   
   output$plotLR_betahCG <- renderPlotly({
     ##Compute the risk
     res21 <- trisomy21(input$age_mother, input$weight_mother, crlbe=input$crlbe, crl=input$crl, nt=input$nt, pappa=input$PAPPA, betahCG=input$freebetahCG, nasalbone=input$nasalbone, background=input$background) 
     
     ##free beta-hCG values from Table 3 of Merz et al. (2016)
     eu <- c(mean=-0.01405, sd=0.61453)
     aneu <- c(mean=0.78042, sd=0.62854)

     y_pred <- res21$y_pred["logbetahCG"]
     y_ref <- c(res21$y_lower["logbetahCG"], res21$y_upper["logbetahCG"])
  
     ##The transformation functions. Identity if not transformed.
     f <- function(x) {
       exp(doe_inv(x, y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2]))
     }
     
     f_inv <- function(x) {
       doe(log(x), y_pred=y_pred, ref_lower=y_ref[1], ref_upper=y_ref[2])
     }
     
     diff_f_inv <- function(x) {
       y <- log(x)
       ref_lower <- y_ref[1] ; ref_upper <- y_ref[2]
       ifelse(y>y_pred,
              1 / (ref_upper - y_pred) * abs(1/x),
              1 / (y_pred - ref_lower) * abs(1/x)
       )
     }
     
     p <- plotLR("betahCG", "βhCG", "βhCG weight corrected [miU/ml]", eu=eu, aneu=aneu, res21=res21,f=f, f_inv=f_inv, diff_f_inv=diff_f_inv)
     plotly::ggplotly(p)
   })
     
}

# Run the application 
shinyApp(ui = ui, server = server)

