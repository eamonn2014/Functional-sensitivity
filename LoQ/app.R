 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls()) 
set.seed(333) # reproducible
library(directlabels)
library(shiny) 
library(shinyjs)  #refresh'
library(shinyWidgets)
library(shinythemes)  # more funky looking apps
library(shinyalert)
library(Hmisc)
library(reshape)
library(rms)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(shinycssloaders)
#library(googleVis)
library(xtable)

options(max.print=1000000)    

fig.width <- 1200
fig.height <- 500
fig.width1 <- 1380
fig.width8 <- 1380
fig.height1 <- 700
fig.width7 <- 700
fig.height7 <- 500
fig.width6 <- 680

## convenience functions
p0 <- function(x) {formatC(x, format="f", digits=0)}
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
p3 <- function(x) {formatC(x, format="f", digits=3)}
p4 <- function(x) {formatC(x, format="f", digits=4)}
p5 <- function(x) {formatC(x, format="f", digits=5)}

logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful

options(width=200)
options(scipen=999)
w=4  # line type
ww=3 # line thickness
wz=1 

# not used,  MSE for linear regression
calc.mse <- function(obs, pred, rsq = FALSE){
    if(is.vector(obs)) obs <- as.matrix(obs)
    if(is.vector(pred)) pred <- as.matrix(pred)
    
    n <- nrow(obs)
    rss <- colSums((obs - pred)^2, na.rm = TRUE)
    if(rsq == FALSE) rss/n else {
        tss <- diag(var(obs, na.rm = TRUE)) * (n - 1)
        1 - rss/tss
    }
}

RR=.37  # used to limit correlations between variables
sd1= 3  # for X covariates sd

# links to Rdata objects uploaded to Git, these are pre run simulations see the save function
# see line 1224 for the save function
pp<- "https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata" # 5000 default log1.5 -1 -.67 -.43
pp2<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/B%205000%20default%20settings%20theta%20log0.5%20-1.68%20-1.39%20%200.71.Rdata"
pp3<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/C%205000%20default%20settings%20theta%20log2%20-3.46%20-1.05%20%201.15%20p1=.75.Rdata"
pp4<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/D_10Ksims_5covariates_p1_0.12_theta_log1.3_covariates_-1.02%20_0.42_0.43_0.61%20_1.01_3_prog.Rdata"


loq <- function (x, y, model, spec, print.plot=1) {
    
    # Define models
    
    if (model %in% 1 ) {mod="Linear Y=a+bX"}  
    if (model %in% 2 ) {mod="Exponential Y=exp(a+bX)"} 
    if (model %in% 3 ) {mod="Reciprocal-Y Y=1/(a+bX)"} 
    if (model %in% 4 ) {mod="Reciprocal-X Y=a+b/X"} 
    if (model %in% 5 ) {mod="Double Reciprocal Y=1/(a+b/X)"} 
    if (model %in% 6 ) {mod="Logarithmic-X Y=a+b(log(X))"} 
    if (model %in% 7 ) {mod="Multiplicative Y=aX^b"} 
    if (model %in% 8 ) {mod="Square Root-X Y=a+b(sqrt(X))"} 
    if (model %in% 9 ) {mod="Square Root-Y Y=(a+bX)^2"} 
    if (model %in% 10) {mod="S-curve Y=exp(a+b/X)"} 
    if (model %in% 11) {mod="Square X and Y Y^2=a+X^2/b"} 
    
    # transformation of data for 11 models
    
    ty1 <- y;       tx1 <- x
    ty2 <- log(y);  tx2 <- x
    ty3 <- 1/y;     tx3 <- x
    ty4 <- y;       tx4 <- 1/x
    ty5 <- 1/y;     tx5 <- 1/x
    ty6 <- y;       tx6 <- log(x)
    ty7 <- log(y);  tx7 <- log(x)
    ty8 <- y;       tx8 <- sqrt(x)
    ty9 <- sqrt(y); tx9 <- x
    ty10 <- log(y); tx10 <- 1/x
    ty11 <- y^2;    tx11 <- x^2
    
    # save the original data
    
    x1 <- x
    y1 <- y
    
    # transform using the selected model (1 - 11)
    
    x <- eval(parse(text=(paste("tx", model, sep="")) )) 
    y <- eval(parse(text=(paste("ty", model, sep="")) )) 
    
    # summary statistics of the independent variable
    
    n <- length(x)
    txbar <- mean(x)
    txstd <- sd(x)
    
    # run regression on the transformed data grab slope intercept
    
    f <- lm(y~x) 
    intercept <- coef(f)[1][[1]]
    slope <- coef(f)[2][[1]]
    
    # R2 on the transformed x and y,  an idea but not used
    
    v1 <- unlist(cor.test(x,y)$estimate^2)[1][[1]]   
    
    # obtain the predictions 
    
    p <- predict.lm(f, interval="confidence")
    
    # transform the predicted values & 95%CI back to original scale 
    
    if (model %in% c(2,7,10)) {p <- exp(p)} 
    if (model %in% c(3,5)  )  {p <- 1/p} 
    if (model %in% c(9)    )  {p <- p^2} 
    if (model %in% c(11)   )  {p <- p^.5} 
    
    # calculate residual squared 
    # residuals original y and transformed back predicted values 
    
    r <- (y1-p[,1])^2 
    
    # R2 on the original x and transformed back predicted y, an idea but not used
    
    v2 <- unlist(cor.test(x1,p[,1])$estimate^2)[1][[1]]
    
    # residual sum of squares, this will be used to judge best model
    
    ssr <- sum(r, na.rm=T) 
    
    # transform the response that we will read back
    
    tyspec <- spec 
    if (model %in% c(2,7,10)) {tyspec <- log(tyspec)} 
    if (model %in% c(3,5)  )  {tyspec <- 1/tyspec} 
    if (model %in% c(9)    )  {tyspec <- sqrt(tyspec)} 
    if (model %in% c(11)   )  {tyspec <- tyspec ^2}  
    
    # grab the residual standard deviation
    
    rsd2<-as.data.frame(anova(f))[2,3]^.5  
    
    # read back on transformed scale
    
    mse <- rsd2^2
    t <- qt(0.975, n-2)
    a <- t^2*mse/((n-1)*txstd^2)-slope^2
    b <- 2*slope*(tyspec-intercept-slope*txbar)
    c <- t^2*mse/n-(tyspec-intercept-slope*txbar) ^2
    txpre <- (tyspec-intercept)/slope
    txup <-  (-b+sqrt(b^2-4*a*c))/(2*a) + txbar
    txlow <- (-b-sqrt(b^2-4*a*c))/(2*a) + txbar
    
    # transform the read back estimates to the original scale
    
    if (model %in% c(4,5,10)) {txpre <- 1/txpre; txup <- 1/txup; txlow <- 1/txlow} 
    if (model %in% c(6,7)  )  {txpre <- exp(txpre); txup <- exp(txup); txlow <- exp(txlow)}  
    if (model %in% c(8)    )  {txpre <- txpre^2;  txup <- txup^2;  txlow <- txlow^2} 
    if (model %in% c(11)   )  {txpre <- txpre^.5; txup <- txup^.5; txlow <- txlow^.5}   
    
    # ensure order of limits is correct
    
    limits <- sort(c(txlow,txup))
    txlow <- limits[1]
    txup <-  limits[2]
    
    # help with plotting
    
    ymin <- min(y)
    ymax <- max(y)
    ystep <- (ymax-ymin)/8
    ymin1 <-  ymin-ystep
    ymax1 <-  ymax+ystep
    
    xmin <- min(x)
    xmax <- max(x)
    xstep <- (xmax-xmin)/8
    xmin1 <-  xmin-xstep
    xmax1 <-  xmax+xstep
    
    # put all pertinent data together, original data and predicted with 95%CI
    
    foo <- data.frame(cbind(x=x1,obsy=y1, pred= p[,1], p2a=p[,2], p3=p[,3]))
    foo <- foo[order(foo$obsy),]
    
    # plot and present the estimated read back
    
    p1 <- ggplot(foo, aes(x=x,y=pred)) + 
        geom_line( ) +
        geom_ribbon(data=foo , aes(ymin= p2a,ymax= p3),alpha=0.2,fill="blue") +
        geom_point(data=foo, aes(x=x ,y=obsy))  
    scale_x_continuous(limits = c(xmin1, xmax1))
    scale_y_continuous(limits = c(xmin1, xmax1))
    
    p <- p1  + geom_hline(yintercept=spec, colour="#990000", linetype="dashed")
    
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13,color="darkred"))
    p <- p + scale_color_manual(values=c("Red","blue"))
    p <- p + theme_bw()
    p <- p + xlab('Independent variable') 
    p <- p + ylab('Dependent variable')
    
    p <- p + ggtitle(paste("Model for the curve #",model," ",mod,", ","\nExploration of model fitting at response of", spec ,", ",p2f(txpre),"and 95% CI: (",p2f(txlow),",",p2f(txup),")","\nResidual sum of squares", p2f(ssr),", Residual standard deviation",p2f(rsd2),
                           # ", R2 transformed x, y ",p2(v1),", R2 x, trans. back pred. y",p2(v2),
                           sep=" ")) + theme(plot.title = element_text(lineheight=1, face="bold", color="black", size=11))
    
    p <- p + labs(x = "Independent variable", y = "Response")
    p <- p + theme(axis.title.y = element_text(size = rel(1.1), angle = 90))
    p <- p + theme(axis.title.x = element_text(size = rel(1.1), angle = 00))
    
    if (print.plot==1) {print(p)}
    
    
    return(ssr)
    
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


css <- "
#large .selectize-input { line-height: 40px; }
#large .selectize-dropdown { line-height: 30px; }"




ui <-  fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2 , paper another option to try
                 # paper
                 useShinyalert(),  # Set up shinyalert
                 setBackgroundColor(
                     color = c( "#2171B5", "#F7FBFF"), 
                     gradient = "linear",
                     direction = "bottom"
                 ),
                 
                 
                 
                 h2("xxxxxxxxxxxxxxxxx"), 
                 
                 h4("xxxxxxxxxxxxxx
                "), 
                 
                 h3("  "), 
                 sidebarLayout(
                     
                     sidebarPanel( width=3 ,
                                   
                                   tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                   
                                   actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/RCT-covariate-adjust-binary-response/master/cov-adj-binary-response/app.R', '_blank')"), 
                                   actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/RCT-covariate-adjust-binary-response/master/cov-adj-binary-response/R%20code.R', '_blank')"),  
                                   actionButton("resample", "Simulate a new sample"),
                                   br(),  
                                   tags$style(".well {background-color:#b6aebd ;}"), 
                                   
                                   # h4("User inputs"),
                                   div(
                                       # font colours for button font
                                       tags$head(
                                           tags$style(HTML('#upload{color:black}'))    
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload2{color:black}'))    
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload3{color:black}'))    
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload4{color:black}'))    
                                       ),
                                       
                                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                                       
                                       # colours for button background
                                       
                                       tags$head(
                                           tags$style(HTML('#upload{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload2{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload3{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#upload4{background-color:chartreuse}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#ab1{background-color:orange}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#resample{background-color:orange}'))
                                       ),
                                       
                                       tags$head(
                                           tags$style(HTML('#resample2{background-color:orange}'))
                                       ),
                                       tags$head(
                                           tags$style(HTML('#sim{background-color:orange}'))
                                       ),
                                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                                       splitLayout(
                                           
                                           textInput('N', 
                                                     div(h5(tags$span(style="color:blue", "N"))), "100"),
                                           textInput('a', 
                                                     div(h5(tags$span(style="color:blue", "intercept"))), "10"),
                                           
                                           textInput('b', 
                                                     div(h5(tags$span(style="color:blue", "slope"))), "1"),
                                           
                                           textInput('sigma', 
                                                     div(h5(tags$span(style="color:blue", "sigma"))), "2")
                                         
                                       ),
                                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                                       
                                 
                                       
                                       
                                       radioButtons(
                                           inputId = "truth",
                                           label =  div(h5(tags$span(style="color:blue","Data generating mechanism :"))),
                                           choiceNames = list(
                                            #   HTML("<font color='blue'>All scenarios</font>"), 
                                               tags$span(style = "color:blue", "Linear Y=a+bX"), 
                                               tags$span(style = "color:blue", "Exponential Y=exp(a+bX)"), 
                                               tags$span(style = "color:blue", "Reciprocal-Y Y=1/(a+bX)"),
                                               tags$span(style = "color:blue", "Reciprocal-X Y=a+b/X"),
                                               tags$span(style = "color:blue", "Double Reciprocal Y=1/(a+b/X)"),
                                               tags$span(style = "color:blue", "Logarithmic-X Y=a+b(log(X))"), 
                                               tags$span(style = "color:blue", "Multiplicative Y=aX^b"), 
                                               tags$span(style = "color:blue", "Square Root-X Y=a+b(sqrt(X))"),
                                               tags$span(style = "color:blue", "Square Root-Y Y=(a+bX)^2"),
                                               tags$span(style = "color:blue", "S-curve Y=exp(a+b/X)"),
                                               tags$span(style = "color:blue", "Square X and Y Y^2=a+X^2/b")
                                               
                                           ),
                                           choiceValues = c( "model1", "model2", "model3",  "model4", "model5", "model6",
                                                             "model7", "model8", "model9",  "model10", "model11")
                                       ),
                                       
                                       
           
                                       
                                       ###https://stackoverflow.com/questions/49616376/r-shiny-radiobuttons-how-to-change-the-colors-of-some-of-the-choices
                                       
                                       radioButtons(
                                           inputId = "dist",
                                           label =  div(h5(tags$span(style="color:blue","Analysis model :"))),
                                           choiceNames = list(
                                               HTML("<font color='blue'>Best scenario</font>"), 
                                               tags$span(style = "color:blue", "Linear Y=a+bX"), 
                                               tags$span(style = "color:blue", "Exponential Y=exp(a+bX)"), 
                                               tags$span(style = "color:blue", "Reciprocal-Y Y=1/(a+bX)"),
                                               tags$span(style = "color:blue", "Reciprocal-X Y=a+b/X"),
                                               tags$span(style = "color:blue", "Double Reciprocal Y=1/(a+b/X)"),
                                               tags$span(style = "color:blue", "Logarithmic-X Y=a+b(log(X))"), 
                                               tags$span(style = "color:blue", "Multiplicative Y=aX^b"), 
                                               tags$span(style = "color:blue", "Square Root-X Y=a+b(sqrt(X))"),
                                               tags$span(style = "color:blue", "Square Root-Y Y=(a+bX)^2"),
                                               tags$span(style = "color:blue", "S-curve Y=exp(a+b/X)"),
                                               tags$span(style = "color:blue", "Square X and Y Y^2=a+X^2/b")
                                               
                                           ),
                                           choiceValues = c( "best","model1", "model2", "model3",  "model4", "model5", "model6",
                                                             "model7", "model8", "model9",  "model10", "model11")
                                       )
                                       
                                   ),
                                   
                     ),
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                     mainPanel(width=9, #eight=4,
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                               navbarPage(       
                                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                                   tags$style(HTML("
                            .navbar-default .navbar-brand {color: orange;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #b6aebd;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
                   ")),
                                   
                                   tabPanel( "xxxxxxxxxxxxxxxxxxxx",
                                             
                                             
                                             
                                             div(class="span7", verbatimTextOutput("dato")),
                                          #   h4(htmlOutput("textWithNumber1a") ),
                                             fluidRow(
                                                 column(width = 6, offset = 0, style='padding:1px;',
                                                     #   shinycssloaders::withSpinner(
                                                     #       div(plotOutput("reg.plotx",  width=fig.width8, height=fig.height7)),
                                                   #     ),
                                               #         div(plotOutput("reg.ploty",  width=fig.width8, height=fig.height7)),
                                                ) ,
                                                 
                                                 
                                                 fluidRow(
                                                     column(width = 6, offset = 0, style='padding:1px;',
                                                            
                                                     ))),#
                                             
                                             h4(paste("xxxxxxxxxxxxxxx

                                                 ")),
                                             
                                             
                                             
                                             
                                             
                                             h4(paste("Table 2 xxxxxxxxxxxxx")),
                                             
                                           #  div( verbatimTextOutput("zz") )  ,
                                             # h4(htmlOutput("textWithNumber99",) ),
                                             # div( verbatimTextOutput("mse.target") )  ,
                                             h4(paste("xxxxxxxxxxxxxxxxx")),
                                        #     div( verbatimTextOutput("betas") )  ,
                                             width = 30 )     ,
                                
                                   
                                   # here have action buttons clicking will load a pre run simulation that can be examined
                                   tabPanel( "2 xxxxxxxxxxxxxxx",
                                             
                                             
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        actionBttn(
                                                            inputId = "upload",
                                                            label = "",
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE,
                                                            no_outline=TRUE
                                                        ),
                                                        
                                                 ),
                                                 
                                                 h4("Hit to load, default settings, the covariate coefficients used were -1 -0.67, -0.43"),
                                             ),
                                             
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        actionBttn(
                                                            inputId = "upload2",
                                                            label = "",
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE,
                                                            no_outline=FALSE
                                                        ),
                                                        
                                                 ),
                                                 
                                                 h4("Hit to load, default settings except that treatment effect is log(0.5). The covariate coefficients used were -1.68 -1.39 0.71."),
                                             ),
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        actionBttn(
                                                            inputId = "upload3",
                                                            label = "",
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE
                                                        ), 
                                                        
                                                 ),
                                                 
                                                 h4("Hit to load, default settings except that treatment effect is log(2), intercept probability 0.7. The covariate coefficients used were -3.46 -1.05 1.15"),
                                             ),
                                             
                                             
                                             
                                             fluidRow(
                                                 
                                                 column(1,
                                                        
                                                        
                                                        actionBttn(
                                                            inputId = "upload4",
                                                            label = "",  
                                                            color = "royal",
                                                            style = "float",
                                                            icon = icon("sliders"),
                                                            block = TRUE
                                                        ),
                                                        
                                                        
                                                        
                                                        
                                                 ),
                                                 h4("Hit to load, 10K simulations, 5 covariates (3 prognostic), treatment effect is log(1.3), intercept prob. 0.12. The covariate coefficients used were -1.02  0.42  0.43  0.61  1.01"),
                                                 
                                             ),
                                             
                                             shinycssloaders::withSpinner(plotOutput("plot1",  width=fig.width8, height=fig.height7),5),
                                             
                                             shinycssloaders::withSpinner(plotOutput("plot2",  width=fig.width8, height=fig.height7),5),
                                             
                                             shinycssloaders::withSpinner(verbatimTextOutput("content1"),type = 5),
                                             
                                             
                                   ),            
                                   
                                   tabPanel("3 xxxxxxxxxxxxxxxxxxxx", value=3, 
                                            
                                            h4("xxxxxxxxxxxxxxxxx
                                               ") ,
                                            tags$hr(),
                                            
                                            h4("1 xxxxxxxxxxxxxxxxxxxxx'") ,
                                            
                                            
                                            #https://github.com/daattali/advanced-shiny/tree/master/select-input-large
                                            tags$style(type='text/css', css),
                                            
                                            div(id = "large",
                                                selectInput("Trueeffect", "", width='35%',
                                                            c("Expected odds ratio" = "ORx",
                                                              "Anticipated proportion of responders in treated" = "p2x" ))
                                            ),
                                            
                                            h4("2 xxxxxxxxxxxxxxxxxxxxxxx'") ,
                                            textInput('pp2', 
                                                      div(h5(tags$span(style="color:blue", ""))), ".35"),
                                            
                                            h4("3 xxxxxxxxxxxxxxxxxxxxxxxx") ,
                                            splitLayout(
                                                
                                                textInput('NN', 
                                                          div(h5(tags$span(style="color:blue", "Total sample size"))), "300"),
                                                
                                                textInput('pp1', 
                                                          div(h5(tags$span(style="color:blue", "Expected proportion of responders in baseline/placebo"))), ".25"),
                                                
                                                textInput('allocation', 
                                                          div(h5(tags$span(style="color:blue", "Randomisation allocation"))), "0.5")
                                                
                                            ),
                                            
                                            tags$hr(),
                                            actionButton("sim","Hit to assess power of the design based on above inputs"),
                                            h4("Power via 499 simulations"),    
                                            withSpinner(verbatimTextOutput("pow1")),
                                            h4("Power via Frank Harrell Hmisc function"),    
                                            withSpinner(verbatimTextOutput("pow2")),
                                            tags$hr(),
                                            actionButton("resample2", "Hit to Simulate another sample based on above inputs"),  
                                            h4("xxxxxxxxxxxxxxxxxxx."),  
                                            
                                            
                                            span(
                                                style = "color: #000000; font-face: bold;",
                                                tableOutput("obs")),
                                            h4(htmlOutput("textWithNumber",) ),
                             
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            withSpinner(verbatimTextOutput("pow3")),
                                            # verbatimTextOutput("pow") %>% withSpinner(color="#0dc5c1"))
                                            
                                            h4(""),                           
                                            
                                   ),
                                   
                                   
                                   tabPanel("4 Notes & references", value=3, 
                                            
                                            h4("xxxxxxxxxxxxxxxxx") ,
                                            
                                            h4("xxxxxxxxxxxxxxx") ,
                                            
                                            
                                            
                                            h4("xxxxxxxxxxxxxxx "),
                                            h4("
                                          xxxxxxxxxxxxxxxxxxx
                                          "),
                                            h4("
                                          xxxxxxxxxxxxxxxxxxx "),
                                            
                                            h4("
                                               xxxxxxxxxxxxxxxxxxx
                                               
                                               "),
                                            
                                            h4("
                                              xxxxxxxxxxxxxxxxxxx'
                                               
                                               "),
                                            
                                            column(width = 12, offset = 0, style='padding:1px;',
                                                   
                                                   tags$hr(),
                                                   div(h4("References:")),  
                                                   tags$a(href = "https://www.bmj.com/content/bmj/340/bmj.c869.full.pdf", tags$span(style="color:blue", "[1] CONSORT 2010 Explanation and Elaboration: updated guidelines for reporting parallel group randomised trials"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://www.linkedin.com/pulse/stop-obsessing-balance-stephen-senn/", tags$span(style="color:blue", "[2] Stephen Senn, Stop obsessing about balance"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://discourse.datamethods.org/t/should-we-ignore-covariate-imbalance-and-stop-presenting-a-stratified-table-one-for-randomized-trials/547/32", tags$span(style="color:blue", "[3] Stephen Senn, point 4, Should we ignore covariate imbalance and stop presenting a stratified table one for randomized trials"),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://twitter.com/f2harrell/status/1298640944405807105",  tags$span(style="color:blue", "[4]  Frank Harrell, twitter, 'unadjusted analysis makes the most severe assumptions of all (that risk factors do not exist)'."),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://statistics.fas.harvard.edu/files/statistics/files/21_stephen_senn.pdf", tags$span(style="color:blue", "[5] Randomisation isn’t perfect but doing better is harder than you think "),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8570", tags$span(style="color:blue", "[6] Graphical calibration curves and the integrated calibration index (ICI) for survival models, Statistics in Medicine. 2020;1–29 "),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://www.sciencedirect.com/science/article/abs/pii/S0002870300900012?via%3Dihub", tags$span(style="color:blue", "[7] Steyerberg, E. W., Bossuyt, P. M. M., & Lee, K. L. (2000). Clinical trials in acute myocardial infarction: Should we adjust for baseline characteristics? American Heart Journal, 139(5), 745–751. doi:10.1016/s0002-8703(00)90001-2"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "http://clinicalpredictionmodels.org/", tags$span(style="color:blue", "[8] Steyerberg, E. W., Clinical Prediction Models, 2019 p459"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://twitter.com/f2harrell/status/1299755896319475712", tags$span(style="color:blue", "[9] Frank Harrell, twitter, Adjusted analysis"),),   
                                                   div(p(" ")),
                                                   tags$a(href = "https://discourse.datamethods.org/t/guidelines-for-covariate-adjustment-in-rcts/2814/2", tags$span(style="color:blue", "[10 Frank Harrell, Guidelines for covariate adjustment in rcts"),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://www.fharrell.com/post/covadj/", tags$span(style="color:blue", "[11] E.Steyerberg explains some of the advantages of conditioning on covariates"),),  
                                                   div(p(" ")),
                                                   tags$a(href = "https://hbiostat.org/doc/bbr.pdf", tags$span(style="color:blue", "[12] Biostatistics for Biomedical Research Frank E Harrell Jr. James C Slaughter Updated August 5, 2020"),),  
                                                   div(p(" ")),
                                                   
                                                   tags$hr()
                                            ) 
                                            
                                   )
                                   
                                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   END NEW   
                               )
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     )
                 ) 
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

server <- shinyServer(function(input, output   ) {
    
    #############################
    
    shinyalert("Welcome! \nFitting models",
               "Have fun!", 
               type = "info")
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is where a new sample is instigated and inputs converted to numeric
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    random.sample <- reactive({
        
        foo <- input$resample
        
        # a <- as.numeric(unlist(strsplit(input$intercept,",")))
        # 
        # b <- as.numeric(unlist(strsplit(input$slope,",")))
        # 
        # sigma <- as.numeric(unlist(strsplit(input$sigma,",")))
        # 
        # N <- as.numeric(unlist(strsplit(input$N,",")))
        # 
        a <- as.numeric(input$a)
        
        b <-as.numeric(input$b)
        
        sigma <- as.numeric(input$sigma)
 
        N <- as.numeric(input$N)
        
        
        return(list(  
            a=a,
            b=b,
            sigma=sigma,
            N=N
        ))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # tab 1 simulate data (covariates and response)  
    # create  response with prognostic covariate
    # create covariates that are not prognostic
    # create a mix of above 2
    # alos look at the difference of the covariates across arms
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 dat <- reactive({
    
    sample <- random.sample()
    
    a=sample$a
    b=sample$b
    sigma=sample$sigma
    N=sample$N
     
    
    x <-  array(runif(N, 1, 100))  # no negative values
    noise <-  rnorm(N,0, sigma)
 
    if (input$truth %in% "model1") {
        y <-  a+ x*b +    noise
    } else if (input$truth %in% "model2") {
        y <-  exp(a+ x*b + noise)
    } else if (input$truth %in% "model3") {
        y <-  1/(a+ x*b +  noise)
    } else if (input$truth %in% "model4") {
        y <-  a + b*(1/x) + noise
    } else if (input$truth %in% "model5") {
        y <-  1/(a + b/x +  noise)
    } else if (input$truth %in% "model6") {
        y <-  a + log(x)*b +    noise   
    } else if (input$truth %in% "model7") {
        y <-  a * x^b +    noise
    } else if (input$truth %in% "model8") {
        y <-  a + sqrt(x)*b +    noise
    } else if (input$truth %in% "model9") {
        y <-  (a + x*b + noise)^2 
    } else if (input$truth %in% "model10") {
        y <-  exp(a+ b/x + noise)
    } else if (input$truth %in% "model11") {
        y <-  ( a + (x^2)/b + noise)^2
    }

    return(list(  y=y))
    
 })
    
    
    
    
    
    
    
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$dato <- renderPrint({
        
      return(dat()$y)
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # simulation code 2 do the same as the first simulation code, but this time correlated covariates are created
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    simul2 <- reactive({
        
       # sample <- random.sample()
     
       
        
        return(list(  
          
        )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # do the same as the first simulation code, but this time imbalanced covariates are created
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # SIMUALTION PLOT
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # collect simulatio standard error estimates from simulation and plot!
    
    output$reg.plotyy <- output$reg.ploty <- renderPlot({         #standard errors
        
        # Get the  data
        
        
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # not used for binary logistic
    output$textWithNumber99 <- renderText({ 
        
        HTML(
            "xxxxxxxxxxxxxxxxx"
        )
        
    })  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # table for simulation summary, quite easy to make a mistake here! too much hard coding really, need to code better
   
    
    output$pow1 <- renderPrint({
        
       # return(twobytwo()$pow)
        
    }) 
    
    
    output$pow2 <- renderPrint({
        
      #  return(twobytwo()$FH)
        
    }) 
    
    output$pow3 <- renderPrint({
        
     #   return(mdata()$fit1)
        
    })
    
    
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loading in user data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


# Run the application 
shinyApp(ui = ui, server = server)