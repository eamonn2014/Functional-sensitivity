#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# White-margined Burrower Bug nymph - Sehirus cinctus 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list=ls()) 
set.seed(333) # reproducible
library(directlabels)
library(shiny) 
library(shinyjs)  
library(shinyWidgets)
library(shinythemes)  # more funky looking apps
library(shinyalert)
library(Hmisc)
library(rms)
library(ggplot2)
library(tidyverse)
library(shinycssloaders)
library(tvthemes)  # nice ggplot addition

options(max.print=1000000)    

fig.width6 <- 1100
fig.height6 <- 600
fig.width8 <- 1380
fig.height7 <- 770

## convenience functions
p0f <- function(x) {formatC(x, format="f", digits=0)}
p1f <- function(x) {formatC(x, format="f", digits=1)}
p2f <- function(x) {formatC(x, format="f", digits=2)}
p3f <- function(x) {formatC(x, format="f", digits=3)}
p4f <- function(x) {formatC(x, format="f", digits=4)}
p5f <- function(x) {formatC(x, format="f", digits=5)}

logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful

options(width=200)
options(scipen=999)

# range of independent variable
lowerV=0
upperV=10

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# function that does all the work!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


loq <- function (x, y, model, spec, print.plot=1, Xspec)  {
  
  # Define analysis models
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
  if (model %in% 12) {mod="Restricted cubic spline (rcs) 4 knots"} 
  
  # transformation of data for 12 models
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
  ty12 <- y;      tx12 <- x 
  
  # transform spec for prediction, note only where x is transformed, above
  if (model %in% 4 ) {Xspec <- 1/Xspec}
  if (model %in% 5 ) {Xspec <- 1/Xspec}
  if (model %in% 6 ) {Xspec <- log(Xspec)}
  if (model %in% 7 ) {Xspec <- log(Xspec)}
  if (model %in% 8 ) {Xspec <- sqrt(Xspec)}
  if (model %in% 10) {Xspec <- 1/Xspec}
  if (model %in% 11) {Xspec <- Xspec^2}
  
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
  # only used in text explanations
  tybar <- mean(y)     
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RCS MODEL
  if (model %in% 12) {       
    
    dat <- data.frame(cbind(y,x))
    ddist <<- datadist(dat)
    options(datadist='ddist')#
    
    f <- ols(y~rcs(x,4), dat)   
    
    # obtain the predictions 
    dat2 <- expand.grid(x=seq(min(x1),max(x1),.001))
    dat2 <- cbind(dat2, predict(f, dat2, se.fit=TRUE))
    dat2$lower <- dat2$linear.predictors - qt(0.975,n-4) * dat2$se.fit     # n-4 as we are using rcs 4 df are used up
    dat2$upper <- dat2$linear.predictors + qt(0.975,n-4) * dat2$se.fit
    
    # if (is.na(spec))  {spec=mean(dat2$linear.predictors, is.finite=TRUE)}
    if (is.na(spec))  {spec=mean(dat$y, is.finite=TRUE)}
    if (is.na(Xspec)) {Xspec=mean(dat$x, is.finite=TRUE)}
    
    #find nearest values to spec using brute force approach
    it <- which.min(abs(dat2$linear.predictors - spec))
    txpre<-dat2[it,]$x 
    
    it <- which.min(abs(dat2$lower - spec))
    txlow<-dat2[it,]$x 
    
    it <- which.min(abs(dat2$upper - spec))
    txup<-dat2[it,]$x 
    
    yspec <- spec
    
    limits <- sort(c(txlow,txup))
    txlow <- limits[1]
    txup <-  limits[2]
    
    # rcs will report the limit as the nearest value to spec if x value is beyond range, so lets report 999
    txpre <- ifelse((txpre >= max(dat$x) |(txpre <= min(dat$x)) ), 999, txpre)
    txlow <- ifelse((txlow >= max(dat$x) |(txlow <= min(dat$x) )), 999, txlow)
    txup <-  ifelse((txup >=  max(dat$x) |(txup  <= min(dat$x) )), 999, txup)
    
    rsd2 <-  anova(f)["ERROR","MS"]^.5
    dfs <- anova(f)["ERROR","d.f."]
    
    XXX  <- predict(f, Xspec, se.fit=TRUE) 
    XXX$L <- XXX$linear.predictors - qt(0.975,n-4) * XXX$se.fit     # n-4 as we are using rcs 4 df are used up
    XXX$U <- XXX$linear.predictors + qt(0.975,n-4) * XXX$se.fit
    pspec <- as.vector(unlist(XXX))
    pspec<- pspec[c(1,3,4)]
    if( pspec[3] < pspec[2] ) {pspec <- pspec[c(1,3,2)] }
    tp <- pspec  # capture this now to help explanation
    
    # predict again for plot, so we have predictions for the actual data
    xx <- predict(f, dat, se.fit=TRUE)
    xx$lower <- xx$linear.predictors - qt(0.975,n-4) * xx$se.fit     # n-4 as we are using rcs 4 df are used up
    xx$upper <- xx$linear.predictors + qt(0.975,n-4) * xx$se.fit
    xx <- as.data.frame(xx)
    
    # sum of squares of resuduals
    r <-   (y1-xx$linear.predictors) 
    r2 <-  (y1-xx$linear.predictors)^2
    ssr <- sum(r2, na.rm=T) 
    
    df2 <- (ssr/dfs)^.5
    
    foo <- as.data.frame(cbind(x=x1,obsy=y1, x2=x,y2=y,pred= xx$linear.predictors, p2a=xx$lower, p3=xx$upper, r=r, rr2=r2, rsd2=rsd2, dfs=dfs,ssr=ssr, df2=df2))
    foo <- foo[order(foo$obsy),]
    xx<- NULL
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# ALL OTHER MODELS
  } else {
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    f <- lm(y~x)  
    # run regression on the transformed data grab slope intercept
    intercept <- coef(f)[1][[1]]
    slope <- coef(f)[2][[1]]
    
    # obtain the predictions 
    p <- predict.lm(f, interval="confidence")
    
    # transform the predicted values & 95%CI back to original scale 
    if (model %in% c(2,7,10)) {p <- exp(p)} 
    if (model %in% c(3,5)  )  {p <- 1/p} 
    if (model %in% c(9)    )  {p <- p^2} 
    if (model %in% c(11)   )  {p <- p^.5} 
    
    if (is.na(Xspec)) {Xspec=mean(x1, is.finite=TRUE)
    
    if (model %in% c(4,5,10)) {Xspec <- 1/Xspec}
    if (model %in% c(6,7)  )  {Xspec <- log(Xspec)  }
    if (model %in% c(8)    )  {Xspec <- Xspec^.5}
    if (model %in% c(11)   )  {Xspec <- Xspec^2 }
    
    }

    tp <- pspec <- predict.lm(f, newdata=data.frame(x=Xspec), interval="confidence")  # tp will be used in explanationary text
    
    if (model %in% c(2,7,10)) {pspec <- exp(pspec)}
    if (model %in% c(3,5)  )  {pspec <- 1/pspec}
    if (model %in% c(9)    )  {pspec <- (pspec)^2}
    if (model %in% c(11)   )  {pspec <- pspec^.5}
    
    if(sum(is.nan(pspec) )==0) { ##if no invalid computation do this:
      if( pspec[3] < pspec[2] ) {pspec <- pspec[c(1,3,2)] }
    }
    
    # residuals original y and transformed back predicted values, residual sum of squares, this will be used to judge best model
    r <- (y1-p[,1]) 
    r2 <-  (y1-p[,1])^2
    ssr <- sum(r2, na.rm=T) 
    
    
    #######################################################################################################################
    # transform the specification that we will read back from, note only those in which y is transformed
    if (is.na(spec)) {spec=mean(pspec[1], is.finite=TRUE)}    
    
    yspec <- spec
    
    
  # else { 
      if (model %in% c(2,7,10)) {spec <- log(spec)}
      if (model %in% c(3,5)  )  {spec <- 1/spec}
      if (model %in% c(9)    )  {spec <- sqrt(spec)}
      if (model %in% c(11)   )  {spec <- spec ^2}
    #}
    
   
    
    # if (model %in% c(2,7,10)) {spec <- exp(spec)} 
    # if (model %in% c(3,5)  )  {spec <- 1/spec} 
    # if (model %in% c(9)    )  {spec <- (spec)^2} 
    # if (model %in% c(11)   )  {spec <- spec^.5}  
    # 
    tyspec <- spec 
    # grab the residual standard deviation
    rsd2 <- as.data.frame(anova(f))[2,3]^.5 
    dfs <- anova(f)["Residuals","Df"]
    
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
    
  
    
    
    
    ########################################################################################################################
    
    
    
    # if i don't do this X vertical spec line on plot for these models will not be located correctly
    if (model %in% c(4,5,10)) {Xspec <- 1/Xspec}
    if (model %in% c(6,7)  )  {Xspec <- exp(Xspec)  }
    if (model %in% c(8)    )  {Xspec <- Xspec^2}
    if (model %in% c(11)   )  {Xspec <- Xspec^.5 }
    
    # ensure order of limits is correct
    limits <- sort(c(txlow,txup))
    txlow <- limits[1]
    txup <-  limits[2]
    
    df2 <- (ssr/dfs)^.5
    
    foo <- data.frame(cbind(x=x1, obsy=y1, x2=x,y2=y,pred= p[,1], p2a=p[,2], p3=p[,3], r=r, rr2=r2, rsd2=rsd2, dfs=dfs,ssr=ssr,df2=df2))
    foo <- foo[order(foo$obsy),]
    
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # help with plotting
  lowerV=floor(min(foo$x)); upperV=ceiling(max(foo$x))
  ymin <- min(y)
  ymax <- max(y)
  ystep <- (ymax-ymin)/8
  ymin1 <-  ymin-ystep
  ymax1 <-  ymax+ystep
  
  # plot and present the estimated read back
  p1 <- ggplot(foo, aes(x=x,y=pred)) + 
    geom_line( ) +
    geom_ribbon(data=foo , aes(ymin= p2a,ymax= p3),alpha=0.2,   fill="green") +
    geom_point(data=foo, aes(x=x ,y=obsy), size=2, color='blue')  #+
  #    scale_x_continuous(limits = c(lowerV, upperV), breaks=seq(lowerV, upperV)) #, by=((upperV-lowerV)/10))) + # breaks = seq(0, 100, by = 20)
  # scale_y_continuous(limits = c(ymin1, ymax1))   
  
  p <- p1  + geom_hline(yintercept=yspec,  colour="#990000", linetype="dashed")
  p <- p   + geom_vline(xintercept=Xspec, colour="#008000", linetype="dashed")
  
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13,color="darkred"))
  p <- p + scale_color_manual(values=c("Red","blue"))
  p <- p + theme_bw()
  # p <- p + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
  p <- p + labs(x = "Independent variable", y = "Response") 
  
  p <- p +  theme(panel.background=element_blank(),
                  plot.title=element_text(size=16), 
                  plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                  legend.text=element_text(size=12),
                  legend.title=element_text(size=14),
                  axis.text.x  = element_text(size=12),
                  axis.text.y  = element_text(size=12),
                  axis.line.x = element_line(color="black"),
                  axis.line.y = element_line(color="black"),
                  plot.caption=element_text(hjust = 0, size = 12),
                  axis.title.y = element_text(size = rel(1.1), angle = 90),
                  axis.title.x = element_text(size = rel(1.1), angle = 00),
                  axis.title = element_text(size = 16, angle = 00)
  )   
  
  p <- p + labs(title = paste0("Fitted analysis model '",mod,"' with 95% confidence and raw data. N = ",length(!is.na(foo$x)),"\nResidual sum of squares = ", p2f(ssr),", residual standard deviation = ",p2f(df2)," \nPredict at input of ", 
                               p2f(Xspec) ,", the estimate of Y is ",
                               p4f(pspec[1])," with 95%CI: (", 
                               p4f(pspec[2]),", ",
                               p4f(pspec[3]),")",
                               "\nRead back at response of ", 
                               p2f(yspec) ,", the estimate of X is ",
                               p2f(txpre)," with 95%CI: (", 
                               p2f(txlow),", ",
                               p2f(txup),")",  
                               sep=" "),
                caption = paste0("If X and/or Y specification is missing, mean of the data is used, dashed lines")
  )  #   +
  
  # tried this package for plots themes but got errors
  #  theme_minimal() +
  #    theme(text = element_text(family = "Cinzel", size = 16),
  #       title = element_text(family = "Cinzel", size = 16)) -> targaryen
  
  if (print.plot==1) {print(p)}
  
  return(list(ssr=ssr,r=r, foo=foo, f=f, mod=mod, rsd2=rsd2, dfs=dfs , tybar=tybar, txbar=txbar, Xspec=Xspec, tp=tp, pspec=pspec ))
  
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
                 
                 h2("Exploring transformations and model fitting"), 
                 
                 h4("An independent variable is generated from a Uniform(0:10) distribution. The response derived from the user inputs and chosen
                 data generating mechanism. The data can be analysed using a selection of models. The best model fit can be selected ('Best scenario'), judged 
                 by the model with minimum sum of square of the residuals. A plot of the model fit is shown, predictions & read back are possible. Tab 2 shows summary stats, tab 3 assumptions and tab 6 upload your own data. See the Wiki."), 
                 
                 h3("  "), 
                 sidebarLayout(
                   
                   sidebarPanel( width=3,
                                 
                                 tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                 
                                 actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),  width =105 ,
                                              onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/LoQ/app.R', '_blank')"), 
                                 actionButton(inputId='ab1', label="R code",   icon = icon("th"),    width =105 ,
                                              onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/Rcode.R', '_blank')"),  
                                 actionButton("resample", "Simulate a new sample", width =180 ),
                                 
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
                                               div(h5(tags$span(style="color:blue", "Intercept"))), "10"),
                                     
                                     textInput('b', 
                                               div(h5(tags$span(style="color:blue", "Slope"))), ".1")
                                   ),
                                   
                                   splitLayout(
                                     textInput('sigma1', 
                                               div(h5(tags$span(style="color:blue", "Residual error"))), ".1"),
                                     
                                     textInput('spec', 
                                               div(h5(tags$span(style="color:blue", "Y specification"))), ""),
                                     
                                     textInput('Xspec', 
                                               div(h5(tags$span(style="color:blue", "X specification"))), "")
                                     
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
                                                       "model7", "model8", "model9",  "model10", "model11"),
                                     selected=c("model1")
                                   ),
                                   
                                   ###https://stackoverflow.com/questions/49616376/r-shiny-radiobuttons-how-to-change-the-colors-of-some-of-the-choices
                                   
                                   radioButtons(
                                     inputId = "ana",
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
                                       tags$span(style = "color:blue", "Square X and Y Y^2=a+X^2/b"),
                                       tags$span(style = "color:blue", "Restricted cubic spline 4 knots")
                                     ),
                                     choiceValues = c( "99","1", "2", "3",  "4", "5", "6",
                                                       "7", "8", "9",  "10", "11", "12")
                                     ,
                                     selected=c("1")
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
                               
                               tabPanel( "1 Model fitting",
                                         
                                         fluidRow(
                                           column(width = 6, offset = 0, style='padding:1px;',
                                                  shinycssloaders::withSpinner(
                                                    div(plotOutput("plot1",  width=fig.width8, height=fig.height7)),
                                                  ),
                                           ) ,
                                           
                                           
                                           fluidRow(
                                             column(width = 6, offset = 0, style='padding:1px;',
                                                    
                                             ))),#
                                         
                                         h4(paste("Figure 1 Model fit showing raw data and prediction with 95% confidence")),
                                         
                                         
                                         
                                         h4(htmlOutput("textWithNumber2",) ),
                                         width = 30  )     ,
                               
                               
                               
                               tabPanel("2 Summary stats", value=3, 
                                        h4(paste("X")),
                                        shinycssloaders::withSpinner(verbatimTextOutput("X"),type = 5),
                                        h4(paste("Y")),
                                        shinycssloaders::withSpinner(verbatimTextOutput("Y"),type = 5),
                                        h4(paste("Table 1 Summary statistics")),
                                        
                                        
                                        h4(htmlOutput("textWithNumber",) ),
                                        
                               ),
                               
                               tabPanel("3 Diagnostics", value=3, 
                                        
                                        shinycssloaders::withSpinner(
                                          div(plotOutput("diagnostics",  width=fig.width8, height=fig.height7)),
                                        ),
                                        
                                        h4("Figure 2. Three residual plots to check for absence of trends in central tendency and in variability"),
                                        p(strong("Upper left panel shows residuals versus fitted on the x-axis. 
                                              Bottom left panel is the QQ plot for checking normality of residuals from the OLS fit.
                                              Top right panel is the histogram for checking normality of residuals from the OLS fit with 
                                              ~N(mean=0, sd=OLS model sigma) curve and true SD superimposed.
                                              ")),
                                        
                               ),
                               
                               tabPanel("4 Models summary", value=3, 
                                        
                                        shinycssloaders::withSpinner(verbatimTextOutput("ssr"),type = 5),
                                        h4(paste("Table 2 Summary of model fits, note model sigma when data generating mechanism and analysis model coincide")),
                                        
                                        shinycssloaders::withSpinner(verbatimTextOutput("ssr2"),type = 5),
                                        h4(paste("Table 3 Models on transformed data")),
                               ),
                               
                               
                               tabPanel("5 Listing", value=3, 
                                        #   h4(paste("Table 4 Means")),
                                        #  shinycssloaders::withSpinner(verbatimTextOutput("dA"),type = 5),
                                        h4(paste("Table 4 Listing")),
                                        shinycssloaders::withSpinner(verbatimTextOutput("d2"),type = 5),
                               ),
                               
                               
                               
                               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NEW
                               tabPanel("6 User upload", fluid = TRUE, width = 4,
                                        
                                        h4(("Upload your own data for analysis. Requires 2 columns of numeric data. Select 'Header' 
                         if your data columns have names. 
                          The top two radio button options are to help load. Of course only the 'Analysis model' radio buttons are needed. Here is a link to example data (download a file and click 'Browse...' to locate and upload for the analysis):")) ,
                                        
                                        
                                        
                                        tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/data_example_1", tags$span(style="color:blue", "Example 1 data for analysis, has a header."),), 
                                        div(p(" ")),
                                        div(p(" ")),
                                        tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/data_example_2", tags$span(style="color:blue", "Example 2 data for analysis, has a header."),), 
                                        div(p(" ")),
                                        div(p(" ")),
                                        tags$a(href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/data_example_3", tags$span(style="color:blue", "Example 3 data for analysis, has a header."),), 
                                        div(p(" ")),
                                        div(p(" ")),
                                        
                                        
                                        
                                        
                                        sidebarLayout(
                                          
                                          # Sidebar panel for inputs ----
                                          sidebarPanel(width=2,
                                                       
                                                       # Input: Select a file ----
                                                       fileInput("file1", "Choose CSV File",
                                                                 multiple = TRUE,
                                                                 accept = c("text/csv",
                                                                            "text/comma-separated-values,text/plain",
                                                                            ".csv")),
                                                       
                                                       # Horizontal line ----
                                                       tags$hr(),
                                                       
                                                       # Input: Checkbox if file has header ----
                                                       checkboxInput("header", "Header", TRUE),
                                                       
                                                       # Input: Select separator ----
                                                       radioButtons("sep", "Separator",
                                                                    choices = c(Comma = ",",
                                                                                Semicolon = ";",
                                                                                Tab = "\t",
                                                                                Whitespace = ""),
                                                                    selected = ""),
                                                       
                                                       # Input: Select quotes ----
                                                       radioButtons("quote", "Quote",
                                                                    choices = c(None = "",
                                                                                "Double Quote" = '"',
                                                                                "Single Quote" = "'"),
                                                                    selected = ''),
                                                       
                                                       # Horizontal line ----
                                                       # tags$hr(),
                                                       
                                                       # Input: Select number of rows to display ----
                                                       # radioButtons("disp", "List all the data or first 6 rows only",
                                                       #              choices = c(Head = "head",
                                                       #                          All = "all"),
                                                       #              selected = "head"),
                                                       
                                                       # Horizontal line ----
                                                       # tags$hr(),
                                                       
                                                       # Input: Select number of rows to display ----
                                                       # radioButtons("what", "Output",
                                                       #              choices = c(Analysis = "Analysis",
                                                       #                          Plot = "plot"),
                                                       #              selected = "Analysis")
                                                       
                                          ),
                                          
                                          # Main panel for displaying outputs ----
                                          mainPanel(
                                            
                                            # Output: Data file ----
                                            
                                            div(plotOutput("plot2", width=fig.width6, height=fig.height6)),
                                            h4(paste("Figure 3. User uploaded data")),  
                                            h4(paste("X")),
                                            shinycssloaders::withSpinner(verbatimTextOutput("Xu"),type = 5),
                                            h4(paste("Y")),
                                            shinycssloaders::withSpinner(verbatimTextOutput("Yu"),type = 5),
                                            
                                            
                                          ),
                                        )
                               ) ,
                               tabPanel("7 User diagnostics", value=3, 
                                        shinycssloaders::withSpinner(
                                          div(plotOutput("diagnosticsu",  width=fig.width8, height=fig.height7)),
                                        ),
                                        
                                        h4("Figure 4. Three residual plots to check for absence of trends in central tendency and in variability"),
                                        p(strong("Upper left panel shows residuals versus fitted on the x-axis. 
                                              Bottom left panel is the QQ plot for checking normality of residuals from the OLS fit.
                                              Top right panel is the histogram for checking normality of residuals from the OLS fit with 
                                              ~N(mean=0, sd=OLS model sigma) curve and true SD superimposed.
                                              ")), 
                                        
                               ),
                               
                               tabPanel("8 User summary", value=3, 
                                        
                                        shinycssloaders::withSpinner(verbatimTextOutput("ssru"),type = 5),
                                        h4(paste("Table 5 Summary of model fits.")),
                                        
                                        shinycssloaders::withSpinner(verbatimTextOutput("ssr2u"),type = 5),
                                        h4(paste("Table 6 Models on transformed data")),
                               ),
                               
                               
                               tabPanel("9 User listing", value=3, 
                                        #h4(paste("Table 8 Means")),
                                        #shinycssloaders::withSpinner(verbatimTextOutput("d4u"),type = 5),
                                        h4(paste("Table 7 Listing")),
                                        shinycssloaders::withSpinner(verbatimTextOutput("d3"),type = 5),
                               ),
                               
                               tabPanel("10 Wiki", value=3, 
                                        
                                        tags$hr(),
                                        
                                        
                                        h4("The panel on the left contains the user inputs. The two empty specifications relate to Y and X, 
                                        when left empty the mean of the Y and the mean of X are used to read back X and predict Y respectively. "),
                                        h4("When we enter a Y specification this goes through an analysis transformation and we read back. 
                                        Now we have an X that we back transform. "),
                                        
                                        h4("If we enter a X or use the mean of X we predict Y and back transform."),
                                        h4("Tab 1 is the model fit based on the selcted radio buttons. Also find a step by step explanation of prediction and the transformation process.
                                        Tab 2 is simple summary stats of the original data (plotted in Figure 1) and a repeat of the step by step explanation of prediction and the transformation process.
                                        
                                        The Diagnostic tab 3 assesses the OLS model fit (using transformed data)."),
                                        
                                        h4("Tab 4 summarises briefly all analysis models. 'Back transformed' sigma is the residual error 
                                        calculated from the column 'ssr' in the listing divded by the degrees of freedom."),
                                        
                                        h4("Now we describe the contents of the data listing tab: X is uniform(N,0,10) and y is then derived from the user
                                        inputs in tandem with the selected 
                                           'Data generating mechanism'. Next we make a transformation of the data determined by the analysis method model. 
                                           See columns 'transformed x' and 'transformed y'. A simple OLS model is fit and prediction made. The prediction 
                                           is then transformed back to the original scale and this is presented next with associated 95% confidence interval. 
                                           This is what we plot along with X and Y mentioned in the first sentence. The residual is Y minus the prediction. 
                                           We square this next, the sum of this is used to judge the best model (see column 'ssr'). 
                                           The column labelled 'sigma' is the residual error from the OLS model. 
                                           The last column is the residual error calculated from the column 'ssr' divded by the degrees of freedom. 
                                           This is the residual that is presented with Figure 1. 
                                           Note the this sigma and the OLS sigma will conincide when the data generating mechanism and analysis 
                                           model coincide and will approximate the user input 'Residual error'"),
                                        h4("Tabs 6-9 
                                        allow the user to upload data, perform and evaluate an analysis. 
                                        Obviously there is no data generating mechanism so the top panel of radio buttons are not required and so have no impact."),
                                        #h4("Note the restricted cubic spine model will be the only model in which the fit may not pass through the mean of the data (defined when data the generating mechanism and analysis method coincide)."),
                                        tags$hr(),
                                        
                                        h4(paste("R code to quickly write a small data set to your desktop...ready to be uploaded to app")),
                                        h4("set.seed(124)"),
                                        h4("n <- sample(1:100,1)"),
                                        h4("x <- sort(rexp(n)*40)"),
                                        h4("y <- 10 + x^.5 + rnorm(n)"),
                                        h4("d <- cbind( x, y )"),
                                        h4("write.table(d, file = paste0( file.path(Sys.getenv('USERPROFILE'),'Desktop'),'/d.txt'),
                                                    sep = ' ', col.names = TRUE,
                                                    quote=FALSE, qmethod = 'double', row.names = FALSE)"),
                                        
                                        tags$hr(),
                                        h4(paste("To do:")),
                                        h4(paste("* Deal with read back when fitted and or limits cross multiple times the  y of interest (specification).")),
                                        h4(paste("* Convert main plot to plotly.")),
                                        
                                        
                               )
                               
                               ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                             )
                             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   )
                 ) 
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
server <- shinyServer(function(input, output   ) {
  
  shinyalert("Welcome! \nPlay with data transformations and fitting models",
             "Explore!", 
             type = "info")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This is where a new sample is instigated and inputs converted to numeric
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  random.sample <- reactive({
    
    foo <- input$resample
    
    a <- as.numeric(input$a)
    
    b <-as.numeric(input$b)
    
    sigma1 <- as.numeric(input$sigma1)
    
    N <- as.numeric(input$N)
    
    return(list(  
      a=a,
      b=b,
      sigma1=sigma1,
      N=N
    ))
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # GENERATE THE DATA
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dat <- reactive({
    
    sample <- random.sample()
    
    a=sample$a
    b=sample$b
    sigma1=sample$sigma1
    N=sample$N
    
    x <-  array(runif(N, lowerV, upperV))  # no negative values
    
    noise <-  rnorm(N,0, sigma1)
    
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
      y <-  ( a + (x^2)/b + noise)^.5
    }
    
    d <- cbind(x,y)
    
    return(list(  y=y, x=x, d=d))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Execute analysis
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  md <- reactive({
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    d <- dat()  # Get the  data
    y <- d$y
    x <- d$x
    
    ssr <- rep(NA,12)
    
    mdata <- list(NA)
    
    if (input$ana %in% 99) {
      
      for (j in 1:12) {
        
        res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0, Xspec=Xspec) # don't print
        ssr[j] <- res$ssr
        
      }
      
      model <- which(ssr==min(ssr)) 
      mdata <- res$foo
      res2 <- loq(x=x, y=y, model=model, spec= spec, print.plot=0,  Xspec=Xspec)  # run best model
      f=res2$f   
      mod<- res2$mod
      
    } else {
      
      res <- loq(x=x, y=y, model=as.numeric(input$ana), spec= spec,  Xspec=Xspec) 
      mdata <- res$foo
      model <- as.numeric(input$ana)
      f=res$f
      mod<- res$mod
      
    }
    
    return(list(  model=model, foo=mdata, f=f, mod=mod))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # collect for listing here
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  md2 <- reactive({
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    d <- dat()  # Get the data
    y <- d$y
    x <- d$x
    
    A <- array(NA, dim=c(12,6))
    M <- list(NA, dim=c(24,1))
    M <- array(NA, dim=c(24,1))
    
    for (j in 1:12) {
      
      res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0,  Xspec=Xspec) # don't print
      
      k <- j*2
      m <- k-1
      
      M[k] <- res$f 
      M[m] <- res$mod
      
      A[j,1] <- j          # model no
      A[j,2] <- res$mod    # model
      
      A[j,3] <- p0f(res$dfs)    # d.f.
      A[j,4] <- p4f(res$ssr)    # sum of squared residuals
      A[j,5] <- p4f(res$rsd2)   # sigma
      A[j,6] <- p4f(sqrt(res$ssr/res$dfs))
      
    }
    
    A <- data.frame( A[,c(1,2)],  apply(A[,c(3,4,5,6)],2, as.numeric))
    
    A <- plyr::arrange(A,A[,4])
    return(list(  ssr=A, M=M))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # X summary stats
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$X <- renderPrint({
    
    d <- md()$foo
    return(print(summary(d$x), digits=6))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Y summary stats
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$Y <- renderPrint({
    
    d <- md()$foo
    return(print(summary(d$obsy), digits=6))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # MAIN PLOT!
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$plot1 <- renderPlot({     
    
    model <- md()$model
    foo <- md()$foo
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    d <- dat()  # Get the  data
    y <- d$y
    x <- d$x
    
    loq(x= x, y= y, model=model, spec= spec, print.plot=1,  Xspec=Xspec) # print plot
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # EXPLANATION OF PREDICTION PROCESS, APPEARS ON TAB1 AND TAB 2
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$textWithNumber2 <- output$textWithNumber <- renderText({ 
    
    ## repeat plot1 code
    model <- md()$model
    foo <- md()$foo
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    d <- dat()  # Get the  data
    y <- d$y
    x <- d$x
    
    res  <- loq(x= x, y= y, model=model, spec= spec, print.plot=0,  Xspec=Xspec) # don't print plot
    
    m <-  as.numeric(gsub("[^0-9.-]", "", input$truth ))
    # a <-  as.numeric(gsub("[^0-9.-]", "", input$ana )) # user uploaded data does not have this, so will get an error
    
    a <- model
    
    if (m %in% 1 ) {mod="Linear Y=a+bX"}  
    if (m %in% 2 ) {mod="Exponential Y=exp(a+bX)"} 
    if (m %in% 3 ) {mod="Reciprocal-Y Y=1/(a+bX)"} 
    if (m %in% 4 ) {mod="Reciprocal-X Y=a+b/X"} 
    if (m %in% 5 ) {mod="Double Reciprocal Y=1/(a+b/X)"} 
    if (m %in% 6 ) {mod="Logarithmic-X Y=a+b(log(X))"} 
    if (m %in% 7 ) {mod="Multiplicative Y=aX^b"} 
    if (m %in% 8 ) {mod="Square Root-X Y=a+b(sqrt(X))"} 
    if (m %in% 9 ) {mod="Square Root-Y Y=(a+bX)^2"} 
    if (m %in% 10) {mod="S-curve Y=exp(a+b/X)"} 
    if (m %in% 11) {mod="Square X and Y Y^2=a+X^2/b"} 
    
    if (a %in% 1 ) {mod2="Linear Y=a+bX"}  
    if (a %in% 2 ) {mod2="Exponential Y=exp(a+bX)"} 
    if (a %in% 3 ) {mod2="Reciprocal-Y Y=1/(a+bX)"} 
    if (a %in% 4 ) {mod2="Reciprocal-X Y=a+b/X"} 
    if (a %in% 5 ) {mod2="Double Reciprocal Y=1/(a+b/X)"} 
    if (a %in% 6 ) {mod2="Logarithmic-X Y=a+b(log(X))"} 
    if (a %in% 7 ) {mod2="Multiplicative Y=aX^b"} 
    if (a %in% 8 ) {mod2="Square Root-X Y=a+b(sqrt(X))"} 
    if (a %in% 9 ) {mod2="Square Root-Y Y=(a+bX)^2"} 
    if (a %in% 10) {mod2="S-curve Y=exp(a+b/X)"} 
    if (a %in% 11) {mod2="Square X and Y Y^2=a+X^2/b"} 
    if (a %in% 12) {mod2="Restricted cubic spline (rcs) 4 knots"} 
    
    # to avoid complications pull this
    if (is.na(Xspec)) {Xspec=mean(res$txbar, is.finite=TRUE)}
    
    HTML(paste0( 
      
      br(), br(),
      
      "Explanation of the prediction and the transformation process:",
      
      br(), br(),
      
      "Step 1 After processing the user inputs and selected data generating mechanism "
      , tags$span(style="color:red",  mod) ,
      
      ", we have our data, with mean X "
      , tags$span(style="color:red",  p4f(mean(d$x))) ,
      " mean Y "
      , tags$span(style="color:red",  p4f(mean(d$y))) ,
      
      br(), br(),
      
      " Step 2 Transform this data according to analysis model "
      , tags$span(style="color:red",  mod2) ,
      
      ", now we have our transformed data, mean X "
      , tags$span(style="color:red",  p4f(mean(res$txbar))) ,
      " mean Y "
      , tags$span(style="color:red",  p4f(mean(res$tybar))),
      
      br(), br(),
      
      " Step 3 We also have our X specification, the mean of the potentially transformed x if no user X specification entered "
      , tags$span(style="color:red",  p4f(mean(Xspec))) ,
      
      br(), br(),
      
      " Step 4 Using our analysis model: "
      , tags$span(style="color:red",  mod2) ,
      " and X specification on the transformed data shown in step 3 let us predict Y = "
      , tags$span(style="color:red",  p4f(res$tp[1])) , 
      ", 95%CI ( "
      , tags$span(style="color:red",  p4f(res$tp[2])) ,
      ", "
      , tags$span(style="color:red",  p4f(res$tp[3])) ,
      " ) ",
      
      br(), br(),
      
      # " Step 5 Now let us back transform (if required) the specification shown in step 3 "
      # , tags$span(style="color:red",  p4f(res$Xspec)) , 
      # 
      # br(), br(),
      
      " Step 5 Finally let us back transform (if required) the prediction shown in step 4 "
      , tags$span(style="color:red",  p4f(res$pspec[1])) , 
      ", 95%CI ( "
      , tags$span(style="color:red",  p4f(res$pspec[2])) ,
      ", "
      , tags$span(style="color:red",  p4f(res$pspec[3])) ,
      " ) "
      
    ))
    
  })      
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DATA LISTING
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$d2 <- renderPrint({
    
    d <- md()$foo
    
    d <- plyr::arrange(d,x)
    names(d) <- c("x","y","transformed x","transformed y","prediction","lower 95%CI", "upper 95%CI", "residual","residual^2","sigma","df","sum of sq residuals (ssr)","ssr/df")
    return(print(d))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DATA LISTING AVERAGED
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$dA <- renderPrint({
    
    d <- md()$foo
    d <- plyr::arrange(d,x)
    dA <- apply(d,2,mean)   # could use colSums
    dA <- p5f(dA)
    
    
    dA <- lapply(dA, as.numeric)
    dA <- as.data.frame(dA)
    names(dA) <- c("x","y","transformed x","transformed y","prediction","lower 95%CI", "upper 95%CI", "residual","residual^2","sigma","df","sum of sq residuals (ssr)","ssr/df")
    return(print(dA, row.names=FALSE))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$d <- renderPrint({
    
    d <- dat()$d
    d <- as.data.frame(d)
    d <- plyr::arrange(d,x)
    
    return(print(d))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SUMMARY OF ALL MODELS
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$ssr <- renderPrint({
    
    d <- md2()$ssr
    d <- as.data.frame(d)
    names(d) <- c("model #","Model description", "d.f.", "Sum of square of residuals","Model Sigma" ,"Back transf. sigma")
    return(print(d, row.names = FALSE))
    
  })  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$ssr2 <- renderPrint({
    
    d <- md2()$M
    
    return(print(d))
    
  })  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DIAGNOSTICS
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$diagnostics<- renderPlot({         
    
    f <- md()$f
    mod <- md()$mod
    model <- md()$model
    
    sigma1 <- as.numeric(input$sigma1)
    residx <-  resid(f)
    fittedx <-    fitted(f)
    
    d <- cbind(residx, fittedx)
    d2 <- as.data.frame(d)
    
    # https://stackoverflow.com/questions/14200027/how-to-adjust-binwidth-in-ggplot2
    hist(residx,breaks="FD")
    breaks <- pretty(range(residx), n = nclass.FD(residx), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    
    yl <- ylab('Residuals')
    
    xl <- xlab("fitted")
    
    p1 <- ggplot(d2 , aes(x=fittedx , y=residx)) + geom_point (   colour="#69b3a2") + yl + xl
    
    p2 <- ggplot(d2 , aes(sample=residx )) + stat_qq(colour="#69b3a2") +
      geom_abline(intercept=mean(residx), slope=sd(residx)  ,  colour="black") +
      xlab('Normal theoretical quantiles')   + yl
    ggtitle( " ")
    
    library(gridExtra)
    library(grid)
    df <- data.frame(Residuals = residx)
    p3 <- ggplot(df, aes(x = Residuals)) + #stat_bin(bins = 30) +
      geom_histogram(aes(y =..density..), 
                     #breaks = seq(-50, 50, by = 2),
                     binwidth=bwidth,
                     colour = "black",
                     fill = "#69b3a2") +
      xlab('Residuals with superimposed sigma')   #+
    
    chk1 <-  as.numeric(gsub("[^0-9.-]", "", input$truth ))
    chk2 <-  as.numeric(gsub("[^0-9.-]", "", input$ana ))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    if (model %in% c(12) ) {
      
      std <-  f$stats["Sigma"][[1]]
      p3 <-  p3 + 
        stat_function(fun = dnorm, args = list(mean = 0, sd =  as.numeric(sigma1)        )) + 
        stat_function(fun = dnorm, args = list(mean = 0, sd =  std    ), col='red') 
      
      grid.arrange(p1,  p3, p2, ncol=2,
                   top = textGrob(paste0(" OLS model fit diagnostics, ",mod,", true sigma (black) ",as.numeric(sigma1)  ,", estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
    } else if (chk1==chk2) {  
      
      std <-  sigma(f) 
      p3 <-  p3 + 
        stat_function(fun = dnorm, args = list(mean = 0, sd =  as.numeric(sigma1)        )) + 
        stat_function(fun = dnorm, args = list(mean = 0, sd =  std    ), col='red') 
      
      grid.arrange(p1,  p3, p2, ncol=2,
                   top = textGrob(paste0(" OLS model fit diagnostics, ",mod,", true sigma (black) ",as.numeric(sigma1)  ,", estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
      # } else if (input$ana =="best") {  
      #   
      #   std <-  sigma(f) 
      #   p3 <-  p3 + 
      #     stat_function(fun = dnorm, args = list(mean = 0, sd =  std    ), col='red') 
      #   
      #   grid.arrange(p1,  p3, p2, ncol=2,
      #                top = textGrob(paste0(" OLS model fit diagnostics, ",mod,", estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
    } else {
      
      std <-  sigma(f) 
      p3 <- p3 + 
        stat_function(fun = dnorm, args = list(mean = 0, sd = std  ), col='red')  
      
      grid.arrange(p1,  p3, p2, ncol=2,
                   top = textGrob(paste0(" OLS model fit diagnostics, ",mod," estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FROM HERE ON WE ANALYSE USER LOADED DATA A LOT OF THE ABOVE CODE IS REUSED
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mdx <- reactive({
    
    df<-NULL
    req(input$file1)
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    d <- as.data.frame(df)
    
    d <- d[,c("x","y")]
    
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    
    y <- d$y
    x <- d$x
    
    ssr <- rep(NA,12)
    
    mdata <- list(NA)
    
    if (input$ana %in% 99) {
      
      for (j in 1:12) {
        
        res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0, Xspec=Xspec) # don't print
        ssr[j] <- res$ssr
        
      }
      
      model <- which(ssr==min(ssr, na.rm=TRUE)) 
      mdata <- res$foo
      res2 <- loq(x=x, y=y, model=model, spec= spec, print.plot=0,  Xspec=Xspec)  # run best model
      f=res2$f   
      mod<- res2$mod
      
    } else {
      
      res <- loq(x=x, y=y, model=as.numeric(input$ana), spec= spec,  Xspec=Xspec) 
      mdata <- res$foo
      model <- as.numeric(input$ana)
      f=res$f
      mod<- res$mod
      
    }
    
    return(list(  model=model, foo=mdata, f=f, mod=mod, x=x,y=y ))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  output$plot2 <- renderPlot({         
    
    model <- mdx()$model
    foo <- mdx()$fo2
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    
    y <- mdx()$y
    x <- mdx()$x
    
    loq(x= x, y= y, model=model, spec= spec, print.plot=1,  Xspec=Xspec) # print plot
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Repeat alot of above for user loaded data  
  md2u <- reactive({
    
    spec <- as.numeric(input$spec)
    Xspec <- as.numeric(input$Xspec)
    
    y <- mdx()$y
    x <- mdx()$x
    
    
    A <- array(NA, dim=c(12,6))
    M <- list(NA, dim=c(24,1))
    M <- array(NA, dim=c(24,1))
    
    for (j in 1:12) {
      
      res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0,  Xspec=Xspec) # don't print
      
      k <- j*2
      m <- k-1
      
      M[k] <- res$f 
      M[m] <- res$mod
      
      A[j,1] <- j          # model no
      A[j,2] <- res$mod    # model
      
      A[j,3] <- p0f(res$dfs)    # d.f.
      A[j,4] <- p4f(res$ssr)    # sum of squared residuals
      A[j,5] <- p4f(res$rsd2)   # sigma
      A[j,6] <- p4f(sqrt(res$ssr/res$dfs))
      
    }
    
    A <- data.frame( A[,c(1,2)],  apply(A[,c(3,4,5,6)],2, as.numeric))
    
    A <- plyr::arrange(A,A[,4])
    return(list(  ssr=A, M=M))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$ssr2u <- renderPrint({
    
    d <- md2u()$M
    
    return(print(d))
    
  })  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$ssru <- renderPrint({
    
    d <- md2u()$ssr
    d <- as.data.frame(d)
    names(d) <- c("model #","Model description", "d.f.", "Sum of square of residuals","Model Sigma" ,"Back transf. sigma")
    return(print(d, row.names = FALSE))
    
  })  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$Xu <- renderPrint({
    
    d <- mdx()$foo
    return(print(summary(d$x), digits=6))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$Yu <- renderPrint({
    
    d <- mdx()$foo
    return(print(summary(d$obsy), digits=6))
    
  }) 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$d3 <- renderPrint({
    
    d <- mdx()$foo
    
    d <- plyr::arrange(d,x)
    names(d) <- c("x","y","transformed x","transformed y","prediction","lower 95%CI", "upper 95%CI", "residual","residual^2","sigma","df","sum of sq residuals (ssr)","ssr/df")
    return(print(d))
    
    
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$d4u <- renderPrint({
    
    d <- mdx()$foo
    d <- plyr::arrange(d,x)
    dA <- apply(d,2,mean)
    dA <- p5f(dA)
    
    
    dA <- lapply(dA, as.numeric)
    dA <- as.data.frame(dA)
    names(dA) <- c("x","y","transformed x","transformed y","prediction","lower 95%CI", "upper 95%CI", "residual","residual^2","sigma","df","sum of sq residuals (ssr)","ssr/df")
    return(print(dA, row.names=FALSE))
    
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$diagnosticsu<- renderPlot({         
    
    f <- mdx()$f
    mod <- mdx()$mod
    model <- mdx()$model
    
    # sigma1 <- as.numeric(input$sigma1)
    residx <-  resid(f)
    fittedx <-    fitted(f)
    
    d <- cbind(residx, fittedx)
    d2 <- as.data.frame(d)
    
    # https://stackoverflow.com/questions/14200027/how-to-adjust-binwidth-in-ggplot2
    # hist(residx,breaks="FD")
    # breaks <- pretty(range(residx), n = nclass.FD(residx), min.n = 1)
    # bwidth <- breaks[2]-breaks[1]
    
    if (length(residx) < 200) {  bwidth <- length(residx) } else { bwidth <- 200 }
    
    yl <- ylab('Residuals')
    
    xl <- xlab("fitted")
    
    p1 <- ggplot(d2 , aes(x=fittedx , y=residx)) + geom_point (   colour="#69b3a2") + yl + xl
    
    p2 <- ggplot(d2 , aes(sample=residx )) + stat_qq(colour="#69b3a2") +
      geom_abline(intercept=mean(residx), slope=sd(residx)  ,  colour="black") +
      xlab('Normal theoretical quantiles')   + yl
    ggtitle( " ")
    
    library(gridExtra)
    library(grid)
    df <- data.frame(Residuals = residx)
    p3 <- ggplot(df, aes(x = Residuals)) + #stat_bin(bins = 30) +
      geom_histogram(aes(y =..density..), 
                     #breaks = seq(-50, 50, by = 2),
                     bins=bwidth,
                     colour = "black",
                     fill = "#69b3a2") +
      xlab('Residuals with superimposed sigma')   #+
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    if (model %in% c(12) ) {
      
      
      std <-  f$stats["Sigma"][[1]]
      p3 <-  p3 + 
        stat_function(fun = dnorm, args = list(mean = 0, sd =  std    ), col='red') 
      
      grid.arrange(p1,  p3, p2, ncol=2,
                   top = textGrob(paste0(" OLS model fit diagnostics, ",mod,", estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
    } else if (input$ana ==99) {  
      
      std <-  sigma(f) 
      p3 <-  p3 + 
        stat_function(fun = dnorm, args = list(mean = 0, sd =  std    ), col='red') 
      
      grid.arrange(p1,  p3, p2, ncol=2,
                   top = textGrob(paste0(" OLS model fit diagnostics, ",mod,", estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
    } else {
      
      std <-  sigma(f) 
      p3 <- p3 + 
        stat_function(fun = dnorm, args = list(mean = 0, sd = std  ), col='red')  
      
      grid.arrange(p1,  p3, p2, ncol=2,
                   top = textGrob(paste0(" OLS model fit diagnostics, ",mod," estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
      
    }
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)