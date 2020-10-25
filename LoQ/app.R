#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# White-margined Burrower Bug nymph - Sehirus cinctus

        rm(list=ls()) 
        set.seed(333) # reproducible
        library(directlabels)
        library(shiny) 
        library(shinyjs)  #refresh'
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
        
        fig.width8 <- 1380
        fig.height7 <- 770

        ## convenience functions
        p0f <- function(x) {formatC(x, format="f", digits=0)}
        p1f <- function(x) {formatC(x, format="f", digits=1)}
        p2f <- function(x) {formatC(x, format="f", digits=2)}
        p3f <- function(x) {formatC(x, format="f", digits=3)}
        p4f <- function(x) {formatC(x, format="f", digits=4)}
        p5f <- function(x) {formatC(x, format="f", digits=5)}
       # p2f <- function(x) {formatC(x, format="f", digits=4)}
        
        logit <- function(p) log(1/(1/p-1))
        expit <- function(x) 1/(1/exp(x) + 1)
        inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
        is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
        
        options(width=200)
        options(scipen=999)
       # w=4  # line type
       # ww=3 # line thickness
       # wz=1 
        
        # range of independent variable
        lowerV=0
        upperV=10
        
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

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
    if (model %in% 12) {mod="Restricted cubic spline 4 knots"} 
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
    ty12 <- y;      tx12 <- x ###############NEW
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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (model %in% 12) {       
      
      dat <- data.frame(cbind(y,x))
      ddist <<- datadist(dat)
      options(datadist='ddist')#
      
      f <- ols(y~rcs(x,4), dat)   
      
      # obtain the predictions 
      dat2 <- expand.grid(x=seq(lowerV,upperV,.001))
      dat2 <- cbind(dat2, predict(f, dat2, se.fit=TRUE))
      dat2$lower <- dat2$linear.predictors - qt(0.975,n-4) * dat2$se.fit     # n-4 as we are using rcs 4 df are used up
      dat2$upper <- dat2$linear.predictors + qt(0.975,n-4) * dat2$se.fit
      
      #find nearest values to spec using brute force approach
      it <- which.min(abs(dat2$linear.predictors - spec))
      txpre<-dat2[it,]$x 
     
      it <- which.min(abs(dat2$lower - spec))
      txlow<-dat2[it,]$x 
      
      it <- which.min(abs(dat2$upper - spec))
      txup<-dat2[it,]$x 
      
      limits <- sort(c(txlow,txup))
      txlow <- limits[1]
      txup <-  limits[2]
      
      # rcs will report the limit as the nearest value to spec if x value is beyond range, so lets report NA 
      txpre <- ifelse((txpre >= upperV )|(txpre <= lowerV ), 999, txpre)
      txlow <- ifelse((txlow >= upperV )|(txlow <= lowerV ), 999, txlow)
      txup <-  ifelse((txup >=  upperV )|(txup  <= lowerV ), 999, txup)
      
      
      rsd2 <-  anova(f)["ERROR","MS"]^.5
   
      
      
      # predict again for plot, so we have predictions for the actual data
      xx <- predict(f, dat, se.fit=TRUE)
      xx$lower <- xx$linear.predictors - qt(0.975,n-4) * xx$se.fit     # n-4 as we are using rcs 4 df are used up
      xx$upper <- xx$linear.predictors + qt(0.975,n-4) * xx$se.fit
      xx <- as.data.frame(xx)
     
      
      # need original length
      
      r <-   (y1-xx$linear.predictors) 
      r2 <-  (y1-xx$linear.predictors)^2
      ssr <- sum(r2, na.rm=T) 
      
      
      foo <- as.data.frame(cbind(x=x1,   obsy=y1, pred= xx$linear.predictors, p2a=xx$lower, p3=xx$upper, r=r, rr2=r2))
      foo <- foo[order(foo$obsy),]
      xx<- NULL
      
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       } else {
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
    f <- lm(y~x)  
    # run regression on the transformed data grab slope intercept
    intercept <- coef(f)[1][[1]]
    slope <- coef(f)[2][[1]]
    
    # R2 on the transformed x and y,  an idea but not used
    # v1 <- unlist(cor.test(x,y)$estimate^2)[1][[1]]   
    
    # obtain the predictions 
    
    p <- predict.lm(f, interval="confidence")
    
    # transform the predicted values & 95%CI back to original scale 
    
    if (model %in% c(2,7,10)) {p <- exp(p)} 
    if (model %in% c(3,5)  )  {p <- 1/p} 
    if (model %in% c(9)    )  {p <- p^2} 
    if (model %in% c(11)   )  {p <- p^.5} 
    
    # calculate residual squared 
    # residuals original y and transformed back predicted values 
    r <- (y1-p[,1]) 
    r2 <-  (y1-p[,1])^2
    ssr <- sum(r2, na.rm=T) 
    
    # R2 on the original x and transformed back predicted y, an idea but not used
    # v2 <- unlist(cor.test(x1,p[,1])$estimate^2)[1][[1]]
    
    # residual sum of squares, this will be used to judge best model
    
    #ssr <- sum(r, na.rm=T) 
    
    # transform the response that we will read back
    
    tyspec <- spec 
    if (model %in% c(2,7,10)) {tyspec <- log(tyspec)} 
    if (model %in% c(3,5)  )  {tyspec <- 1/tyspec} 
    if (model %in% c(9)    )  {tyspec <- sqrt(tyspec)} 
    if (model %in% c(11)   )  {tyspec <- tyspec ^2}  
    
    # grab the residual standard deviation
    
    rsd2 <- as.data.frame(anova(f))[2,3]^.5  
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
    # put all pertinent data together, original data and predicted with 95%CI
    
    limits <- sort(c(txlow,txup))
    txlow <- limits[1]
    txup <-  limits[2]
    
    foo <- data.frame(cbind(x=x1, obsy=y1, pred= p[,1], p2a=p[,2], p3=p[,3], r=r, rr2=r2))
    foo <- foo[order(foo$obsy),]
    
    
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    
 # plot and present the estimated read back
    
    p1 <- ggplot(foo, aes(x=x,y=pred)) + 
        geom_line( ) +
        geom_ribbon(data=foo , aes(ymin= p2a,ymax= p3),alpha=0.2,   fill="green") +
        geom_point(data=foo, aes(x=x ,y=obsy), size=2, color='blue')  
    scale_x_continuous(limits = c(xmin1, xmax1))
    scale_y_continuous(limits = c(ymin1, ymax1))  # this was x?
    
    p <- p1  + geom_hline(yintercept=spec, colour="#990000", linetype="dashed")
    
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13,color="darkred"))
    p <- p + scale_color_manual(values=c("Red","blue"))
    p <- p + theme_bw()
    p <- p + xlab('Independent variable') 
   p <- p + ylab('Dependent variable')
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
              plot.caption=element_text(hjust = 0, size = 11),
              axis.title.y = element_text(size = rel(1.1), angle = 90),
              axis.title.x = element_text(size = rel(1.1), angle = 00),
              axis.title = element_text(size = 16, angle = 00)
        )   
        
    p <- p + labs(title = paste("Figure of the fitted model '",mod,"' with 95% confidence and raw data. \nExploration of model fitting at response (spec) of", 
                                spec ,", the estimate of x is",
                                p2f(txpre),"and 95% CI: (", 
                                p2f(txlow),",",
                                p2f(txup),")","\nResidual sum of squares", p2f(ssr),", Residual standard deviation",p2f(rsd2),
                                   sep=" "),
                  #  subtitle = paste("Model for the curve #",model," ",mod,""),
                    caption = paste0("We are interested in the independent variable value when y = ",p4f(spec),"")
    )  #   +
 
  #  theme_minimal() +
  #    theme(text = element_text(family = "Cinzel", size = 16),
     #       title = element_text(family = "Cinzel", size = 16)) -> targaryen
 
        
    if (print.plot==1) {print(p)}

    return(list(ssr=ssr,r=r, foo=foo, f=f, mod=mod))
    
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
                 
                 h4("An independent variable is generated using a uniform(0:10) distribution. Using the user inputs a response is derived from a selection 
                 of data generating mechanisms. The data can then be analysed using a selection of models. The best model fit can be selected ('Best scenario' button), judged 
                 by the model with the minimum sum of square of the residuals. A plot of the model fit is presented, then on tab 2 model assumptions are evaluated. 
                    The final tab presents a listing of the data."), 
                 
                 h3("  "), 
                 sidebarLayout(
                     
                     sidebarPanel( width=3 ,
                                   
                                   tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                   
                                   actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/LoQ/app.R', '_blank')"), 
                                   actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                                onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/Rcode.R', '_blank')"),  
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
                                           
                                           textInput('sigma1', 
                                                     div(h5(tags$span(style="color:blue", "sigma1"))), "2"),
                                           
                                           textInput('spec', 
                                                     div(h5(tags$span(style="color:blue", "spec"))), "0")
                                         
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
                                           choiceValues = c( "best","1", "2", "3",  "4", "5", "6",
                                                             "7", "8", "9",  "10", "11", "12")
                                           ,
                                           selected=c("11")
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

                                             h4(paste("Figure 1 Model Fit")),
                                             
                                             width = 30 )     ,
          
                                   tabPanel("2 Diagnostics", value=3, 
                                            
                                             shinycssloaders::withSpinner(
                                                div(plotOutput("diagnostics",  width=fig.width8, height=fig.height7)),
                                            ),

                                            h4("Figure 2. Three residual plots to check for absence of trends in central tendency and in variability"),
                                             p(strong("Upper left panel shows residuals versus fitted on the x-axis. 
                                              Bottom left panel is the QQ plot for checking normality of residuals from the LS fit.
                                              Top right panel is the histogram for checking normality of residuals from the LS fit with 
                                              ~N(mean=0, sd=GLS model sigma) curve and true SD superimposed.
                                              ")),
                  
                                   ),
                                            
                                            tabPanel("3 Data listing", value=3, 
                                                     shinycssloaders::withSpinner(verbatimTextOutput("d2"),type = 5),
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
    
    shinyalert("Welcome! \nPlay with data transformations and fitting models",
               "Have fun!", 
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
        
        spec <- as.numeric(input$spec)
        
        return(list(  
            a=a,
            b=b,
            sigma1=sigma1,
            N=N,
            spec=spec
        ))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     dat <- reactive({
        
       # generate data , we don
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
     md <- reactive({
    
         spec <- as.numeric(input$spec)
         
         d <- dat()  # Get the  data
         y <- d$y
         x <- d$x
         
         ssr <- rep(NA,12)
         
         mdata <- list(NA)
         
         if (input$ana %in% "best") {
             
             for (j in 1:12) {
                 
                 res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0) # don't print
                 ssr[j] <- res$ssr
                 
             }
             
              model <- which(ssr==min(ssr)) 
              mdata <- res$foo
              res2 <- loq(x=x, y=y, model=model, spec= spec, print.plot=0)  # run best model
              f=res2$f   
              mod<- res2$mod
              
         } else {
             
             res <- loq(x=x, y=y, model=as.numeric(input$ana), spec= spec) 
             mdata <- res$foo
             model <- as.numeric(input$ana)
             f=res$f
             mod<- res$mod
             
         }
    
         return(list(  model=model, foo=mdata, f=f, mod=mod))
    
     })
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     output$plot1 <- renderPlot({         #standard errors
        
        model <- md()$model
        foo <- md()$foo
        
        spec <- as.numeric(input$spec)
                
        d <- dat()  # Get the  data
        y <- d$y
        x <- d$x
        
        loq(x= x, y= y, model=model, spec= spec, print.plot=1) # print plot
                 
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     output$d2 <- renderPrint({
         
         d <- md()$foo
         
       d <- plyr::arrange(d,x)
         names(d) <- c("x","y","prediction","lower 95%CI", "upper 95%CI", "residual","residual^2")
         return(print(d))
         
     }) 
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         output$d <- renderPrint({
        
       d <- dat()$d
       d <- as.data.frame(d)
       d <- plyr::arrange(d,x)
       
       return(print(d))
        
    }) 
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         output$diagnostics<- renderPlot({         
             
        
             f <- md()$f
             mod <- md()$mod
             model <- md()$model
             
             sigma1 <- as.numeric(input$sigma1)
             residx <-  resid(f)
             fittedx <-    fitted(f)
             
             d <- cbind(residx, fittedx)
             d2 <- as.data.frame(d)
             
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
             p3 <- ggplot(df, aes(x = Residuals)) +
                 geom_histogram(aes(y =..density..),
                                #breaks = seq(-50, 50, by = 2),
                                colour = "black",
                                fill = "#69b3a2") +
                 xlab('Residuals with superimposed sigma')   #+
               
               
               if (model %in% 12) {         
                 std <-  f$stats["Sigma"][[1]]
                 p3 <-  p3 + 
                   stat_function(fun = dnorm, args = list(mean = 0, sd =  as.numeric(sigma1)        )) + 
                   stat_function(fun = dnorm, args = list(mean = 0, sd =  std    ), col='red') 
                 
                 
              
                 
               } else {
                 std <-  sigma(f) 
                 p3 <- p3 + 
                   stat_function(fun = dnorm, args = list(mean = 0, sd = as.numeric(sigma1)   )) + 
                   stat_function(fun = dnorm, args = list(mean = 0, sd = std  ), col='red')  
                 
               }
             
             grid.arrange(p1,  p3, p2, ncol=2,
                          top = textGrob(paste0(" LS model fit diagnostics, ",mod,", true sigma (black) ",as.numeric(sigma1)  ,", estimated sigma (red) ", p4f(std),""),gp=gpar(fontsize=20,font=3)))
             
               
             #          stat_function(fun = dnorm, args = list(mean = 0, sd = as.numeric(input$sigma)   )) +# sigma(f)
             # stat_function(fun = dnorm, args = list(mean = 0, sd = sigma(f)    ), col='red') # 
             # 
             # grid.arrange(p1,  p3, p2, ncol=2,
             #     top = textGrob(paste0(" LS model fit diagnostics, ",mod,", true sigma (black) ",as.numeric(input$sigma) ,", estimated sigma (red) ", p4(sigma(f)),""),gp=gpar(fontsize=20,font=3)))
             # 
             
         })
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

})
 
# Run the application 
shinyApp(ui = ui, server = server)