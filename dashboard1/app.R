# load the required packages
library(shiny)
require(shinydashboard)
library(ggplot2)
library(dplyr)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



set.seed(333) # reproducible
library(directlabels)
library(shiny) 
# library(shinyjs)  
# library(shinyWidgets)
# library(shinythemes)  # more funky looking apps
library(shinyalert)
library(Hmisc)
library(rms)
library(ggplot2)
library(tidyverse)
#library(shinycssloaders)
#library(tvthemes)  # nice ggplot addition
library(scales) # For the trans_format function
options(max.print=1000000)    

fig.width6 <- 1100
fig.height6 <- 600
fig.width8 <- 1380
fig.height7 <- 770

## convenience functions, see next code for better approach
# p0f <- function(x) {formatC(x, format="f", digits=0)}
# p1f <- function(x) {formatC(x, format="f", digits=1)}
# p2f <- function(x) {formatC(x, format="f", digits=2)}
# p3f <- function(x) {formatC(x, format="f", digits=3)}
# p4f <- function(x) {formatC(x, format="f", digits=4)}
# p5f <- function(x) {formatC(x, format="f", digits=5)}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# trying new formatting approach
formatz <- function(x){
    
    if (!is.na(x)  ) {
        
        if (abs(x) > 1000000) {
            
            formatC(x, format="f", digits=0,  big.mark=",")
            
        } else if (abs(x) > 1000 & abs(x) <=1000000) {
            
            formatC(x, format="f", digits=2,  big.mark=",")
            
        } else if (abs(x) >= 1 & abs(x) <=1000) {
            
            formatC(x, format="f", digits=2,  big.mark=",")
            
        } else if (abs(x) < 1) {
            
            formatC(x, format="f", digits=7,  big.mark=",")
            
        } else {
            
            formatC(x, format="f", digits=2,  big.mark=",")
            
        }
        
    }
    
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful

options(width=200)
options(scipen=999)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# function to create minor lines to match log tick values https://r-graphics.org/recipe-axes-axis-log-ticks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

breaks_5log10 <- function(x) {
    low <- floor(log10(min(x)/5))
    high <- ceiling(log10(max(x)/5))
    
    c(2:9 %o% 10^(low:high))
}

breaks_log10 <- function(x) {
    low <- floor(log10(min(x)))
    high <- ceiling(log10(max(x)))
    
    10^(seq.int(low, high))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# function to create ticks, not used
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# https://stackoverflow.com/questions/46412250/ggplot2-displaying-unlabeled-tick-marks-between-labeled-tick-marks
# breaks <- seq(lowerV, upperV, 10)
# breaks <- seq(lowerV, upperV, 1)
# labels <- as.character(breaks)
# labels[!(breaks %% 5 == 0)] <- ''
# tick.sizes <- rep(.5, length(breaks))
# tick.sizes[(breaks %% lowerV == 0)] <- 1

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
    
    # transformation of data for 12 models, create all and select below 
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
    
  
    
    if (model %in% c(4,5,10)) {Xspec <- 1/Xspec}
    if (model %in% c(6,7)  )  {Xspec <- log(Xspec)  }
    if (model %in% c(8)    )  {Xspec <- Xspec^.5}
    if (model %in% c(11)   )  {Xspec <- Xspec^2 }
    
    # save the original data
    x1 <- x
    y1 <- y
    
    # assign the selected transformed data to x and y 
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
        
        f <- ols(y~rcs(x,4), dat)   # OLS
        
        # obtain the predictions 
        dat2 <- expand.grid(x=seq(min(x1),max(x1),.001))
        dat2 <- cbind(dat2, predict(f, dat2, se.fit=TRUE))
        dat2$lower <- dat2$linear.predictors - qt(0.975,n-4) * dat2$se.fit     # n-4 as we are using rcs 4 df are used up
        dat2$upper <- dat2$linear.predictors + qt(0.975,n-4) * dat2$se.fit
        
        # assign specs if not entered by user
        if (is.na(spec))  {spec=mean(dat$y, is.finite=TRUE)}
        if (is.na(Xspec)) {Xspec=mean(dat$x, is.finite=TRUE)}
        
        # find nearest values to spec using brute force approach
        it <- which.min(abs(dat2$linear.predictors - spec))
        txpre<-dat2[it,]$x 
        
        it <- which.min(abs(dat2$lower - spec))
        txlow<-dat2[it,]$x 
        
        it <- which.min(abs(dat2$upper - spec))
        txup<-dat2[it,]$x 
        
        # assign spec to object
        yspec <- spec
        
        # ensure order is correct
        limits <- sort(c(txlow,txup))
        txlow <- limits[1]
        txup <-  limits[2]
        
        # rcs will report the limit as the nearest value to spec if x value is beyond range, so lets report 999
        txpre <- ifelse((txpre >= max(dat$x) |(txpre <= min(dat$x)) ), 999, txpre)
        txlow <- ifelse((txlow >= max(dat$x) |(txlow <= min(dat$x) )), 999, txlow)
        txup <-  ifelse((txup >=  max(dat$x) |(txup  <= min(dat$x) )), 999, txup)
        
        rsd2 <-  anova(f)["ERROR","MS"]^.5
        dfs <-   anova(f)["ERROR","d.f."]
        
        # get a prediction for the x spec
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
        
        # sum of squares of residuals
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
        
        # run regression on the transformed data grab slope intercept
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
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Now the X specification
        if (is.na(Xspec)) {
            
            Xspec=mean(x1, is.finite=TRUE)
            
            if (model %in% c(4,5,10)) {Xspec <- 1/Xspec}
            if (model %in% c(6,7)  )  {Xspec <- log(Xspec)  }
            if (model %in% c(8)    )  {Xspec <- Xspec^.5}
            if (model %in% c(11)   )  {Xspec <- Xspec^2 }
        }
        
        tp <- pspec <- predict.lm(f, newdata=data.frame(x=Xspec), interval="confidence")  # tp will be used in explanationary text
        
        # transform back
        if (model %in% c(2,7,10)) {pspec <- exp(pspec)}
        if (model %in% c(3,5)  )  {pspec <- 1/pspec}
        if (model %in% c(9)    )  {pspec <- (pspec)^2}
        if (model %in% c(11)   )  {pspec <- pspec^.5}
        
        if(sum(is.nan(pspec) )==0) { ##if no invalid computation do this:
            if( pspec[3] < pspec[2] ) {pspec <- pspec[c(1,3,2)] }
        }
        
        # p created above, residuals using original y and transformed back predicted values, residual sum of squares, this will be used to judge best model
        r <- (y1-p[,1]) 
        r2 <-  (y1-p[,1])^2
        ssr <- sum(r2, na.rm=T) 
        
        #######################################################################################################################
        # transform the specification that we will read back from, note only those in which y is transformed
        if (is.na(spec)) {spec=mean(pspec[1], is.finite=TRUE)}    
        
        yspec <- spec
        
        if (model %in% c(2,7,10)) {spec <- log(spec)}
        if (model %in% c(3,5)  )  {spec <- 1/spec}
        if (model %in% c(9)    )  {spec <- sqrt(spec)}
        if (model %in% c(11)   )  {spec <- spec ^2}
        
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
    
    lowerV <- floor(min(foo$x)/1)*1
    upperV <- ceiling(max(foo$x)/.01)*.01
    ymin <- min(y)
    ymax <- max(y)
    ystep <- (ymax-ymin)/8
    ymin1 <-  ymin-ystep
    ymax1 <-  ymax+ystep
    
    breaks <- seq(lowerV, upperV, 1)
    breaks <- seq(lowerV, upperV, 1)
    labels <- as.character(breaks)
    labels[!(breaks %% 5 == 0)] <- ''
    tick.sizes <- rep(.5, length(breaks))
    tick.sizes[(breaks %% lowerV == 0)] <- 1
    
    
    # plot and present the estimated read back
    p1 <- ggplot(foo, aes(x=x,y=pred)) +  
        geom_line( ) +
        geom_ribbon(data=foo , aes(ymin= p2a,ymax= p3),alpha=0.2,   fill="green") +
        geom_point(data=foo, aes(x=x ,y=obsy), size=2, color='blue')  #+
    #    scale_x_continuous(limits = c(lowerV, upperV), breaks=seq(lowerV, upperV)) #, by=((upperV-lowerV)/10))) + # breaks = seq(0, 100, by = 20)
    # scale_y_continuous(limits = c(ymin1, ymax1))   
    
    p <- p1  + geom_hline(yintercept=yspec,  colour="#990000", linetype="dashed")
    p <- p   + geom_vline(xintercept=Xspec, colour="#008000", linetype="dashed")
    
    #p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13,color="darkred"))
    p <- p + scale_color_manual(values=c("Red","blue"))
    p <- p + theme_bw()
    # p <- p + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    # p <- p + scale_y_continuous(trans="log", breaks=c(0.001, 0.1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6))
    
    # p <- p + scale_y_continuous(minor_breaks = waiver()        )
    
    #p <- p + scale_x_continuous(breaks = breaks, labels = labels, limits = c(lowerV, upperV)) 
    #p <- p + scale_x_continuous(breaks = waiver() , labels = waiver() ) 
    p <- p + scale_y_continuous(labels = scales::comma) 
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
                    axis.title = element_text(size = 16, angle = 00),
                    #panel.grid.minor = element_line(colour="gainsboro", size=0.3 , linetype = 'solid'),
                    # panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "gray88"),
                    #  panel.grid = element_blank(), axis.ticks.x = element_line(size = tick.sizes)
                    panel.grid.minor.x = element_line( size=0.5 ,  linetype = 'solid', colour = "gray93"),
                    panel.grid.minor.y = element_line( size=0.5 ,  linetype = 'solid', colour = "gray88"),
                    panel.grid.major = element_line(size =0.5,   linetype = 'solid', colour = "gray88"),            
                    panel.grid = element_blank(), axis.ticks.x = element_line(size = tick.sizes)
                    
                    
                    
    )   
    
    p <- p + labs(title = paste0("Fitted analysis model '",mod,"' with 95% confidence and raw data. N = ",length(!is.na(foo$x)),"\nResidual sum of squares = ", formatz(ssr),", residual standard deviation = ",formatz(df2)," \nPredict at input of ", 
                                 formatz(Xspec) ,", the estimate of Y is ",
                                 formatz(pspec[1])," with 95%CI: (", 
                                 formatz(pspec[2]),", ",
                                 formatz(pspec[3]),")",
                                 "\nRead back at response of ", 
                                 formatz(yspec) ,", the estimate of X is ",
                                 formatz(txpre)," with 95%CI: (", 
                                 formatz(txlow),", ",
                                 formatz(txup),")",  
                                 sep=" "),
                  caption = paste0("If X specification is missing, the mean of X is used. If Y specification is missing the prediction of X used as the Y value to read back from (dashed lines).")
    )  #   +
    
    # tried this package for plots themes but got errors
    #  theme_minimal() +
    #    theme(text = element_text(family = "Cinzel", size = 16),
    #       title = element_text(family = "Cinzel", size = 16)) -> targaryen
    
    if (print.plot==1) {print(p)}
    
    return(list(ssr=ssr,r=r, foo=foo, f=f, mod=mod, rsd2=rsd2, dfs=dfs , tybar=tybar, txbar=txbar, Xspec=Xspec, tp=tp, pspec=pspec ))
    
} 
# has log transformation of y axis
loq1 <- function (x, y, model, spec, print.plot=1, Xspec)  {
    
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
        
        if (is.na(Xspec)) {
            
            Xspec=mean(x1, is.finite=TRUE)
            
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
        
        if (model %in% c(2,7,10)) {spec <- log(spec)}
        if (model %in% c(3,5)  )  {spec <- 1/spec}
        if (model %in% c(9)    )  {spec <- sqrt(spec)}
        if (model %in% c(11)   )  {spec <- spec ^2}
        
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
    
    lowerV <- floor(min(foo$x)/1)*1
    upperV <- ceiling(max(foo$x)/.01)*.01
    ymin <- min(y)
    ymax <- max(y)
    ystep <- (ymax-ymin)/8
    ymin1 <-  ymin-ystep
    ymax1 <-  ymax+ystep
    
    breaks <- seq(lowerV, upperV, 10)
    breaks <- seq(lowerV, upperV, 1)
    labels <- as.character(breaks)
    labels[!(breaks %% 5 == 0)] <- ''
    tick.sizes <- rep(.5, length(breaks))
    tick.sizes[(breaks %% lowerV == 0)] <- 1
    
    
    
    # plot and present the estimated read back
    p1 <- ggplot(foo, aes(x=x,y=pred)) + 
        geom_line( ) +
        geom_ribbon(data=foo , aes(ymin= p2a,ymax= p3),alpha=0.2,   fill="green") +
        geom_point(data=foo, aes(x=x ,y=obsy), size=2, color='blue')  #+
    #    scale_x_continuous(limits = c(lowerV, upperV), breaks=seq(lowerV, upperV)) #, by=((upperV-lowerV)/10))) + # breaks = seq(0, 100, by = 20)
    # scale_y_continuous(limits = c(ymin1, ymax1))   
    
    p <- p1  + geom_hline(yintercept=yspec,  colour="#990000", linetype="dashed")
    p <- p   + geom_vline(xintercept=Xspec,  colour="#008000", linetype="dashed")
    
    # p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=13,color="darkred"))
    p <- p + scale_color_manual(values=c("Red","blue"))
    p <- p + theme_bw()
    
    # p <- p + scale_y_continuous(trans="log10", breaks=lseq(),
    #                          minor_breaks = breaks_5log10,
    #                             )  # see lseq function
    # 
    p <- p + scale_y_log10(breaks = breaks_log10,
                           minor_breaks = breaks_5log10,
                           labels = trans_format(log10, math_format(10^.x))) 
    
    p <- p + annotation_logticks(sides = "lr")
    
    #p <- p + scale_x_continuous(breaks = breaks, labels = labels, limits = c(0,100)) 
    #p <- p + scale_x_continuous(breaks = waiver() , labels = waiver() ) 
    # memory error if I use this!
    # p <- p +scale_x_log10(
    #   breaks = scales::trans_breaks("log10", function(x) 10^x),
    #   labels = scales::trans_format("log10", scales::math_format(10^.x))
    # )
    # 
    # p <- p + annotation_logticks(sides = "lr")   + theme(panel.grid.minor = element_blank())
    #p <- p + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
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
                    axis.title = element_text(size = 16, angle = 00),
                    panel.grid.minor.x = element_line( size=0.5 ,  linetype = 'solid', colour = "gray93"),
                    panel.grid.minor.y = element_line( size=0.5 ,  linetype = 'solid', colour = "gray88"),
                    panel.grid.major = element_line(size =0.5,   linetype = 'solid', colour = "gray88"),
                    panel.grid = element_blank(), axis.ticks.x = element_line(size = tick.sizes)
                    
                    
                    
                    
                    
                    
                    
    )   
    
    p <- p + labs(title = paste0("Fitted analysis model '",mod,"' with 95% confidence and raw data. N = ",length(!is.na(foo$x)),"\nResidual sum of squares = ", formatz(ssr),", residual standard deviation = ",formatz(df2)," \nPredict at input of ", 
                                 formatz(Xspec) ,", the estimate of Y is ",
                                 formatz(pspec[1])," with 95%CI: (", 
                                 formatz(pspec[2]),", ",
                                 formatz(pspec[3]),")",
                                 "\nRead back at response of ", 
                                 formatz(yspec) ,", the estimate of X is ",
                                 formatz(txpre)," with 95%CI: (", 
                                 formatz(txlow),", ",
                                 formatz(txup),")",  
                                 sep=" "),
                  caption = paste0("If X specification is missing, the mean of X is used. If Y specification is missing the prediction of X used as the Y value to read back from (dashed lines).")
    )  #   +
    
    # tried this package for plots themes but got errors
    #  theme_minimal() +
    #    theme(text = element_text(family = "Cinzel", size = 16),
    #       title = element_text(family = "Cinzel", size = 16)) -> targaryen
    
    if (print.plot==1) {print(p)}
    
    return(list(ssr=ssr,r=r, foo=foo, f=f, mod=mod, rsd2=rsd2, dfs=dfs , tybar=tybar, txbar=txbar, Xspec=Xspec, tp=tp, pspec=pspec ))
    
} 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

recommendation <- read.csv('C:\\Users\\Lenovo\\Documents/recommendation.csv',stringsAsFactors = F,header=T)

#Dashboard header carrying the title of the dashboard
header <- dashboardHeader(title = "Transformations")  

#Sidebar content of the dashboard
sidebar <- dashboardSidebar(width=300,
                            
                            sidebarMenu(
                                id = "tabs",
                                br(),
                                
                                menuItem("Code", icon = icon("bar-chart-o"),
                                         menuSubItem("Shiny",  
                                                     icon = icon("send",lib='glyphicon'), 
                                                     href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/LoQ/app.R"),
                                         
                                         
                                         menuSubItem("R",  
                                                     icon = icon("send",lib='glyphicon'), 
                                                     href = "https://raw.githubusercontent.com/eamonn2014/Functional-sensitivity/master/Rcode.R")),
                                
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                
                                
                                splitLayout(
                                    
                                    tags$div(
                                        textInput(inputId="N", label='N', width = '90%' , value="100"),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='a', label='Intercept', width = '90%' , "10"),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='b', label='Slope', width = '90%' , ".1"),
                                    )
                                    
                                ),
                                
                                splitLayout(
                                    
                                    tags$div(
                                        textInput(inputId="sigma1", label='Sigma', width = '90%' , value=".4"),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='spec', label='Y spec', width = '90%' , ""),
                                    ),
                                    
                                    tags$div(
                                        textInput(inputId='Xspec', label='X spec', width = '90%' , ""),
                                    )
                                    
                                ),
                                
                                
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                menuItem("Data generating mechanism ", icon = icon("bar-chart-o"),
                                         menuSubItem( 
                                             
                                             radioButtons(
                                                 inputId = "truth",
                                                 label =  div(h5(tags$span(style="color:white",""))),
                                                 choiceNames = list(
                                                     #   HTML("<font color='white'>All scenarios</font>"), 
                                                     tags$span(style = "color:white", "Linear Y=a+bX"), 
                                                     tags$span(style = "color:white", "Exponential Y=exp(a+bX)"), 
                                                     tags$span(style = "color:white", "Reciprocal-Y Y=1/(a+bX)"),
                                                     tags$span(style = "color:white", "Reciprocal-X Y=a+b/X"),
                                                     tags$span(style = "color:white", "Double Reciprocal Y=1/(a+b/X)"),
                                                     tags$span(style = "color:white", "Logarithmic-X Y=a+b(log(X))"), 
                                                     tags$span(style = "color:white", "Multiplicative Y=aX^b"), 
                                                     tags$span(style = "color:white", "Square Root-X Y=a+b(sqrt(X))"),
                                                     tags$span(style = "color:white", "Square Root-Y Y=(a+bX)^2"),
                                                     tags$span(style = "color:white", "S-curve Y=exp(a+b/X)"),
                                                     tags$span(style = "color:white", "Square X and Y Y^2=a+X^2/b")
                                                     
                                                 ),
                                                 choiceValues = c( "model1", "model2", "model3",  "model4", "model5", "model6",
                                                                   "model7", "model8", "model9",  "model10", "model11"),
                                                 selected=c("model2")
                                             )
                                         )  # end sub
                                         
                                ),  # end main
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                
                                menuItem("Analysis transformation", icon = icon("bar-chart-o"),
                                         
                                         menuSubItem( 
                                             
                                             radioButtons(
                                                 inputId = "ana",
                                                 label =  div(h5(tags$span(style="color:white"," "))),
                                                 choiceNames = list(
                                                     HTML("<font color='white'>Best scenario</font>"), 
                                                     tags$span(style = "color:white", "Linear Y=a+bX"), 
                                                     tags$span(style = "color:white", "Exponential Y=exp(a+bX)"), 
                                                     tags$span(style = "color:white", "Reciprocal-Y Y=1/(a+bX)"),
                                                     tags$span(style = "color:white", "Reciprocal-X Y=a+b/X"),
                                                     tags$span(style = "color:white", "Double Reciprocal Y=1/(a+b/X)"),
                                                     tags$span(style = "color:white", "Logarithmic-X Y=a+b(log(X))"), 
                                                     tags$span(style = "color:white", "Multiplicative Y=aX^b"), 
                                                     tags$span(style = "color:white", "Square Root-X Y=a+b(sqrt(X))"),
                                                     tags$span(style = "color:white", "Square Root-Y Y=(a+bX)^2"),
                                                     tags$span(style = "color:white", "S-curve Y=exp(a+b/X)"),
                                                     tags$span(style = "color:white", "Square X and Y Y^2=a+X^2/b"),
                                                     tags$span(style = "color:white", "Restricted cubic spline 4 knots")
                                                 ),
                                                 choiceValues = c( "99","1", "2", "3",  "4", "5", "6",
                                                                   "7", "8", "9",  "10", "11", "12")
                                                 ,
                                                 selected=c("2")
                                             )
                                         )    
                                         
                                )
                                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            )
                            
                            
)


frow1 <- fluidRow(
    valueBoxOutput("value1")
    ,valueBoxOutput("value2")
    ,valueBoxOutput("value3")
)

frow2 <- fluidRow(
    
    box(
        title = "Revenue per Account"
        ,status = "primary"
        ,solidHeader = TRUE 
        ,collapsible = TRUE 
        ,plotOutput("plot1", height = "750px")
    )
    
    ,box(
        title = "Revenue per Product"
        ,status = "primary"
        ,solidHeader = TRUE 
        ,collapsible = TRUE 
        ,plotOutput("revenuebyRegion", height = "750px")
    ) 
    
)




# combine the two fluid rows to make the body
body <- dashboardBody(frow1, frow2)
#completing the ui part with dashboardPage


ui <- dashboardPage(title = 'This is my Page title', header, sidebar, body, skin='blue')


# create the server functions for the dashboard  
server <- function(input, output) { 
    
    
    
    
    #some data manipulation to derive the values of KPI boxes
    total.revenue <- sum(recommendation$Revenue)
    sales.account <- recommendation %>% group_by(Account) %>% summarise(value = sum(Revenue)) %>% filter(value==max(value))
    prof.prod <- recommendation %>% group_by(Product) %>% summarise(value = sum(Revenue)) %>% filter(value==max(value))
    
    
    #creating the valueBoxOutput content
    output$value1 <- renderValueBox({
        valueBox(
            formatC(sales.account$value, format="d", big.mark=',')
            ,paste('Residual sum of squares:',sales.account$Account)
            ,icon = icon("stats",lib='glyphicon')
            ,color = "purple")
        
        
    })
    
    
    
    output$value2 <- renderValueBox({
        
        valueBox(
            formatC(total.revenue, format="d", big.mark=',')
            ,'Model sigma'
            ,icon = icon("gbp",lib='glyphicon')
            ,color = "green")
        
    })
    
    
    
    output$value3 <- renderValueBox({
        
        valueBox(
            formatC(prof.prod$value, format="d", big.mark=',')
            ,paste('N:',prof.prod$Product)
            ,icon = icon("menu-hamburger",lib='glyphicon')
            ,color = "yellow")
        
    })
    
    #creating the plotOutput content
    
    output$revenuebyPrd <- renderPlot({
        ggplot(data = recommendation, 
               aes(x=Product, y=Revenue, fill=factor(Region))) + 
            geom_bar(position = "dodge", stat = "identity") + ylab("Revenue (in Euros)") + 
            xlab("Product") + theme(legend.position="bottom" 
                                    ,plot.title = element_text(size=15, face="bold")) + 
            ggtitle("Revenue by Product") + labs(fill = "Region")
    })
    
    
    output$revenuebyRegion <- renderPlot({
        ggplot(data = recommendation, 
               aes(x=Account, y=Revenue, fill=factor(Region))) + 
            geom_bar(position = "dodge", stat = "identity") + ylab("Revenue (in Euros)") + 
            xlab("Account") + theme(legend.position="bottom" 
                                    ,plot.title = element_text(size=15, face="bold")) + 
            ggtitle("Revenue by Region") + labs(fill = "Region")
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
        
        # x <-  array(runif(N, lowerV, upperV))  # no negative values
        x <-  array(runif(N, 0, 100))  # no negative values
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
            res2 <- loq(x=x, y=y, model=model, spec= spec, print.plot=1,  Xspec=Xspec)  # run best model
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
    # md2 <- reactive({
    #     
    #     spec <- as.numeric(input$spec)
    #     Xspec <- as.numeric(input$Xspec)
    #     d <- dat()  # Get the data
    #     y <- d$y
    #     x <- d$x
    #     
    #     A <- array(NA, dim=c(12,6))
    #     M <- list(NA, dim=c(24,1))
    #     M <- array(NA, dim=c(24,1))
    #     
    #     for (j in 1:12) {
    #         
    #         res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0,  Xspec=Xspec) # don't print
    #         
    #         k <- j*2
    #         m <- k-1
    #         
    #         M[k] <- res$f 
    #         M[m] <- res$mod
    #         
    #         A[j,1] <- j          # model no
    #         A[j,2] <- res$mod    # model
    #         
    #         A[j,3] <- p0f(res$dfs)    # d.f.
    #         A[j,4] <- (res$ssr)    # sum of squared residuals
    #         A[j,5] <- (res$rsd2)   # sigma
    #         A[j,6] <- (sqrt(res$ssr/res$dfs))
    #         
    #     }
    #     
    #     A <- data.frame( A[,c(1,2)],  apply(A[,c(3,4,5,6)],2, as.numeric))
    #     
    #     A <- plyr::arrange(A,A[,4])
    #     
    #     A <- data.frame( A[,c(1,2,3)],  apply(A[,c(4,5,6)],2, formatz))
    #     
    #     return(list(  ssr=A, M=M))
    #     
    # })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # X summary stats
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # output$X <- renderPrint({
    #     
    #     d <- md()$foo
    #     return(print(summary(d$x), digits=6))
    #     
    # }) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Y summary stats
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # output$Y <- renderPrint({
    #     
    #     d <- md()$foo
    #     return(print(summary(d$obsy), digits=6))
    #     
    # }) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAIN PLOT! updated with log transformation  option
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot1 <-renderPlot({     
        
        model <- md()$model
        foo <- md()$foo
        
        spec <- as.numeric(input$spec)
        Xspec <- as.numeric(input$Xspec)
        d <- dat()  # Get the  data
        y <- d$y
        x <- d$x
        
       # if(values$uno)
            loq(x= x, y= y, model=model, spec= spec, print.plot=1,  Xspec=Xspec) # print plot
        #else 
         #   if(values$dos)
           #     loq1(x= x, y= y, model=model, spec= spec, print.plot=1,  Xspec=Xspec) # print plot
        #else
            return()  
        
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}


shinyApp(ui, server)