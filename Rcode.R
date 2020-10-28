#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  rm(list=ls()) 
  set.seed(333) # reproducible
  
  library(rms)
  library(ggplot2)
  library(tidyverse)
  library(tvthemes)
  
  options(max.print=1000000)    
  
  ## convenience functions
  pf0 <- function(x) {formatC(x, format="f", digits=0)}
  p1f <- function(x) {formatC(x, format="f", digits=1)}
  p2f <- function(x) {formatC(x, format="f", digits=2)}
  p3f <- function(x) {formatC(x, format="f", digits=3)}
  p4f <- function(x) {formatC(x, format="f", digits=4)}
  p5f <- function(x) {formatC(x, format="f", digits=5)}
  p4f <- function(x) {formatC(x, format="f", digits=4)}
  
  logit <- function(p) log(1/(1/p-1))
  expit <- function(x) 1/(1/exp(x) + 1)
  inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
  is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
  
  options(width=200)
  options(scipen=999)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create data # range of independent variable
  lowerV=0
  upperV=10
  N=100
  a=10 
  b=.1 
  sigma=.1
  spec=15
  Xspec=5
  model="model11"    

  x <-  array(runif(N, lowerV, upperV))  # no negative values            
  noise <-  rnorm(N,0, sigma)           # residual error

  # create response according to data generation mechanism
  
  if (model %in% "model1") {
    y <-  a+ x*b +    noise
  } else if (model %in% "model2") {
    y <-  exp(a+ x*b + noise)
  } else if (model %in% "model3") {
    y <-  1/(a+ x*b +  noise)
  } else if (model %in% "model4") {
    y <-  a + b*(1/x) + noise
  } else if (model %in% "model5") {
    y <-  1/(a + b/x +  noise)
  } else if (model %in% "model6") {
    y <-  a + log(x)*b +    noise   
  } else if (model %in% "model7") {
    y <-  a * x^b +    noise
  } else if (model %in% "model8") {
    y <-  a + sqrt(x)*b +    noise
  } else if (model %in% "model9") {
    y <-  (a + x*b + noise)^2 
  } else if (model %in% "model10") {
    y <-  exp(a+ b/x + noise)
  } else if (model %in% "model11") {
    y <-  ( a + (x^2)/b + noise)^.5
  }
  
  d <- as.data.frame(cbind(x,y))
  y <- d$y
  x <- d$x
  
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # function that does all the work!
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  loq <- function (x, y, model, spec, print.plot=1, Xspec) {
    
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
    
    # transform spec for prediction, not only where x is trandsformed, above
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
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RCS MODEL
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
      
      if (is.na(spec))  {spec=mean(dat2$linear.predictors, is.finite=TRUE)}
      if (is.na(Xspec)) {Xspec=mean(dat$x, is.finite=TRUE)}
      
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
      
      # rcs will report the limit as the nearest value to spec if x value is beyond range, so lets report 999
      txpre <- ifelse((txpre >= upperV )|(txpre <= lowerV ), 999, txpre)
      txlow <- ifelse((txlow >= upperV )|(txlow <= lowerV ), 999, txlow)
      txup <-  ifelse((txup >=  upperV )|(txup  <= lowerV ), 999, txup)
      
      rsd2 <-  anova(f)["ERROR","MS"]^.5
      dfs <- anova(f)["ERROR","d.f."]
      
      XXX  <- predict(f, Xspec, se.fit=TRUE) 
      XXX$L <- XXX$linear.predictors - qt(0.975,n-4) * XXX$se.fit     # n-4 as we are using rcs 4 df are used up
      XXX$U <- XXX$linear.predictors + qt(0.975,n-4) * XXX$se.fit
      pspec <- as.vector(unlist(XXX))
      
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
      
      if (is.na(Xspec)) {Xspec=mean(x, is.finite=TRUE)}
      
      pspec <- predict.lm(f, newdata=data.frame(x=Xspec), interval="confidence")
      
      if (model %in% c(2,7,10)) {pspec <- exp(pspec)}
      if (model %in% c(3,5)  )  {pspec <- 1/pspec}
      if (model %in% c(9)    )  {pspec <- (pspec)^2}
      if (model %in% c(11)   )  {pspec <- pspec^.5}
      
      if( pspec[3] < pspec[2] ) {pspec <- pspec[c(1,3,2)] }
      
      # residuals original y and transformed back predicted values, residual sum of squares, this will be used to judge best model
      r <- (y1-p[,1]) 
      r2 <-  (y1-p[,1])^2
      ssr <- sum(r2, na.rm=T) 
      
      # transform the specification that we will read back from, note only those in which y is transformed
      if (is.na(spec)) {spec=mean(y, is.finite=TRUE)}    
      else { 
        if (model %in% c(2,7,10)) {spec <- log(spec)}
        if (model %in% c(3,5)  )  {spec <- 1/spec}
        if (model %in% c(9)    )  {spec <- sqrt(spec)}
        if (model %in% c(11)   )  {spec <- spec ^2}
      }
      
      tyspec <- spec     
      if (model %in% c(2,7,10)) {spec <- exp(spec)} 
      if (model %in% c(3,5)  )  {spec <- 1/spec} 
      if (model %in% c(9)    )  {spec <- (spec)^2} 
      if (model %in% c(11)   )  {spec <- spec^.5}  
      
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
    ymin <- min(y)
    ymax <- max(y)
    ystep <- (ymax-ymin)/8
    ymin1 <-  ymin-ystep
    ymax1 <-  ymax+ystep
    
    # plot and present the estimated read back
    p1 <- ggplot(foo, aes(x=x,y=pred)) + 
      geom_line( ) +
      geom_ribbon(data=foo , aes(ymin= p2a,ymax= p3),alpha=0.2,   fill="green") +
      geom_point(data=foo, aes(x=x ,y=obsy), size=2, color='blue')  +
      scale_x_continuous(limits = c(lowerV, upperV), breaks=lowerV: upperV)  
    scale_y_continuous(limits = c(ymin1, ymax1))   
    
    p <- p1  + geom_hline(yintercept=spec,  colour="#990000", linetype="dashed")
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
    
    p <- p + labs(title = paste0("Fitted model '",mod,"' with 95% confidence and raw data. Residual sum of squares = ", p2f(ssr),", residual standard deviation = ",p2f(df2)," \nPredict at input of ", 
                                 p2f(Xspec) ,", the estimate of Y is ",
                                 p4f(pspec[1])," with 95%CI: (", 
                                 p4f(pspec[2]),", ",
                                 p4f(pspec[3]),")",
                                 "\nRead back at response of ", 
                                 p2f(spec) ,", the estimate of X is ",
                                 p2f(txpre)," with 95%CI: (", 
                                 p2f(txlow),", ",
                                 p2f(txup),")",  
                                 sep=" "),
                  caption = paste0("If spec is missing, mean of the data is used")
    )  #   +
    
    # tried this package for plots themes but got errors
    #  theme_minimal() +
    #    theme(text = element_text(family = "Cinzel", size = 16),
    #       title = element_text(family = "Cinzel", size = 16)) -> targaryen
    
    if (print.plot==1) {print(p)}
    
    return(list(ssr=ssr,r=r, foo=foo, f=f, mod=mod, rsd2=rsd2, dfs=dfs  ))
    
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ssr <- rep(NA,12)
  
  mdata <- list(NA)
  model=2  # data generating
  spec=mean(y)
  Xspec=mean(x)
  
  res <- loq(x=x, y=y, model=2, spec=NA  , print.plot=1, Xspec=NA ) # don't
  
  
  # best model
    for (j in 1:12) {
      
      res <- loq(x=x, y=y, model=j, spec= NA, print.plot=0, Xspec=NA) # don't print
      ssr[j] <- res$ssr
      
    }
    
    model <- which(ssr==min(ssr)) 
    mdata <- res$foo
    res2 <- loq(x=x, y=y, model=model, spec=  NA, print.plot=1,  Xspec=NA)  # run best model
     #
    
    
    ## select any model
    res <- loq(x=x, y=y, model=12, spec= NA,  Xspec=NA, print.plot=1) 
    head(res$foo)
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
  
  for (j in 1:12) {    ###############NEW
    
    res <- loq(x=x, y=y, model=j, spec= spec, print.plot=0) # don't print
    ssr[j] <- res$ssr
    
  }
  
  model <- which(ssr==min(ssr))     # capture which model has lowest residual sum of squares
  mdata <- res$foo                  # capture data
  res2 <- loq(x=x, y=y, model=model, spec= spec, print.plot=0)  # run best model
  f=res2$f                          # linear model captured here
  mod<- res2$mod                    # linear model text description
  
  
 ###### of choosing
  
  res <- loq(x=x, y=y, model=2, spec= spec) 
  mdata <- res$foo
  model <- as.numeric(ana)
  f=res$f
  mod<- res$mod
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run model again based on above

res <-  loq(x= x, y= y, model=model, spec= spec, print.plot=1) # print plot

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# capture data

names(res$foo) <- c("x","y","prediction","lower 95%CI", "upper 95%CI", "residual","residual^2")
foo <- plyr::arrange(res$foo,x)
foo

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# capture model, object f and investigate residuals

f <- res$f

resid <- r <- resid(f)
fitted <- fitted(f)
d <- cbind(resid, fitted)
d2 <- as.data.frame(d)

yl <- ylab('Residuals')

xl <- xlab("time")

p1 <- ggplot(d2 , aes(x=fitted , y=resid)) + geom_point (   colour="#69b3a2") + yl

p2 <- ggplot(d2 , aes(sample=resid )) + stat_qq(colour="#69b3a2") +
  geom_abline(intercept=mean(r), slope=sd(r)  ,  colour="black") +
  xlab('Residuals')   +
  ggtitle( " ")

library(gridExtra)
library(grid)
df <- data.frame(Residuals = r)
p3 <- ggplot(df, aes(x = Residuals)) +
  geom_histogram(aes(y =..density..),
                 #breaks = seq(-50, 50, by = 2),
                 colour = "black",
                 fill = "#69b3a2") +
  xlab('Residuals with superimposed sigma')   

###############NEW
if (! model %in% 12) {          
  p3 <- p3 + stat_function(fun = dnorm, args = list(mean = 0, sd = as.numeric(sigma)   )) + 
    stat_function(fun = dnorm, args = list(mean = 0, sd = sigma(f)    ), col='red')  
  std <-  sigma(f) 
} else {
  
  p3 <-  p3 + stat_function(fun = dnorm, args = list(mean = 0, sd = as.numeric(sigma)   )) + 
    stat_function(fun = dnorm, args = list(mean = 0, sd =  f$stats["Sigma"][[1]]    ), col='red') 
  std <-  f$stats["Sigma"][[1]]
}

grid.arrange(p1,  p3, p2, ncol=2,
             top = textGrob(paste0(" LS model fit diagnostics, ",mod,", true sigma (black) ",as.numeric(sigma)  ,", estimated sigma (red) ", p4(std),""),gp=gpar(fontsize=20,font=3)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
