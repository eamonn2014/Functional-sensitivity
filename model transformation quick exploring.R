


set.seed(333)
N <-1000
a<- 10
b <- 1
sigma <- 2

x <-  array(runif(N, 0, 10))
noise <-  rnorm(N,0, sigma)

y <-  a+ x*b +    noise
lm(y ~ x)

y <-  exp(a+ x*b + noise)
lm(log(y) ~ x)


y <-  1/(a+ x*b +  noise)
lm(1/y ~ x)


#y <-   a + b/(x+  noise)
y <-   a + b*(1/x) + noise
X <- 1/x
lm(y ~ X)   


y <-   1/(a + b/x +  noise)
Y <- 1/y ; X <-1/x
lm(Y ~ X)  


y <-  a + log(x)*b +    noise   
lm(y ~ log(x))


y <-  a * x^b +    noise
Y ~ log(y) ; X ~ log(x)
lm(Y ~ X)


y <-  a + sqrt(x)*b +    noise
lm(y ~ sqrt(x))

 
y <-  (a + x*b + noise)^2 
lm(sqrt(y) ~ (x)) 


y <-  exp(a+ b/x + noise)
Y ~ log(y) ; X ~ 1/(x)
lm(Y ~ X)


y <-  ( a + (x^2)/b + noise)^.5
Y ~ (y)^2 ; X ~x^2
lm(Y ~ X)
