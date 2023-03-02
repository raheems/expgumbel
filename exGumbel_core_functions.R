
# Created: May 13,2015 by Enayetur Raheem 
# Contact: Enayetur.Raheem@unco.edu
# Last modified: 

# Code for Expontiated Gumbel Distribution
# Jointly with Sanku Dey, Shaikat Mukherjee, and Enayetur Raheem

# Generating random numbers from transmuted Rayleigh distribution

library(fitdistrplus)
library(nleqslv)
library(goftest)
library(boot)

# shape = lambda
# scale = sigma

# Random number from EG distribution

regumbel <- function(n, shape = 1, scale = 1){
if (shape > 0 & scale > 0) 
{
  u <- runif(n)
    x <- -scale * log(-log(1 - (1 - u)^(1/shape)))
    x
} 
else stop ("negative")
}

regumbel(1, shape = 2, scale = 1)
# regumbel(10, shape = -1, scale = 4)

# PDF of EG distribution
degumbel <- function(x, shape = 1, scale = .5){
  (shape/scale) * (1 - exp(-exp(-x/scale)))^(shape - 1) * exp(- x/scale) * exp(-exp(-x/scale))
}

# PDF test
degumbel(2)


# CDF of EG distribution
pegumbel <- function(q, shape = 1, scale = 1){
  1 - (1 - exp(-exp(-q/scale)))^shape
}

# CDF test: q=0.3665129 gives .5
pegumbel(0.3665129)


# quantile for EG: qqegulbel()
qegumbel <- function(p, shape = 1, scale = 1){
  -scale * log(-log(1 - (1 - p)^(1/shape)))
}

# Quantile test p=.5 gives 0.3665129
qegumbel(.5)

## Hazard function of EG distribution

hegumbel <- function(x, shape = 1, scale = 1){
  (shape * exp(-x/scale) * exp(-exp(-x/scale))) / (scale * (1 - exp(-exp(-x/scale))))
}

# 
# 
# ############### EG with location ####################
# 
# PDF of EG distribution
deGumbel <- function(x, shape = 1, scale = 2, location = 5){
  eterm <- exp(-(x - location)/scale)
  (shape/scale) * (1 - exp(-eterm))^(shape - 1) * eterm * exp(-eterm)
}

# PDF test
deGumbel(2)



# CDF of EG distribution
peGumbel <- function(q, shape = 1, scale = 2, location = 5){
  eterm <- exp(-(q - location)/scale)
  1 - (1 - exp(-eterm))^shape
}

# CDF test: q=0.3665129 gives .5
peGumbel(5.733026)


# quantile for EG: qqegulbel()
qeGumbel <- function(p, shape = 1, scale = 2, location = 5){
  location + -scale * log(-log(1 - (1 - p)^(1/shape)))
}

# Quantile test p=.5 gives 0.3665129
qeGumbel(.5)



# quantile for EG: qqegulbel()
reGumbel <- function(n, shape = 1, scale = 2, location = 5){
  u <- runif(n)
  location + -scale * log(-log(1 - (1 - u)^(1/shape)))
}

reGumbel(10, location = 0)

# fn(x=dat1, n = 1, shape=1, scale=1)
# 
# dat1 <- regumbel(100)
# integrate(fn, -Inf, 100, n=2, shape=5, scale=4) 

############################################
## GUMBEL DISTRIBUTION
############################################
# 
# The density, cdf, and random number generator for 
# Type II Gumbel distribution 

# Gumbel PDF
dGumbel <- function(x, alpha, beta) {
  alpha * beta * (x^(-alpha - 1)) * exp(-beta * x^(-alpha))
}

# 
# Gumbel  CDF (distribution function)
# alpha, beta real numbers

pGumbel <- function(q, alpha, beta){
  exp(-beta * q^(-alpha))
}
 
# q=1.442695 gives 0.5 
pGumbel(1.442695, alpha=1, beta=1)
# 

### Gumbel Quantile function###
qGumbel <- function(p, alpha, beta) {
  (beta/(log(1/p)))^(1/alpha)
}

# p=.5 should give a value close to 1.83
qGumbel(.5, alpha=1, beta=1)


# Gumbel Random number generator: quantile function
rGumbel <- function(n, alpha, beta) {
  u <- runif(n) 
  (beta/(log(1/u)))^(1/alpha)
}

rGumbel(1, alpha = 1, beta = 1)


############################################
## Frechet DISTRIBUTION
############################################
# 
# The density, cdf, and random number generator for 
# Two-parameter Frechet distribution 

# PDF
dfrechet <- function(x, alpha, beta) {
  (alpha/beta) * (beta / x)^(alpha + 1) * exp(-(beta/x)^alpha)
}

# 
# CDF (distribution function)
# alpha, beta real numbers

pfrechet <- function(q, alpha, beta){
  exp(-(beta/q)^alpha)
}

# q=1.201122 gives 0.5 
pfrechet(2.402245, alpha=2, beta=2)
# 

# QF
qfrechet <- function(p, alpha, beta) {
  beta / ((- log(p))^(1 / alpha))
}

# p=.5 should give a value close to 1.83
qfrechet(.5, alpha=2, beta=2)


# Random number generator: quantile function
rfrechet <- function(n, alpha, beta) {
  u <- runif(n) 
  beta / ((- log(u))^(1 / alpha))
}

rfrechet(1, alpha = 1, beta = 1)


############################################
## Generalized ExtremeValue DISTRIBUTION
############################################
# 
# The density, cdf, and random number generator for 
# Three-parameter GEV
# alpha, sigma, mu (mu is the location parameter)
# consider mu = 0

# PDF
dgev <- function(x, shape, scale, location) {
  part.shape <- shape * ((x - location)/scale)
  part1 <- 1/scale
  part2 <- exp(-((1 + part.shape)^(-1/shape)))
  part3 <- (1 + part.shape)^(-(1 + 1/shape))
  return(part1 * part2 * part3)
}

dgev(wdat, shape=.22, scale=2, location=2.5)

# 
# CDF (distribution function)

pgev <- function(q, shape, scale, location){
  exp(-(1 + shape * ((q - location)/scale))^(-1/shape))
}

# q=5.442695 gives 0.5 
pgev(5.442695, shape=1, scale=1, location=5)


# QF
qgev <- function(p, shape, scale, location) {
  location + (scale/shape) * ((-log(p))^(-shape) - 1)
}

# p=.5 should give a value close to 5.44
qgev(.5, shape=1, scale = 1, location=5)


# Random number generator: quantile function
rgev <- function(n, shape, scale, location) {
  u <- runif(n) 
  location + (scale/shape) * ((-log(u))^(-shape) - 1)
}

rgev(1, shape = 1, scale = 1, location = 5)



########################################
###### Maximum Likelihood Estimation ###
########################################

# log-liklihood function for ExGumbel distribution

fun.ml <- function(param, x){
  # ML method for Ex Gumbel distribution
  n<-length(x)
  shape <- param[1]
  scale <- param[2]
  location <- param[3]
  n * log(shape) - n * log(scale) - (1/scale) * sum((x - location)) + (shape - 1) * sum(log(1 - exp(- exp(-(x - location)/scale)))) - sum(exp(-(x - location)/scale))
}

x <- reGumbel(100, shape=100, scale=200, location=5)
fun.ml(c(1, 2, 4), x)

# Function for 2-parameter ExGumbel distribution
# for simulation

fun.ml2 <- function(param, x){
  # ML method for Ex Gumbel distribution
  n<-length(x)
  shape <- param[1]
  scale <- param[2]
  n * log(shape) - n * log(scale) - (1/scale) * sum((x)) + (shape - 1) * sum(log(1 - exp(- exp(-(x)/scale)))) - sum(exp(-(x)/scale))
}

x <- regumbel(100, shape=1, scale=2)
fun.ml2(c(1, 2), x)


####################################
## Least Squares Estimation (LSE)
####################################

fun.lse <- function(param, x){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  val <- NULL
  for (i in 1:n) {
    val[i] <- (1 - (1 - exp(-exp(-x[i]/scale)))^shape - i/(n + 1))^2
  }
  sum(val)
}


####################################
## Weighted Least Squares Estimation (WLSE)
####################################

fun.wlse <- function(param, x){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  val <- NULL
  for (i in 1:n) {
    val[i] <- ((n + 1)^2 * (n + 2))/(n - i + 1) * (1 - (1 - exp(-exp(-x[i]/scale)))^shape - i/(n + 1))^2
  }
  sum(val)
}


####################################
## Percentile Estimators (PCE)
####################################

# This function has been modified 
# for Ex Gumbel Model

fun.pce <- function(param, x){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  val <- NULL
  for (i in 1:n) {
    #p.i <- pegumbel(x[i], shape = shape, scale = scale)
    p.i <- i/(n+1)
    inner.val <- 1 - (1 - p.i)^(1/shape)
    inner.val <- ifelse(inner.val == 1, 0.999, inner.val)
    val[i] <- (x[i] + (scale * log(-log(inner.val))))^2
  }
  sum(val)
}

      

####################################
## Method of Maximum Product Spacing (MPS)
####################################


fun.mps <- function(param, x){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(break.ties(x))
  val <- NULL
  for (i in 1:n) {
    val[i] <- pegumbel(x[i], shape = shape, scale = scale)
  }
  d <- c(0, val, 1)
  d <- d[-1] - d[- (n + 2)]
  d <- d[d!=0] # remove the zeros
  sum(log(d))/ (length(d) + 1)
}
# 
# dat1 <- regumbel(100, shape = 4, scale = 3)
# optim(par = c(3, 2), fn=fun.mps, 
#              #method = "L-BFGS-B",
#              #lower = c(0.01, 0.01),
#              #upper = c(Inf, Inf),
#              control=list(fnscale=-1), 
#              x=dat1
# )



####################################
## Cramer von Mises (CVM)
####################################

fun.cvm <- function(param, x){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  val <- NULL
  for (i in 1:n) {
    val[i] <- (pegumbel(x[i], shape = shape, scale = scale) - (2 * i - 1)/ (2 * n))^2
  }
  1/(12 * n) + sum(val)
}



####################################
## Anderson Darling Estimator (AD)
####################################

fun.ad <- function(param, x){
  x <- break.ties(x)  # Remove ties first
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  val <- NULL
  for (i in 1:n) {
    val[i] <- (2 * i - 1) * (log(pegumbel(x[i], shape = shape, scale = scale)) + log(1 - pegumbel(x[n + 1 - i],  shape = shape, scale = scale)))
  }
  n <- length(na.omit(val)) # remove the NaNs
  -n - sum(na.omit(val))/n  # based on valid "val"s
}



####################################
## Right Tailed Anderson Darling Estimator (RTAD) 
####################################

fun.rtad <- function(param, x){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  val1 <- NULL
  val2 <- NULL
  for (i in 1:n) {
    val1[i] <- pegumbel(x[i], shape = shape, scale = scale)
    val2[i] <- (2*i - 1) * log(1 - pegumbel(x[n + 1 - i],  shape = shape, scale = scale))
  }
  (n/2) - (2 * sum(val1)) - (sum(val2)/n)
}




####################################
## Method of Moment Estimator (MME) 
####################################


## nth Moment
nth.moment.egumbel.fn <- function(x, n, shape, scale) {
 val1 <- x^n 
  val1 <- ifelse((is.nan(val1) | is.infinite(val1)), 0, val1)
 val2 <- (shape - 1) * log((1 - exp(- exp(-x/scale))))
 val2 <- ifelse((is.nan(val2) | is.infinite(val2)), 0, val2)
 val3 <- exp(-x/scale)
 val3 <- ifelse((is.nan(val3) | is.infinite(val3)), 0, val3)
 val4 <- exp(-exp(-x/scale))
 val4 <- ifelse((is.nan(val4) | is.infinite(val4)), 0, val4)
 #return(ifelse((x == 0 | x == Inf), 0, val1 * val2 * val3 * val4))
 return(val1 * val2 * val3 * val4)
}

# Workout with 
shape = 0.05; scale = 0.65; n = 1
x <- -5:100

val1 <- x^n 
val1
val1 <- ifelse((is.nan(val1) | is.infinite(val1)), 0, val1)
val1

val2 <- (shape - 1) * log((1 - exp(- exp(-x/scale))))
val2
val2 <- ifelse((is.nan(val2) | is.infinite(val2)), 0, val2)
val2

val3 <- exp(-x/scale)
val3
val3 <- ifelse((is.nan(val3) | is.infinite(val3)), 0, val3)
val3

val4 <- exp(-exp(-x/scale))
val4
val4 <- ifelse((is.nan(val4) | is.infinite(val4)), 0, val4)
val4

rm(shape, scale, n, x, val1, val2, val3, val4)

# See what values of x produce valid values for the the function

nth.moment.egumbel.fn(-5:100, n=1, shape=0.05, scale = 0.65)
nth.moment.egumbel.fn(-5:100, n=2, shape=0.05, scale = 0.65)

# nth.moment.egumbel.fn(dat1, 1, 4, 2)

nth.moment.egumbel <- function(n, low, up, shape, scale) {
  int <-integrate(nth.moment.egumbel.fn, low=low, up=up, n = n, shape = shape, scale = scale)
  (shape/scale) * int$val
}

# First moment about origin
nth.moment.egumbel(n = 1, low=0, up=20, shape = .05, scale = 0.65)

# Second moment about origin
nth.moment.egumbel(n = 2, low=0, up=20, shape = 0.05, scale = 0.65)


fun.mme <- function(param, x, low, up){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  euler <- 0.577215
  x <- sort(x)
  l1 <- sum(x)/n 
  l2 <- sum(x^2)/n
  L1 <- nth.moment.egumbel(n = 1, shape = shape, scale = scale, low=low, up=up)
  L2 <- nth.moment.egumbel(n = 2, shape = shape, scale = scale, low=low, up=up)
  
  # Writing the main function  
  est <- numeric(2)
  est[1] <- L1 - l1
  est[2] <- L2 - l2
  est
}

# # example
# dat1 <- regumbel(1000, shape=10, scale = 2)
# fun.mme(c(3, 2), dat1, low=-5, up=Inf)
# 
# # Using BB package
# BBsolve(c(9, 1.5), fun.mme, x = dat1, low=-5, up=Inf)

####################################
## Method of Modified Moment Estimator (MMM) 
####################################
# 
# fun.mmm <- function(param, x){
#   n <- length(x)
#   alpha <- param[1]
#   lambda <- param[2]
#   x <- sort(x)
#   l1 <- mean(x) 
#   l2 <- var(x)
#   L1 <- (1/sqrt(alpha)) * gamma(1 + 1/2) * (1 - lambda + (lambda/sqrt(2)))
#   L2 <- (1/alpha) * (1 - lambda + (lambda/2)) - L1^2
#   
#   # Writing the main function  
#   est <- numeric(2)
#   est[1] <- L1 - l1
#   est[2] <- L2 - l2
#   est
# }



####################################
## L-Moment Estimator (LME) 
####################################

# Computation of E(X22)
int21 <-function(x, shape, scale){
  int.res <- log(exp(-x/scale)) * exp(-exp(-x/scale)) * (1 - exp(-exp(-x/scale)))^(2*shape - 1)
  int.res[is.na(int.res)] <- 0
  int.res
} 

#Testing
# int21(dat1, 3, 2)
integrate(int21, 0, Inf, shape=3, scale=2)$val

int22 <- function(x, shape, scale){
  int.res <- log(exp(-x/scale)) * exp(-exp(-x/scale)) * (1 - exp(-exp(-x/scale)))^(shape - 1)
  int.res[is.na(int.res)] <- 0
  int.res
}
# Testing
# int22(dat1, 3, 2)
integrate(int22, 0, Inf, shape=3, scale=2)$val

fun.lme <- function(param, x, low, up){
  n <- length(x)
  shape <- param[1]
  scale <- param[2]
  x <- sort(x)
  l1 <- sum(x)/n 
  val <- NULL
  for (i in 1:n) {
    val[i] <- (i - 1) * x[i]
  }
  l2 <- (2/(n * (n - 1))) *  sum(val) - l1
  
  # L1 is integration of nth moment for n=1
  L1 <- nth.moment.egumbel(n = 1, low=low, up=up, shape = shape, scale = scale)
  int21.obj <- integrate(int21, low = low, up = up, shape=shape, scale=scale)
  int22.obj <- integrate(int22, low = low, up = up, shape=shape, scale=scale)
  L2 <- 2 * shape * scale * int21.obj$val - shape * scale * int22.obj$val
  
  # Writing the main function  
  est <- numeric(2)
  est[1] <- L1 - l1
  est[2] <- L2 - l2
  est
}

# Testing the function fun.lme()
dat1 <- regumbel(100, shape = 3, scale = 2)
fun.lme(c(2.5, 1.5), dat1, low=-5, up=20)



###### Issues with tied observations ###
#### Perturb data to break ties ########

break.ties <- function(x){
  while(any(duplicated(x)==TRUE)){
    n <- length(x)
    unq <- x[duplicated(x)==FALSE]
    dup <-x[duplicated(x)==TRUE]
    r1 <- runif(length(dup), min=.001, max=.01)
    x <- c(unq, dup+r1)
  }
 x
}

# 
# x2 <- c(1, 2, 2, 4, 4, 5)
# break.ties(x2)
# 
# x3 <- rep(2, 10)
# x3
# break.ties(x3)
