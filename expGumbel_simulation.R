
# Exponentiated Gumbel Distribution
# Simulation
# Date: May 13, 2015
# Date modified: Dec 23, 2015


library(boot)
library(fitdistrplus)
library(BB)
library(nleqslv)

# Generate 70 observations from EGumbel 
# distribution with shape = 5 and scale = 2

set.seed(2)
d1 <- regumbel(100, shape = 5, scale = 2)
summary(d1)

# Plot the histogram and empirical cdf
plotdist(d1, histo = TRUE, demp = TRUE)

# Descriptive summary of distribution
descdist(d1, boot = 1000)


feg <- fitdist(d1, distr ="egumbel", start=list(shape = 4, scale = 1), optim.method = "L-BFGS-B", lower = c(.001, .001), upper = c(Inf, Inf))
summary(feg)

# Fitting additional distributions
fn <- fitdist(d1, "norm")

# Plotting these distributions
plot.legend <- c("EG", "Normal")
denscomp(list(feg, fn), legendtext = plot.legend)
qqcomp(list(feg, fn), legendtext = plot.legend)
cdfcomp(list(feg, fn), legendtext = plot.legend)
ppcomp(list(feg, fn), legendtext = plot.legend)


########################################
####### FUNCTIONS FOR SIMULATION #######
########################################
# True parameters: 
# lambda = .5, 1, 3, sigma = .5
# lambda = .5, 1, 3, sigma = 3
# n = 20, 50, 100, 200

#  MLE DONE#

# Simulation for 2-parameter ExGumbe distribution
# Notice hte fn=fun.ml2 option

run.mle <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n*R, shape = shape, scale = scale), ncol=R)

fun.sim <- function(data){
# Optimizing the loglik to estimate the alpha and lambda
  out <- optim(par = c(start[1], start[2]), fn=fun.ml2, x=data,
               #method = "L-BFGS-B",
               #lower = c(0.01, 0.01),
               #upper = c(Inf, Inf),
               control=list(fnscale=-1)
  )
out$par  
}

out.est <- apply(data.sim, 2, fun.sim)
out.est <- as.data.frame(t(out.est))
colnames(out.est) <- c("shape", "scale")

# Biases
bias.shape <- sum(out.est$shape - shape)/R
bias.scale <- sum(out.est$scale - scale)/R

# Root MSEs
rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)

F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale = scale))

F.emp <- matrix(0, nrow=n, ncol=R)

for(j in 1: R) {
  for (i in 1:n) {
    par <- as.matrix(out.est[j, ])
    F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale = par[2])
  }
}

Dobs <- sum(abs(F.true - F.emp))/(n * R)
Dmax <- mean(apply(abs(F.true - F.emp), 2, max))

list(true.par = c(shape, scale),
     est = c(mean(out.est$shape), mean(out.est$scale)),
     bias.shape = bias.shape,
     rmse.shape = rmse.shape,
     bias.scale = bias.scale,
     rmse.scale = rmse.scale,
     Dobs = Dobs,
     Dmax = Dmax, 
     n=n, R=R)
}

run.mle(n = 20, R= 100, param=c(2, 3), start=c(1.5, 2.5))


## LSE ##

run.lse <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n*R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.lse, 
                 method = "L-BFGS-B",
                 lower = c(0.01, 0.01),
                 upper = c(Inf, Inf),
                 #control=list(fnscale=-1), 
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

run.lse(n = 100, R= 10, param=c(3, 2), start=c(2, 1))



## WLSE ## 

run.wlse <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.wlse, 
                 method = "L-BFGS-B",
                 lower = c(0.01, 0.01),
                 upper = c(Inf, Inf),
                 #control=list(fnscale=-1), 
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}
ini = c(1.25, .4)
run.wlse(n = 100, R= 10, param=c(3, 2), start=c(2, 1))


## PCE ## Done but doubtful results

run.pce <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.pce, 
                 method = "L-BFGS-B",
                 lower = c(0.01, 0.01),
                 upper = c(Inf, Inf),
                 #control=list(fnscale=-1), 
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

run.pce(n = 100, R= 100, param=c(5, 3), start=c(4, 2.5))


## MPS ## Done
# Maximization problem: fnscale=-1

run.mps <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.mps, 
                 method = "L-BFGS-B",
                 lower = c(0.01, 0.01),
                 upper = c(5, 5),
                 control=list(fnscale=-1), 
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

run.mps(n = 20, R = 100, param=c(4, 2), start=c(4, 3))



## CVM ##

run.cvm <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.cvm, 
                 method = "L-BFGS-B",
                 lower = c(0.1, 0.1),
                 upper = c(Inf, Inf),
                 #control=list(fnscale=-1), 
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

ini <- c(3, 3)
run.cvm(n = 100, R = 100, param=c(4, 2), start=ini)


## AD ##
# Run with Nelder-Mead

run.ad <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.ad, 
                 #method = "L-BFGS-B",
                 #lower = c(0.1, 0.1),
                 #upper = c(10, 10),
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

ini <- c(3, 1.5)
run.ad(n = 30, R = 100, param=c(5, 2), start=ini)


## RTAD ##

run.rtad <- function(n, R = 10, param=c(1, 1), start) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
    out <- optim(par = c(start[1], start[2]), fn=fun.rtad, 
                 #method = "L-BFGS-B",
                 #lower = c(0.01, 0.01),
                 #upper = c(10, 10),
                 #control=list(fnscale=-1), 
                 x=data
    )
    out$par  
  }
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}
ini <- c(4, 2)
run.rtad(n = 40, R = 100, param=c(5, 3), start = ini)



#  MME # 
# Results are questionnable 
# library(nleqslv)
# library(BB)

run.mme <- function(n, R = 10, param=c(1, 1), start, low, up) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n*R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
#     # Using BB package
#     out <- BBsolve(c(start[1], start[2]), 
#                    fun.mme, 
#                    x = data, low = -5, up = Inf
#                    )
#     out$par
    out<- nleqslv(c(start[1], start[2]), fun.mme, jac=NULL, low=low, up= up, data, control=list(btol=.01, delta="cauchy"))
    out$x
  }
  
  # fun.sim(data.sim[,1])
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

# For (1, .7), use start=c(1.25, .4)
ini = c(5, 2)
run.mme(n = 20, R= 100, param=c(5, 2), start = ini, low=1, up=100)


# LME (Using nleqslv package)
# library(nleqslv)

run.lme <- function(n, R = 10, param=c(1, 1), start, low, up) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
#     # Using BB package
#     out <- BBsolve(c(start[1], start[2]), 
#                    fun.lme, low=low, up=up,
#                    x = data
#     )
#      out$par
    out<- nleqslv(c(start[1], start[2]), fun.lme, jac=NULL, low=low, up=up, data, control=list(btol=.01,delta="cauchy", maxit = 200, allowSingular=T))
    out$x
  }
  
    
  # fun.sim(data.sim[,1])
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

# For (1, .7), use start=c(1.25, .4)
ini = c(2, 1)
run.lme(n = 100, R= 10, param=c(2, 1), start=ini, low = -2, up = 10)

# LME (Using BB package)
# library(BB)

run.lme.bb <- function(n, R = 10, param=c(1, 1), start, low, up) {
  shape <- param[1]
  scale <- param[2]
  set.seed(2)
  data.sim <- matrix(regumbel(n * R, shape = shape, scale = scale), ncol=R)
  
  fun.sim <- function(data){
        # Using BB package
        out <- BBsolve(c(start[1], start[2]), 
                       fun.lme, low=low, up=up,
                       x = data
        )
         out$par
  }
  
  
  # fun.sim(data.sim[,1])
  out.est <- apply(data.sim, 2, fun.sim)
  out.est <- as.data.frame(t(out.est))
  colnames(out.est) <- c("shape", "scale")
  
  # Biases
  bias.shape <- sum(out.est$shape - shape)/R
  bias.scale <- sum(out.est$scale - scale)/R
  
  # Root MSEs
  rmse.shape  <- sqrt(sum((out.est$shape - shape)^2)/R)
  rmse.scale  <- sqrt(sum((out.est$scale - scale)^2)/R)
  
  F.true <- apply(data.sim, 2, function(x) pegumbel(x, shape=shape, scale=scale))
  
  F.emp <- matrix(0, nrow=n, ncol=R)
  
  for(j in 1: R) {
    for (i in 1:n) {
      par <- as.matrix(out.est[j, ])
      F.emp[, j] <- pegumbel(data.sim[,j], shape=par[1], scale=par[2])
    }
  }
  
  Dobs <- sum(abs(F.true - F.emp))/(n * R)
  Dmax <- mean(apply(abs(F.true - F.emp), 2, max))
  
  list(true.par = c(shape, scale),
       est = c(mean(out.est$shape), mean(out.est$scale)),
       bias.shape = bias.shape,
       rmse.shape = rmse.shape,
       bias.scale = bias.scale,
       rmse.scale = rmse.scale,
       Dobs = Dobs,
       Dmax = Dmax, 
       n=n, R=R)
}

# For (1, .7), use start=c(1.25, .4)
ini = c(2, 1)
run.lme.bb(n = 100, R= 10, param=c(2, 1), start=ini, low = 1, up = 10)



#################################################
############ RUNNING SIMULATION #################
############ lambda = .5, sigma =.5 ##############
#################################################

# MLE (done for EG)
# run.mle(n = 20, R= 100, param=c(.5,.5), start=c(.4, .4))

ini = c(.4, .4)
mle.n20l.5s.5 <- run.mle(n = 20, R=5000, param=c(.5,.5), start=ini)
mle.n50l.5s.5 <- run.mle(n = 50, R= 5000, param=c(.5,.5), start=ini)
mle.n100l.5s.5 <- run.mle(n = 100, R= 5000, param=c(.5,.5), start=ini)
mle.n200l.5s.5 <- run.mle(n = 200, R= 5000, param=c(.5,.5), start=ini)


# MME (done for EG)
mme.n20l.5s.5 <- run.mme(n=20, R=5000, param=c(.5,.5), start=ini, low=1, up=100)
mme.n50l.5s.5 <- run.mme(n=50, R=5000, param=c(.5,.5), start=ini, low=1, up=100)
mme.n100l.5s.5 <- run.mme(n=100, R=5000, param=c(.5,.5), start=ini, low=1, up=100)
mme.n200l.5s.5 <- run.mme(n=200, R=5000, param=c(.5,.5), start=ini, low=1, up=100)

# LME (DOES NOT WORK)
lme.n20l.5s.5 <- run.lme(n=20, R=5000, param=c(.5,.5), start=ini, low=0, up=10)
lme.n50l.5s.5 <- run.lme(n=50, R=5000, param=c(.5,.5), start=ini, low=1, up=100)
lme.n100l.5s.5 <- run.lme(n=100, R=5000, param=c(.5,.5), start=ini, low=1, up=100)
lme.n200l.5s.5 <- run.lme(n=200, R=5000, param=c(.5,.5), start=ini, low=1, up=100)

# LSE (Done)
lse.n20l.5s.5 <- run.lse(n=20, R= 5000, param=c(.5,.5), start=ini)
lse.n50l.5s.5 <- run.lse(n=50, R= 5000, param=c(.5,.5), start=ini)
lse.n100l.5s.5 <- run.lse(n=100, R= 5000, param=c(.5,.5), start=ini)
lse.n200l.5s.5<- run.lse(n=200, R= 5000, param=c(.5,.5), start=ini)

# WLSE (Done)
wlse.n20l.5s.5 <- run.wlse(n = 20, R= 5000, param=c(.5,.5), start=ini)
wlse.n50l.5s.5 <- run.wlse(n = 50, R= 5000, param=c(.5,.5), start=ini)
wlse.n100l.5s.5 <- run.wlse(n = 100, R= 5000, param=c(.5,.5), start=ini)
wlse.n200l.5s.5 <- run.wlse(n = 200, R= 5000, param=c(.5,.5), start=ini)

# PCE (Done)
pce.n20l.5s.5 <- run.pce(n = 20, R= 5000, param=c(.5,.5), start=ini)
pce.n50l.5s.5 <- run.pce(n = 50, R= 5000, param=c(.5,.5), start= ini)
pce.n100l.5s.5 <- run.pce(n = 100, R= 5000, param=c(.5,.5), start=ini)
pce.n200l.5s.5 <- run.pce(n = 200, R= 5000, param=c(.5,.5), start=ini)

# MPS (Done)
mps.n20l.5s.5 <- run.mps(n=20, R=5000, param=c(.5,.5), start=ini)
mps.n50l.5s.5 <- run.mps(n=50, R=5000, param=c(.5,.5), start=ini)
mps.n100l.5s.5 <- run.mps(n=100, R=5000, param=c(.5,.5), start=ini)
mps.n200l.5s.5 <- run.mps(n=200, R=5000, param=c(.5,.5), start=ini)

# CVM (Done)
cvm.n20l.5s.5 <- run.cvm(n = 20, R= 5000, param=c(.5,.5), start=ini)
cvm.n50l.5s.5 <- run.cvm(n = 50, R= 5000, param=c(.5,.5), start=ini)
cvm.n100l.5s.5 <- run.cvm(n = 100, R= 5000, param=c(.5,.5), start=ini)
cvm.n200l.5s.5 <- run.cvm(n = 200, R= 5000, param=c(.5,.5), start=ini)

# AD (Done)
ad.n20l.5s.5 <- run.ad(n = 20, R= 5000, param=c(.5,.5), start=ini)
ad.n50l.5s.5 <- run.ad(n = 50, R= 5000, param=c(.5,.5), start=ini)
ad.n100l.5s.5 <- run.ad(n = 100, R= 5000, param=c(.5,.5), start=ini)
ad.n200l.5s.5 <- run.ad(n = 200, R= 5000, param=c(.5,.5), start=ini)

# RTAD (Done)
rtad.n20l.5s.5 <- run.rtad(n = 20, R= 5000, param=c(.5,.5), start = ini)
rtad.n50l.5s.5 <- run.rtad(n = 50, R= 5000, param=c(.5,.5), start = ini)
rtad.n100l.5s.5 <- run.rtad(n = 100, R= 5000, param=c(.5,.5), start = ini)
rtad.n200l.5s.5 <- run.rtad(n = 200, R= 5000, param=c(.5,.5), start = ini)


save.image("EG June 23 2015l.5s.5.RData")

sink("output_1.txt", type = c("output", "message"))

# Combining the estimates n = 20
n20l.5s.5 <- data.frame(
  mle = c(mle.n20l.5s.5$bias.shape, mle.n20l.5s.5$rmse.shape, mle.n20l.5s.5$bias.scale, mle.n20l.5s.5$rmse.scale, mle.n20l.5s.5$Dobs, mle.n20l.5s.5$Dmax),
  mme = c(mme.n20l.5s.5$bias.shape,mme.n20l.5s.5$rmse.shape, mme.n20l.5s.5$bias.scale, mme.n20l.5s.5$rmse.scale, mme.n20l.5s.5$Dobs, mme.n20l.5s.5$Dmax),
  lse = c(lse.n20l.5s.5$bias.shape,lse.n20l.5s.5$rmse.shape, lse.n20l.5s.5$bias.scale, lse.n20l.5s.5$rmse.scale, lse.n20l.5s.5$Dobs, lse.n20l.5s.5$Dmax),
  wls = c(wlse.n20l.5s.5$bias.shape,wlse.n20l.5s.5$rmse.shape, wlse.n20l.5s.5$bias.scale, wlse.n20l.5s.5$rmse.scale, wlse.n20l.5s.5$Dobs, wlse.n20l.5s.5$Dmax),
  pce = c(pce.n20l.5s.5$bias.shape,pce.n20l.5s.5$rmse.shape, pce.n20l.5s.5$bias.scale, pce.n20l.5s.5$rmse.scale, pce.n20l.5s.5$Dobs, pce.n20l.5s.5$Dmax),
  mps = c(mps.n20l.5s.5$bias.shape,mps.n20l.5s.5$rmse.shape, mps.n20l.5s.5$bias.scale, mps.n20l.5s.5$rmse.scale, mps.n20l.5s.5$Dobs, mps.n20l.5s.5$Dmax),
  cvm = c(cvm.n20l.5s.5$bias.shape,cvm.n20l.5s.5$rmse.shape, cvm.n20l.5s.5$bias.scale, cvm.n20l.5s.5$rmse.scale, cvm.n20l.5s.5$Dobs, cvm.n20l.5s.5$Dmax),
  ad = c(ad.n20l.5s.5$bias.shape,ad.n20l.5s.5$rmse.shape, ad.n20l.5s.5$bias.scale, ad.n20l.5s.5$rmse.scale, ad.n20l.5s.5$Dobs, ad.n20l.5s.5$Dmax),
  rad = c(rtad.n20l.5s.5$bias.shape,rtad.n20l.5s.5$rmse.shape, rtad.n20l.5s.5$bias.scale, rtad.n20l.5s.5$rmse.scale, rtad.n20l.5s.5$Dobs, rtad.n20l.5s.5$Dmax)
)


n20 <- n20l.5s.5
round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))



# Combining the estimates n = 50
n50l.5s.5 <- data.frame(
  mle = c(mle.n50l.5s.5$bias.shape,mle.n50l.5s.5$rmse.shape, mle.n50l.5s.5$bias.scale, mle.n50l.5s.5$rmse.scale, mle.n50l.5s.5$Dobs, mle.n50l.5s.5$Dmax),
  mme = c(mme.n50l.5s.5$bias.shape,mme.n50l.5s.5$rmse.shape, mme.n50l.5s.5$bias.scale, mme.n50l.5s.5$rmse.scale, mme.n50l.5s.5$Dobs, mme.n50l.5s.5$Dmax),
  lse = c(lse.n50l.5s.5$bias.shape,lse.n50l.5s.5$rmse.shape, lse.n50l.5s.5$bias.scale, lse.n50l.5s.5$rmse.scale, lse.n50l.5s.5$Dobs, lse.n50l.5s.5$Dmax),
  wls = c(wlse.n50l.5s.5$bias.shape,wlse.n50l.5s.5$rmse.shape, wlse.n50l.5s.5$bias.scale, wlse.n50l.5s.5$rmse.scale, wlse.n50l.5s.5$Dobs, wlse.n50l.5s.5$Dmax),
  pce = c(pce.n50l.5s.5$bias.shape,pce.n50l.5s.5$rmse.shape, pce.n50l.5s.5$bias.scale, pce.n50l.5s.5$rmse.scale, pce.n50l.5s.5$Dobs, pce.n50l.5s.5$Dmax),
  mps = c(mps.n50l.5s.5$bias.shape,mps.n50l.5s.5$rmse.shape, mps.n50l.5s.5$bias.scale, mps.n50l.5s.5$rmse.scale, mps.n50l.5s.5$Dobs, mps.n50l.5s.5$Dmax),
  cvm = c(cvm.n50l.5s.5$bias.shape,cvm.n50l.5s.5$rmse.shape, cvm.n50l.5s.5$bias.scale, cvm.n50l.5s.5$rmse.scale, cvm.n50l.5s.5$Dobs, cvm.n50l.5s.5$Dmax),
  ad = c(ad.n50l.5s.5$bias.shape,ad.n50l.5s.5$rmse.shape, ad.n50l.5s.5$bias.scale, ad.n50l.5s.5$rmse.scale, ad.n50l.5s.5$Dobs, ad.n50l.5s.5$Dmax),
  rad = c(rtad.n50l.5s.5$bias.shape,rtad.n50l.5s.5$rmse.shape, rtad.n50l.5s.5$bias.scale, rtad.n50l.5s.5$rmse.scale, rtad.n50l.5s.5$Dobs, rtad.n50l.5s.5$Dmax)
)

n50 <- n50l.5s.5
round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))


# Combining the estimates n = 100
n100l.5s.5 <- data.frame(
  mle = c(mle.n100l.5s.5$bias.shape,mle.n100l.5s.5$rmse.shape, mle.n100l.5s.5$bias.scale, mle.n100l.5s.5$rmse.scale, mle.n100l.5s.5$Dobs, mle.n100l.5s.5$Dmax),
  mme = c(mme.n100l.5s.5$bias.shape,mme.n100l.5s.5$rmse.shape, mme.n100l.5s.5$bias.scale, mme.n100l.5s.5$rmse.scale, mme.n100l.5s.5$Dobs, mme.n100l.5s.5$Dmax),
  lse = c(lse.n100l.5s.5$bias.shape,lse.n100l.5s.5$rmse.shape, lse.n100l.5s.5$bias.scale, lse.n100l.5s.5$rmse.scale, lse.n100l.5s.5$Dobs, lse.n100l.5s.5$Dmax),
  wls = c(wlse.n100l.5s.5$bias.shape,wlse.n100l.5s.5$rmse.shape, wlse.n100l.5s.5$bias.scale, wlse.n100l.5s.5$rmse.scale, wlse.n100l.5s.5$Dobs, wlse.n100l.5s.5$Dmax),
  pce = c(pce.n100l.5s.5$bias.shape,pce.n100l.5s.5$rmse.shape, pce.n100l.5s.5$bias.scale, pce.n100l.5s.5$rmse.scale, pce.n100l.5s.5$Dobs, pce.n100l.5s.5$Dmax),
  mps = c(mps.n100l.5s.5$bias.shape,mps.n100l.5s.5$rmse.shape, mps.n100l.5s.5$bias.scale, mps.n100l.5s.5$rmse.scale, mps.n100l.5s.5$Dobs, mps.n100l.5s.5$Dmax),
  cvm = c(cvm.n100l.5s.5$bias.shape,cvm.n100l.5s.5$rmse.shape, cvm.n100l.5s.5$bias.scale, cvm.n100l.5s.5$rmse.scale, cvm.n100l.5s.5$Dobs, cvm.n100l.5s.5$Dmax),
  ad = c(ad.n100l.5s.5$bias.shape,ad.n100l.5s.5$rmse.shape, ad.n100l.5s.5$bias.scale, ad.n100l.5s.5$rmse.scale, ad.n100l.5s.5$Dobs, ad.n100l.5s.5$Dmax),
  rad = c(rtad.n100l.5s.5$bias.shape,rtad.n100l.5s.5$rmse.shape, rtad.n100l.5s.5$bias.scale, rtad.n100l.5s.5$rmse.scale, rtad.n100l.5s.5$Dobs, rtad.n100l.5s.5$Dmax)
)

n100 <- n100l.5s.5
round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))

# Combining the estimates n = 200
n200l.5s.5 <- data.frame(
  mle = c(mle.n200l.5s.5$bias.shape,mle.n200l.5s.5$rmse.shape, mle.n200l.5s.5$bias.scale, mle.n200l.5s.5$rmse.scale, mle.n200l.5s.5$Dobs, mle.n200l.5s.5$Dmax),
  mme = c(mme.n200l.5s.5$bias.shape,mme.n200l.5s.5$rmse.shape, mme.n200l.5s.5$bias.scale, mme.n200l.5s.5$rmse.scale, mme.n200l.5s.5$Dobs, mme.n200l.5s.5$Dmax),
  lse = c(lse.n200l.5s.5$bias.shape,lse.n200l.5s.5$rmse.shape, lse.n200l.5s.5$bias.scale, lse.n200l.5s.5$rmse.scale, lse.n200l.5s.5$Dobs, lse.n200l.5s.5$Dmax),
  wls = c(wlse.n200l.5s.5$bias.shape,wlse.n200l.5s.5$rmse.shape, wlse.n200l.5s.5$bias.scale, wlse.n200l.5s.5$rmse.scale, wlse.n200l.5s.5$Dobs, wlse.n200l.5s.5$Dmax),
  pce = c(pce.n200l.5s.5$bias.shape,pce.n200l.5s.5$rmse.shape, pce.n200l.5s.5$bias.scale, pce.n200l.5s.5$rmse.scale, pce.n200l.5s.5$Dobs, pce.n200l.5s.5$Dmax),
  mps = c(mps.n200l.5s.5$bias.shape,mps.n200l.5s.5$rmse.shape, mps.n200l.5s.5$bias.scale, mps.n200l.5s.5$rmse.scale, mps.n200l.5s.5$Dobs, mps.n200l.5s.5$Dmax),
  cvm = c(cvm.n200l.5s.5$bias.shape,cvm.n200l.5s.5$rmse.shape, cvm.n200l.5s.5$bias.scale, cvm.n200l.5s.5$rmse.scale, cvm.n200l.5s.5$Dobs, cvm.n200l.5s.5$Dmax),
  ad = c(ad.n200l.5s.5$bias.shape,ad.n200l.5s.5$rmse.shape, ad.n200l.5s.5$bias.scale, ad.n200l.5s.5$rmse.scale, ad.n200l.5s.5$Dobs, ad.n200l.5s.5$Dmax),
  rad = c(rtad.n200l.5s.5$bias.shape,rtad.n200l.5s.5$rmse.shape, rtad.n200l.5s.5$bias.scale, rtad.n200l.5s.5$rmse.scale, rtad.n200l.5s.5$Dobs, rtad.n200l.5s.5$Dmax)
)

n200 <- n200l.5s.5
round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))

sink()



# End #

#################################################
############ RUNNING SIMULATION #################
############ lambda = 1, sigma =.5 ##############
#################################################

par = c(1, .5)
ini = c(1.1, .4)

# MLE (Done for EG)

mle.n20l1.1s.5 <- run.mle(n=20, R=5000, param=par, start=ini)
mle.n50l1.1s.5 <- run.mle(n=50, R=5000, param=par, start=ini)
mle.n100l1.1s.5 <- run.mle(n=100, R=5000, param=par, start=ini)
mle.n200l1.1s.5 <- run.mle(n=200, R=5000, param=par, start=ini)


# MME (Done for EG)

mme.n20l1.1s.5 <- run.mme(n=20, R=5000, param=par, start=ini,low=-5, up=30)
mme.n50l1.1s.5 <- run.mme(n=50, R=5000, param=par, start=ini,low=-5, up=30)
mme.n100l1.1s.5 <- run.mme(n=100, R=5000, param=par, start=ini,low=-5, up=30)
mme.n200l1.1s.5 <- run.mme(n=200, R=5000, param=par, start=ini,low=-5, up=30)

# # LME (Does not work)
# 
# lme.n20l1.1s.5 <- run.lme(n=20, R=10000, param=c(1,.5), start=ini,low=-5, up=30)
# lme.n50l1.1s.5 <- run.lme(n = 50, R= 10000, param=c(1,.5), start=ini)
# lme.n100l1.1s.5 <- run.lme(n = 100, R= 10000, param=c(1,.5), start=ini)
# lme.n200l1.1s.5 <- run.lme(n = 200, R= 10000, param=c(1,.5), start=ini)

# LSE (Done with EG)

lse.n20l1.1s.5 <- run.lse(n = 20, R= 5000, param=par, start=ini)
lse.n50l1.1s.5 <- run.lse(n = 50, R= 5000, param=par, start=ini)
lse.n100l1.1s.5 <- run.lse(n = 100, R= 5000, param=par, start=ini)
lse.n200l1.1s.5 <- run.lse(n = 200, R= 5000, param=par, start=ini)

# WLSE (Done with EG)

wlse.n20l1.1s.5 <- run.wlse(n = 20, R=5000, param=par, start=ini)
wlse.n50l1.1s.5 <- run.wlse(n = 50, R=5000, param=par, start=ini)
wlse.n100l1.1s.5 <- run.wlse(n = 100, R=5000, param=par, start=ini)
wlse.n200l1.1s.5 <- run.wlse(n = 200, R=5000, param=par, start=ini)

# PCE (Done with EG)

pce.n20l1.1s.5 <- run.pce(n= 20, R= 5000, param=par, start=ini)
pce.n50l1.1s.5 <- run.pce(n= 50, R= 5000, param=par, start= ini)
pce.n100l1.1s.5 <- run.pce(n=100, R=5000, param=par, start=ini)
pce.n200l1.1s.5 <- run.pce(n= 200, R=5000, param=par, start=ini)

# MPS (Done with EG)

mps.n20l1.1s.5 <- run.mps(n = 20, R= 5000, param=par, start=ini)
mps.n50l1.1s.5 <- run.mps(n = 50, R= 5000, param=par, start=ini)
mps.n100l1.1s.5 <- run.mps(n = 100, R=5000, param=par, start=ini)
mps.n200l1.1s.5 <- run.mps(n = 200, R=5000, param=par, start=ini)

# CVM (Done with EG)

cvm.n20l1.1s.5 <- run.cvm(n = 20, R= 5000, param=par, start=ini)
cvm.n50l1.1s.5 <- run.cvm(n = 50, R= 5000, param=par, start=ini)
cvm.n100l1.1s.5 <- run.cvm(n = 100, R= 5000, param= par, start=ini)
cvm.n200l1.1s.5 <- run.cvm(n = 200, R= 5000, param=par, start=ini)

# AD (Done with EG)

ad.n20l1.1s.5 <- run.ad(n = 20, R= 5000, param=par, start=ini)
ad.n50l1.1s.5 <- run.ad(n = 50, R= 5000, param=par, start=ini)
ad.n100l1.1s.5 <- run.ad(n = 100, R= 5000, param=par, start=ini)
ad.n200l1.1s.5 <- run.ad(n = 200, R= 5000, param=par, start=ini)

# RTAD (Done with EG)

rtad.n20l1.1s.5 <- run.rtad(n = 20, R=5000, param=par, start = ini)
rtad.n50l1.1s.5 <- run.rtad(n = 50, R=5000, param=par, start = ini)
rtad.n100l1.1s.5 <- run.rtad(n = 100, R=5000, param=par, start = ini)
rtad.n200l1.1s.5 <- run.rtad(n = 200, R=5000, param=par, start = ini)


save.image("EG JUne 24 2015.l1.1s.5.RData")

sink("output_2.txt", type = c("output", "message"))

# Combining the estimates n = 20
n20l1.1s.5 <- data.frame(
  mle = c(mle.n20l1.1s.5$bias.shape,mle.n20l1.1s.5$rmse.shape, mle.n20l1.1s.5$bias.scale, mle.n20l1.1s.5$rmse.scale, mle.n20l1.1s.5$Dobs, mle.n20l1.1s.5$Dmax),
  mme = c(mme.n20l1.1s.5$bias.shape,mme.n20l1.1s.5$rmse.shape, mme.n20l1.1s.5$bias.scale, mme.n20l1.1s.5$rmse.scale, mme.n20l1.1s.5$Dobs, mme.n20l1.1s.5$Dmax),
  lse = c(lse.n20l1.1s.5$bias.shape,lse.n20l1.1s.5$rmse.shape, lse.n20l1.1s.5$bias.scale, lse.n20l1.1s.5$rmse.scale, lse.n20l1.1s.5$Dobs, lse.n20l1.1s.5$Dmax),
  wls = c(wlse.n20l1.1s.5$bias.shape,wlse.n20l1.1s.5$rmse.shape, wlse.n20l1.1s.5$bias.scale, wlse.n20l1.1s.5$rmse.scale, wlse.n20l1.1s.5$Dobs, wlse.n20l1.1s.5$Dmax),
  pce = c(pce.n20l1.1s.5$bias.shape,pce.n20l1.1s.5$rmse.shape, pce.n20l1.1s.5$bias.scale, pce.n20l1.1s.5$rmse.scale, pce.n20l1.1s.5$Dobs, pce.n20l1.1s.5$Dmax),
  mps = c(mps.n20l1.1s.5$bias.shape,mps.n20l1.1s.5$rmse.shape, mps.n20l1.1s.5$bias.scale, mps.n20l1.1s.5$rmse.scale, mps.n20l1.1s.5$Dobs, mps.n20l1.1s.5$Dmax),
  cvm = c(cvm.n20l1.1s.5$bias.shape,cvm.n20l1.1s.5$rmse.shape, cvm.n20l1.1s.5$bias.scale, cvm.n20l1.1s.5$rmse.scale, cvm.n20l1.1s.5$Dobs, cvm.n20l1.1s.5$Dmax),
  ad = c(ad.n20l1.1s.5$bias.shape,ad.n20l1.1s.5$rmse.shape, ad.n20l1.1s.5$bias.scale, ad.n20l1.1s.5$rmse.scale, ad.n20l1.1s.5$Dobs, ad.n20l1.1s.5$Dmax),
  rad = c(rtad.n20l1.1s.5$bias.shape,rtad.n20l1.1s.5$rmse.shape, rtad.n20l1.1s.5$bias.scale, rtad.n20l1.1s.5$rmse.scale, rtad.n20l1.1s.5$Dobs, rtad.n20l1.1s.5$Dmax)
)

n20 <- n20l1.1s.5
round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))


# Combining the estimates n = 50
n50l1.1s.5 <- data.frame(
  mle = c(mle.n50l1.1s.5$bias.shape,mle.n50l1.1s.5$rmse.shape, mle.n50l1.1s.5$bias.scale, mle.n50l1.1s.5$rmse.scale, mle.n50l1.1s.5$Dobs, mle.n50l1.1s.5$Dmax),
  mme = c(mme.n50l1.1s.5$bias.shape,mme.n50l1.1s.5$rmse.shape, mme.n50l1.1s.5$bias.scale, mme.n50l1.1s.5$rmse.scale, mme.n50l1.1s.5$Dobs, mme.n50l1.1s.5$Dmax),
  lse = c(lse.n50l1.1s.5$bias.shape,lse.n50l1.1s.5$rmse.shape, lse.n50l1.1s.5$bias.scale, lse.n50l1.1s.5$rmse.scale, lse.n50l1.1s.5$Dobs, lse.n50l1.1s.5$Dmax),
  wls = c(wlse.n50l1.1s.5$bias.shape,wlse.n50l1.1s.5$rmse.shape, wlse.n50l1.1s.5$bias.scale, wlse.n50l1.1s.5$rmse.scale, wlse.n50l1.1s.5$Dobs, wlse.n50l1.1s.5$Dmax),
  pce = c(pce.n50l1.1s.5$bias.shape,pce.n50l1.1s.5$rmse.shape, pce.n50l1.1s.5$bias.scale, pce.n50l1.1s.5$rmse.scale, pce.n50l1.1s.5$Dobs, pce.n50l1.1s.5$Dmax),
  mps = c(mps.n50l1.1s.5$bias.shape,mps.n50l1.1s.5$rmse.shape, mps.n50l1.1s.5$bias.scale, mps.n50l1.1s.5$rmse.scale, mps.n50l1.1s.5$Dobs, mps.n50l1.1s.5$Dmax),
  cvm = c(cvm.n50l1.1s.5$bias.shape,cvm.n50l1.1s.5$rmse.shape, cvm.n50l1.1s.5$bias.scale, cvm.n50l1.1s.5$rmse.scale, cvm.n50l1.1s.5$Dobs, cvm.n50l1.1s.5$Dmax),
  ad = c(ad.n50l1.1s.5$bias.shape,ad.n50l1.1s.5$rmse.shape, ad.n50l1.1s.5$bias.scale, ad.n50l1.1s.5$rmse.scale, ad.n50l1.1s.5$Dobs, ad.n50l1.1s.5$Dmax),
  rad = c(rtad.n50l1.1s.5$bias.shape,rtad.n50l1.1s.5$rmse.shape, rtad.n50l1.1s.5$bias.scale, rtad.n50l1.1s.5$rmse.scale, rtad.n50l1.1s.5$Dobs, rtad.n50l1.1s.5$Dmax)
)

n50 <- n50l1.1s.5
round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))


# Combining the estimates n = 100
n100l1.1s.5 <- data.frame(
  mle = c(mle.n100l1.1s.5$bias.shape,mle.n100l1.1s.5$rmse.shape, mle.n100l1.1s.5$bias.scale, mle.n100l1.1s.5$rmse.scale, mle.n100l1.1s.5$Dobs, mle.n100l1.1s.5$Dmax),
  mme = c(mme.n100l1.1s.5$bias.shape,mme.n100l1.1s.5$rmse.shape, mme.n100l1.1s.5$bias.scale, mme.n100l1.1s.5$rmse.scale, mme.n100l1.1s.5$Dobs, mme.n100l1.1s.5$Dmax),
  lse = c(lse.n100l1.1s.5$bias.shape,lse.n100l1.1s.5$rmse.shape, lse.n100l1.1s.5$bias.scale, lse.n100l1.1s.5$rmse.scale, lse.n100l1.1s.5$Dobs, lse.n100l1.1s.5$Dmax),
  wls = c(wlse.n100l1.1s.5$bias.shape,wlse.n100l1.1s.5$rmse.shape, wlse.n100l1.1s.5$bias.scale, wlse.n100l1.1s.5$rmse.scale, wlse.n100l1.1s.5$Dobs, wlse.n100l1.1s.5$Dmax),
  pce = c(pce.n100l1.1s.5$bias.shape,pce.n100l1.1s.5$rmse.shape, pce.n100l1.1s.5$bias.scale, pce.n100l1.1s.5$rmse.scale, pce.n100l1.1s.5$Dobs, pce.n100l1.1s.5$Dmax),
  mps = c(mps.n100l1.1s.5$bias.shape,mps.n100l1.1s.5$rmse.shape, mps.n100l1.1s.5$bias.scale, mps.n100l1.1s.5$rmse.scale, mps.n100l1.1s.5$Dobs, mps.n100l1.1s.5$Dmax),
  cvm = c(cvm.n100l1.1s.5$bias.shape,cvm.n100l1.1s.5$rmse.shape, cvm.n100l1.1s.5$bias.scale, cvm.n100l1.1s.5$rmse.scale, cvm.n100l1.1s.5$Dobs, cvm.n100l1.1s.5$Dmax),
  ad = c(ad.n100l1.1s.5$bias.shape,ad.n100l1.1s.5$rmse.shape, ad.n100l1.1s.5$bias.scale, ad.n100l1.1s.5$rmse.scale, ad.n100l1.1s.5$Dobs, ad.n100l1.1s.5$Dmax),
  rad = c(rtad.n100l1.1s.5$bias.shape,rtad.n100l1.1s.5$rmse.shape, rtad.n100l1.1s.5$bias.scale, rtad.n100l1.1s.5$rmse.scale, rtad.n100l1.1s.5$Dobs, rtad.n100l1.1s.5$Dmax)
)

n100 <- n100l1.1s.5
round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))

# Combining the estimates n = 200
n200l1.1s.5 <- data.frame(
  mle = c(mle.n200l1.1s.5$bias.shape,mle.n200l1.1s.5$rmse.shape, mle.n200l1.1s.5$bias.scale, mle.n200l1.1s.5$rmse.scale, mle.n200l1.1s.5$Dobs, mle.n200l1.1s.5$Dmax),
  mme = c(mme.n200l1.1s.5$bias.shape,mme.n200l1.1s.5$rmse.shape, mme.n200l1.1s.5$bias.scale, mme.n200l1.1s.5$rmse.scale, mme.n200l1.1s.5$Dobs, mme.n200l1.1s.5$Dmax),
  lse = c(lse.n200l1.1s.5$bias.shape,lse.n200l1.1s.5$rmse.shape, lse.n200l1.1s.5$bias.scale, lse.n200l1.1s.5$rmse.scale, lse.n200l1.1s.5$Dobs, lse.n200l1.1s.5$Dmax),
  wls = c(wlse.n200l1.1s.5$bias.shape,wlse.n200l1.1s.5$rmse.shape, wlse.n200l1.1s.5$bias.scale, wlse.n200l1.1s.5$rmse.scale, wlse.n200l1.1s.5$Dobs, wlse.n200l1.1s.5$Dmax),
  pce = c(pce.n200l1.1s.5$bias.shape,pce.n200l1.1s.5$rmse.shape, pce.n200l1.1s.5$bias.scale, pce.n200l1.1s.5$rmse.scale, pce.n200l1.1s.5$Dobs, pce.n200l1.1s.5$Dmax),
  mps = c(mps.n200l1.1s.5$bias.shape,mps.n200l1.1s.5$rmse.shape, mps.n200l1.1s.5$bias.scale, mps.n200l1.1s.5$rmse.scale, mps.n200l1.1s.5$Dobs, mps.n200l1.1s.5$Dmax),
  cvm = c(cvm.n200l1.1s.5$bias.shape,cvm.n200l1.1s.5$rmse.shape, cvm.n200l1.1s.5$bias.scale, cvm.n200l1.1s.5$rmse.scale, cvm.n200l1.1s.5$Dobs, cvm.n200l1.1s.5$Dmax),
  ad = c(ad.n200l1.1s.5$bias.shape,ad.n200l1.1s.5$rmse.shape, ad.n200l1.1s.5$bias.scale, ad.n200l1.1s.5$rmse.scale, ad.n200l1.1s.5$Dobs, ad.n200l1.1s.5$Dmax),
  rad = c(rtad.n200l1.1s.5$bias.shape,rtad.n200l1.1s.5$rmse.shape, rtad.n200l1.1s.5$bias.scale, rtad.n200l1.1s.5$rmse.scale, rtad.n200l1.1s.5$Dobs, rtad.n200l1.1s.5$Dmax)
)

n200 <- n200l1.1s.5
round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))


sink()


# End #

#################################################
############ RUNNING SIMULATION #################
############ lambda = .5, sigma=3 ##############
#################################################

ini = c(.5, 2.5)
par=c(.5, 3)

# MLE (Done for EG)
mle.n20l.5s3 <-run.mle(n = 20, R= 5000, param=par, start=ini)
mle.n50l.5s3  <- run.mle(n = 50, R= 5000, param=par, start=ini)
mle.n100l.5s3  <- run.mle(n = 100, R= 5000, param=par, start=ini)
mle.n200l.5s3  <- run.mle(n = 200, R= 5000, param=par, start=ini)


# MME (Done for EG)
mme.n20l.5s3<-run.mme(n=20,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n50l.5s3<-run.mme(n=50,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n100l.5s3<-run.mme(n=100,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n200l.5s3<-run.mme(n=200,R=5000, param=par, start=ini,low=1, up=Inf)

# LME (Does not work for EG)
# lme.n20l.5s3  <- run.lme(n = 20,R=5000, param=par, start=ini,low=1, up=Inf)
# lme.n50l.5s3  <- run.lme(n = 50, R= 5000, param=par, start=ini)
# lme.n100l.5s3  <- run.lme(n = 100, R= 5000, param=par, start=ini)
# lme.n200l.5s3  <- run.lme(n = 200, R= 5000, param=par, start=ini)

# LSE (Done for EG)
lse.n20l.5s3  <- run.lse(n = 20, R= 5000, param=par, start=ini)
lse.n50l.5s3  <- run.lse(n = 50, R= 5000, param=par, start=ini)
lse.n100l.5s3  <- run.lse(n = 100, R= 5000, param=par, start=ini)
lse.n200l.5s3  <- run.lse(n = 200, R= 5000, param=par, start=ini)

# WLSE(Done for EG)
wlse.n20l.5s3  <- run.wlse(n = 20, R= 5000, param=par, start=ini)
wlse.n50l.5s3  <- run.wlse(n = 50, R= 5000, param=par, start=ini)
wlse.n100l.5s3  <- run.wlse(n = 100, R= 5000, param=par, start=ini)
wlse.n200l.5s3  <- run.wlse(n = 200, R= 5000, param=par, start=ini)

# PCE(Done for EG)
pce.n20l.5s3  <- run.pce(n = 20, R= 5000, param=par, start=ini)
pce.n50l.5s3  <- run.pce(n = 50, R= 5000, param=par, start= ini)
pce.n100l.5s3  <- run.pce(n = 100, R= 5000, param=par, start=ini)
pce.n200l.5s3  <- run.pce(n = 200, R= 5000, param=par, start=ini)

# MPS(Done for EG)
mps.n20l.5s3  <- run.mps(n = 20, R= 5000, param=par, start=ini)
mps.n50l.5s3  <- run.mps(n = 50, R= 5000, param=par, start=ini)
mps.n100l.5s3  <- run.mps(n = 100, R= 5000, param=par, start=ini)
mps.n200l.5s3  <- run.mps(n = 200, R= 5000, param=par, start=ini)

# CVM(Done for EG)
cvm.n20l.5s3  <- run.cvm(n = 20, R= 5000, param=par, start=ini)
cvm.n50l.5s3  <- run.cvm(n = 50, R= 5000, param=par, start=ini)
cvm.n100l.5s3  <- run.cvm(n = 100, R= 5000, param=par, start=ini)
cvm.n200l.5s3  <- run.cvm(n = 200, R= 5000, param=par, start=ini)

# AD(Done for EG)
ad.n20l.5s3  <- run.ad(n = 20, R= 5000, param=par, start=ini)
ad.n50l.5s3  <- run.ad(n = 50, R= 5000, param=par, start=ini)
ad.n100l.5s3  <- run.ad(n = 100, R= 5000, param=par, start=ini)
ad.n200l.5s3  <- run.ad(n = 200, R= 5000, param=par, start=ini)

# RTAD(Done for EG)
rtad.n20l.5s3  <- run.rtad(n = 20, R= 5000, param=par, start = ini)
rtad.n50l.5s3  <- run.rtad(n = 50, R= 5000, param=par, start = ini)
rtad.n100l.5s3  <- run.rtad(n = 100, R= 5000, param=par, start = ini)
rtad.n200l.5s3  <- run.rtad(n = 200, R= 5000, param=par, start = ini)


save.image("June 28 2015.RData")

sink("output_3.txt", type = c("output", "message"))

# Combining the estimates n = 20
n20 <- data.frame(
  mle = c(mle.n20l.5s3$bias.shape,mle.n20l.5s3$rmse.shape, mle.n20l.5s3$bias.scale, mle.n20l.5s3$rmse.scale, mle.n20l.5s3$Dobs, mle.n20l.5s3$Dmax),
  mme = c(mme.n20l.5s3$bias.shape,mme.n20l.5s3$rmse.shape, mme.n20l.5s3$bias.scale, mme.n20l.5s3$rmse.scale, mme.n20l.5s3$Dobs, mme.n20l.5s3$Dmax),
  lse = c(lse.n20l.5s3$bias.shape,lse.n20l.5s3$rmse.shape, lse.n20l.5s3$bias.scale, lse.n20l.5s3$rmse.scale, lse.n20l.5s3$Dobs, lse.n20l.5s3$Dmax),
  wls = c(wlse.n20l.5s3$bias.shape,wlse.n20l.5s3$rmse.shape, wlse.n20l.5s3$bias.scale, wlse.n20l.5s3$rmse.scale, wlse.n20l.5s3$Dobs, wlse.n20l.5s3$Dmax),
  pce = c(pce.n20l.5s3$bias.shape,pce.n20l.5s3$rmse.shape, pce.n20l.5s3$bias.scale, pce.n20l.5s3$rmse.scale, pce.n20l.5s3$Dobs, pce.n20l.5s3$Dmax),
  mps = c(mps.n20l.5s3$bias.shape,mps.n20l.5s3$rmse.shape, mps.n20l.5s3$bias.scale, mps.n20l.5s3$rmse.scale, mps.n20l.5s3$Dobs, mps.n20l.5s3$Dmax),
  cvm = c(cvm.n20l.5s3$bias.shape,cvm.n20l.5s3$rmse.shape, cvm.n20l.5s3$bias.scale, cvm.n20l.5s3$rmse.scale, cvm.n20l.5s3$Dobs, cvm.n20l.5s3$Dmax),
  ad = c(ad.n20l.5s3$bias.shape,ad.n20l.5s3$rmse.shape, ad.n20l.5s3$bias.scale, ad.n20l.5s3$rmse.scale, ad.n20l.5s3$Dobs, ad.n20l.5s3$Dmax),
  rad = c(rtad.n20l.5s3$bias.shape,rtad.n20l.5s3$rmse.shape, rtad.n20l.5s3$bias.scale, rtad.n20l.5s3$rmse.scale, rtad.n20l.5s3$Dobs, rtad.n20l.5s3$Dmax)
)

round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))


# Combining the estimates n = 50
n50 <- data.frame(
  mle = c(mle.n50l.5s3$bias.shape,mle.n50l.5s3$rmse.shape, mle.n50l.5s3$bias.scale, mle.n50l.5s3$rmse.scale, mle.n50l.5s3$Dobs, mle.n50l.5s3$Dmax),
  mme = c(mme.n50l.5s3$bias.shape,mme.n50l.5s3$rmse.shape, mme.n50l.5s3$bias.scale, mme.n50l.5s3$rmse.scale, mme.n50l.5s3$Dobs, mme.n50l.5s3$Dmax),
  lse = c(lse.n50l.5s3$bias.shape,lse.n50l.5s3$rmse.shape, lse.n50l.5s3$bias.scale, lse.n50l.5s3$rmse.scale, lse.n50l.5s3$Dobs, lse.n50l.5s3$Dmax),
  wls = c(wlse.n50l.5s3$bias.shape,wlse.n50l.5s3$rmse.shape, wlse.n50l.5s3$bias.scale, wlse.n50l.5s3$rmse.scale, wlse.n50l.5s3$Dobs, wlse.n50l.5s3$Dmax),
  pce = c(pce.n50l.5s3$bias.shape,pce.n50l.5s3$rmse.shape, pce.n50l.5s3$bias.scale, pce.n50l.5s3$rmse.scale, pce.n50l.5s3$Dobs, pce.n50l.5s3$Dmax),
  mps = c(mps.n50l.5s3$bias.shape,mps.n50l.5s3$rmse.shape, mps.n50l.5s3$bias.scale, mps.n50l.5s3$rmse.scale, mps.n50l.5s3$Dobs, mps.n50l.5s3$Dmax),
  cvm = c(cvm.n50l.5s3$bias.shape,cvm.n50l.5s3$rmse.shape, cvm.n50l.5s3$bias.scale, cvm.n50l.5s3$rmse.scale, cvm.n50l.5s3$Dobs, cvm.n50l.5s3$Dmax),
  ad = c(ad.n50l.5s3$bias.shape,ad.n50l.5s3$rmse.shape, ad.n50l.5s3$bias.scale, ad.n50l.5s3$rmse.scale, ad.n50l.5s3$Dobs, ad.n50l.5s3$Dmax),
  rad = c(rtad.n50l.5s3$bias.shape,rtad.n50l.5s3$rmse.shape, rtad.n50l.5s3$bias.scale, rtad.n50l.5s3$rmse.scale, rtad.n50l.5s3$Dobs, rtad.n50l.5s3$Dmax)
)

round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))


# Combining the estimates n = 100
n100 <- data.frame(
  mle = c(mle.n100l.5s3$bias.shape,mle.n100l.5s3$rmse.shape, mle.n100l.5s3$bias.scale, mle.n100l.5s3$rmse.scale, mle.n100l.5s3$Dobs, mle.n100l.5s3$Dmax),
  mme = c(mme.n100l.5s3$bias.shape,mme.n100l.5s3$rmse.shape, mme.n100l.5s3$bias.scale, mme.n100l.5s3$rmse.scale, mme.n100l.5s3$Dobs, mme.n100l.5s3$Dmax),
  lse = c(lse.n100l.5s3$bias.shape,lse.n100l.5s3$rmse.shape, lse.n100l.5s3$bias.scale, lse.n100l.5s3$rmse.scale, lse.n100l.5s3$Dobs, lse.n100l.5s3$Dmax),
  wls = c(wlse.n100l.5s3$bias.shape,wlse.n100l.5s3$rmse.shape, wlse.n100l.5s3$bias.scale, wlse.n100l.5s3$rmse.scale, wlse.n100l.5s3$Dobs, wlse.n100l.5s3$Dmax),
  pce = c(pce.n100l.5s3$bias.shape,pce.n100l.5s3$rmse.shape, pce.n100l.5s3$bias.scale, pce.n100l.5s3$rmse.scale, pce.n100l.5s3$Dobs, pce.n100l.5s3$Dmax),
  mps = c(mps.n100l.5s3$bias.shape,mps.n100l.5s3$rmse.shape, mps.n100l.5s3$bias.scale, mps.n100l.5s3$rmse.scale, mps.n100l.5s3$Dobs, mps.n100l.5s3$Dmax),
  cvm = c(cvm.n100l.5s3$bias.shape,cvm.n100l.5s3$rmse.shape, cvm.n100l.5s3$bias.scale, cvm.n100l.5s3$rmse.scale, cvm.n100l.5s3$Dobs, cvm.n100l.5s3$Dmax),
  ad = c(ad.n100l.5s3$bias.shape,ad.n100l.5s3$rmse.shape, ad.n100l.5s3$bias.scale, ad.n100l.5s3$rmse.scale, ad.n100l.5s3$Dobs, ad.n100l.5s3$Dmax),
  rad = c(rtad.n100l.5s3$bias.shape,rtad.n100l.5s3$rmse.shape, rtad.n100l.5s3$bias.scale, rtad.n100l.5s3$rmse.scale, rtad.n100l.5s3$Dobs, rtad.n100l.5s3$Dmax)
)

round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))

# Combining the estimates n = 200
n200 <- data.frame(
  mle = c(mle.n200l.5s3$bias.shape,mle.n200l.5s3$rmse.shape, mle.n200l.5s3$bias.scale, mle.n200l.5s3$rmse.scale, mle.n200l.5s3$Dobs, mle.n200l.5s3$Dmax),
  mme = c(mme.n200l.5s3$bias.shape,mme.n200l.5s3$rmse.shape, mme.n200l.5s3$bias.scale, mme.n200l.5s3$rmse.scale, mme.n200l.5s3$Dobs, mme.n200l.5s3$Dmax),
  lse = c(lse.n200l.5s3$bias.shape,lse.n200l.5s3$rmse.shape, lse.n200l.5s3$bias.scale, lse.n200l.5s3$rmse.scale, lse.n200l.5s3$Dobs, lse.n200l.5s3$Dmax),
  wls = c(wlse.n200l.5s3$bias.shape,wlse.n200l.5s3$rmse.shape, wlse.n200l.5s3$bias.scale, wlse.n200l.5s3$rmse.scale, wlse.n200l.5s3$Dobs, wlse.n200l.5s3$Dmax),
  pce = c(pce.n200l.5s3$bias.shape,pce.n200l.5s3$rmse.shape, pce.n200l.5s3$bias.scale, pce.n200l.5s3$rmse.scale, pce.n200l.5s3$Dobs, pce.n200l.5s3$Dmax),
  mps = c(mps.n200l.5s3$bias.shape,mps.n200l.5s3$rmse.shape, mps.n200l.5s3$bias.scale, mps.n200l.5s3$rmse.scale, mps.n200l.5s3$Dobs, mps.n200l.5s3$Dmax),
  cvm = c(cvm.n200l.5s3$bias.shape,cvm.n200l.5s3$rmse.shape, cvm.n200l.5s3$bias.scale, cvm.n200l.5s3$rmse.scale, cvm.n200l.5s3$Dobs, cvm.n200l.5s3$Dmax),
  ad = c(ad.n200l.5s3$bias.shape,ad.n200l.5s3$rmse.shape, ad.n200l.5s3$bias.scale, ad.n200l.5s3$rmse.scale, ad.n200l.5s3$Dobs, ad.n200l.5s3$Dmax),
  rad = c(rtad.n200l.5s3$bias.shape,rtad.n200l.5s3$rmse.shape, rtad.n200l.5s3$bias.scale, rtad.n200l.5s3$rmse.scale, rtad.n200l.5s3$Dobs, rtad.n200l.5s3$Dmax)
)


round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))


sink()


# End #


#################################################
############ RUNNING SIMULATION #################
############ lambda = 1, sigma=3 ##############
#################################################

ini = c(.8, 2.5)
par=c(1, 3)

# MLE (Done for EG)
mle.n20l1s3 <-run.mle(n = 20, R= 5000, param=par, start=ini)
mle.n50l1s3  <- run.mle(n = 50, R= 5000, param=par, start=ini)
mle.n100l1s3  <- run.mle(n = 100, R= 5000, param=par, start=ini)
mle.n200l1s3  <- run.mle(n = 200, R= 5000, param=par, start=ini)


# MME 
mme.n20l1s3<-run.mme(n=20,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n50l1s3<-run.mme(n=50,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n100l1s3<-run.mme(n=100,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n200l1s3<-run.mme(n=200,R=5000, param=par, start=ini,low=1, up=Inf)

# LME (Does not work for EG)
# lme.n20l1s3  <- run.lme(n = 20,R=5000, param=par, start=ini,low=1, up=Inf)
# lme.n50l1s3  <- run.lme(n = 50, R= 5000, param=par, start=ini)
# lme.n100l1s3  <- run.lme(n = 100, R= 5000, param=par, start=ini)
# lme.n200l1s3  <- run.lme(n = 200, R= 5000, param=par, start=ini)

# LSE 
lse.n20l1s3  <- run.lse(n = 20, R= 5000, param=par, start=ini)
lse.n50l1s3  <- run.lse(n = 50, R= 5000, param=par, start=ini)
lse.n100l1s3  <- run.lse(n = 100, R= 5000, param=par, start=ini)
lse.n200l1s3  <- run.lse(n = 200, R= 5000, param=par, start=ini)

# WLSE
wlse.n20l1s3  <- run.wlse(n = 20, R= 5000, param=par, start=ini)
wlse.n50l1s3  <- run.wlse(n = 50, R= 5000, param=par, start=ini)
wlse.n100l1s3  <- run.wlse(n = 100, R= 5000, param=par, start=ini)
wlse.n200l1s3  <- run.wlse(n = 200, R= 5000, param=par, start=ini)

# PCE
pce.n20l1s3  <- run.pce(n = 20, R= 5000, param=par, start=ini)
pce.n50l1s3  <- run.pce(n = 50, R= 5000, param=par, start= ini)
pce.n100l1s3  <- run.pce(n = 100, R= 5000, param=par, start=ini)
pce.n200l1s3  <- run.pce(n = 200, R= 5000, param=par, start=ini)

# MPS
mps.n20l1s3  <- run.mps(n = 20, R= 5000, param=par, start=ini)
mps.n50l1s3  <- run.mps(n = 50, R= 5000, param=par, start=ini)
mps.n100l1s3  <- run.mps(n = 100, R= 5000, param=par, start=ini)
mps.n200l1s3  <- run.mps(n = 200, R= 5000, param=par, start=ini)

# CVM
cvm.n20l1s3  <- run.cvm(n = 20, R= 5000, param=par, start=ini)
cvm.n50l1s3  <- run.cvm(n = 50, R= 5000, param=par, start=ini)
cvm.n100l1s3  <- run.cvm(n = 100, R= 5000, param=par, start=ini)
cvm.n200l1s3  <- run.cvm(n = 200, R= 5000, param=par, start=ini)

# AD
ad.n20l1s3  <- run.ad(n = 20, R= 5000, param=par, start=ini)
ad.n50l1s3  <- run.ad(n = 50, R= 5000, param=par, start=ini)
ad.n100l1s3  <- run.ad(n = 100, R= 5000, param=par, start=ini)
ad.n200l1s3  <- run.ad(n = 200, R= 5000, param=par, start=ini)

# RTAD
rtad.n20l1s3  <- run.rtad(n = 20, R= 5000, param=par, start = ini)
rtad.n50l1s3  <- run.rtad(n = 50, R= 5000, param=par, start = ini)
rtad.n100l1s3  <- run.rtad(n = 100, R= 5000, param=par, start = ini)
rtad.n200l1s3  <- run.rtad(n = 200, R= 5000, param=par, start = ini)


save.image("June 28 2015.RData")

sink("output_4.txt", type = c("output", "message"))

# Combining the estimates n = 20
n20l1s3 <- data.frame(
  mle = c(mle.n20l1s3$bias.shape,mle.n20l1s3$rmse.shape, mle.n20l1s3$bias.scale, mle.n20l1s3$rmse.scale, mle.n20l1s3$Dobs, mle.n20l1s3$Dmax),
  mme = c(mme.n20l1s3$bias.shape,mme.n20l1s3$rmse.shape, mme.n20l1s3$bias.scale, mme.n20l1s3$rmse.scale, mme.n20l1s3$Dobs, mme.n20l1s3$Dmax),
  lse = c(lse.n20l1s3$bias.shape,lse.n20l1s3$rmse.shape, lse.n20l1s3$bias.scale, lse.n20l1s3$rmse.scale, lse.n20l1s3$Dobs, lse.n20l1s3$Dmax),
  wls = c(wlse.n20l1s3$bias.shape,wlse.n20l1s3$rmse.shape, wlse.n20l1s3$bias.scale, wlse.n20l1s3$rmse.scale, wlse.n20l1s3$Dobs, wlse.n20l1s3$Dmax),
  pce = c(pce.n20l1s3$bias.shape,pce.n20l1s3$rmse.shape, pce.n20l1s3$bias.scale, pce.n20l1s3$rmse.scale, pce.n20l1s3$Dobs, pce.n20l1s3$Dmax),
  mps = c(mps.n20l1s3$bias.shape,mps.n20l1s3$rmse.shape, mps.n20l1s3$bias.scale, mps.n20l1s3$rmse.scale, mps.n20l1s3$Dobs, mps.n20l1s3$Dmax),
  cvm = c(cvm.n20l1s3$bias.shape,cvm.n20l1s3$rmse.shape, cvm.n20l1s3$bias.scale, cvm.n20l1s3$rmse.scale, cvm.n20l1s3$Dobs, cvm.n20l1s3$Dmax),
  ad = c(ad.n20l1s3$bias.shape,ad.n20l1s3$rmse.shape, ad.n20l1s3$bias.scale, ad.n20l1s3$rmse.scale, ad.n20l1s3$Dobs, ad.n20l1s3$Dmax),
  rad = c(rtad.n20l1s3$bias.shape,rtad.n20l1s3$rmse.shape, rtad.n20l1s3$bias.scale, rtad.n20l1s3$rmse.scale, rtad.n20l1s3$Dobs, rtad.n20l1s3$Dmax)
)

n20 <- n20l1s3
round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))


# Combining the estimates n = 50
n50l1s3 <- data.frame(
  mle = c(mle.n50l1s3$bias.shape,mle.n50l1s3$rmse.shape, mle.n50l1s3$bias.scale, mle.n50l1s3$rmse.scale, mle.n50l1s3$Dobs, mle.n50l1s3$Dmax),
  mme = c(mme.n50l1s3$bias.shape,mme.n50l1s3$rmse.shape, mme.n50l1s3$bias.scale, mme.n50l1s3$rmse.scale, mme.n50l1s3$Dobs, mme.n50l1s3$Dmax),
  lse = c(lse.n50l1s3$bias.shape,lse.n50l1s3$rmse.shape, lse.n50l1s3$bias.scale, lse.n50l1s3$rmse.scale, lse.n50l1s3$Dobs, lse.n50l1s3$Dmax),
  wls = c(wlse.n50l1s3$bias.shape,wlse.n50l1s3$rmse.shape, wlse.n50l1s3$bias.scale, wlse.n50l1s3$rmse.scale, wlse.n50l1s3$Dobs, wlse.n50l1s3$Dmax),
  pce = c(pce.n50l1s3$bias.shape,pce.n50l1s3$rmse.shape, pce.n50l1s3$bias.scale, pce.n50l1s3$rmse.scale, pce.n50l1s3$Dobs, pce.n50l1s3$Dmax),
  mps = c(mps.n50l1s3$bias.shape,mps.n50l1s3$rmse.shape, mps.n50l1s3$bias.scale, mps.n50l1s3$rmse.scale, mps.n50l1s3$Dobs, mps.n50l1s3$Dmax),
  cvm = c(cvm.n50l1s3$bias.shape,cvm.n50l1s3$rmse.shape, cvm.n50l1s3$bias.scale, cvm.n50l1s3$rmse.scale, cvm.n50l1s3$Dobs, cvm.n50l1s3$Dmax),
  ad = c(ad.n50l1s3$bias.shape,ad.n50l1s3$rmse.shape, ad.n50l1s3$bias.scale, ad.n50l1s3$rmse.scale, ad.n50l1s3$Dobs, ad.n50l1s3$Dmax),
  rad = c(rtad.n50l1s3$bias.shape,rtad.n50l1s3$rmse.shape, rtad.n50l1s3$bias.scale, rtad.n50l1s3$rmse.scale, rtad.n50l1s3$Dobs, rtad.n50l1s3$Dmax)
)

n50 <- n50l1s3
round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))


# Combining the estimates n = 100
n100l1s3 <- data.frame(
  mle = c(mle.n100l1s3$bias.shape,mle.n100l1s3$rmse.shape, mle.n100l1s3$bias.scale, mle.n100l1s3$rmse.scale, mle.n100l1s3$Dobs, mle.n100l1s3$Dmax),
  mme = c(mme.n100l1s3$bias.shape,mme.n100l1s3$rmse.shape, mme.n100l1s3$bias.scale, mme.n100l1s3$rmse.scale, mme.n100l1s3$Dobs, mme.n100l1s3$Dmax),
  lse = c(lse.n100l1s3$bias.shape,lse.n100l1s3$rmse.shape, lse.n100l1s3$bias.scale, lse.n100l1s3$rmse.scale, lse.n100l1s3$Dobs, lse.n100l1s3$Dmax),
  wls = c(wlse.n100l1s3$bias.shape,wlse.n100l1s3$rmse.shape, wlse.n100l1s3$bias.scale, wlse.n100l1s3$rmse.scale, wlse.n100l1s3$Dobs, wlse.n100l1s3$Dmax),
  pce = c(pce.n100l1s3$bias.shape,pce.n100l1s3$rmse.shape, pce.n100l1s3$bias.scale, pce.n100l1s3$rmse.scale, pce.n100l1s3$Dobs, pce.n100l1s3$Dmax),
  mps = c(mps.n100l1s3$bias.shape,mps.n100l1s3$rmse.shape, mps.n100l1s3$bias.scale, mps.n100l1s3$rmse.scale, mps.n100l1s3$Dobs, mps.n100l1s3$Dmax),
  cvm = c(cvm.n100l1s3$bias.shape,cvm.n100l1s3$rmse.shape, cvm.n100l1s3$bias.scale, cvm.n100l1s3$rmse.scale, cvm.n100l1s3$Dobs, cvm.n100l1s3$Dmax),
  ad = c(ad.n100l1s3$bias.shape,ad.n100l1s3$rmse.shape, ad.n100l1s3$bias.scale, ad.n100l1s3$rmse.scale, ad.n100l1s3$Dobs, ad.n100l1s3$Dmax),
  rad = c(rtad.n100l1s3$bias.shape,rtad.n100l1s3$rmse.shape, rtad.n100l1s3$bias.scale, rtad.n100l1s3$rmse.scale, rtad.n100l1s3$Dobs, rtad.n100l1s3$Dmax)
)

n100 <- n100l1s3
round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))

# Combining the estimates n = 200
n200l1s3 <- data.frame(
  mle = c(mle.n200l1s3$bias.shape,mle.n200l1s3$rmse.shape, mle.n200l1s3$bias.scale, mle.n200l1s3$rmse.scale, mle.n200l1s3$Dobs, mle.n200l1s3$Dmax),
  mme = c(mme.n200l1s3$bias.shape,mme.n200l1s3$rmse.shape, mme.n200l1s3$bias.scale, mme.n200l1s3$rmse.scale, mme.n200l1s3$Dobs, mme.n200l1s3$Dmax),
  lse = c(lse.n200l1s3$bias.shape,lse.n200l1s3$rmse.shape, lse.n200l1s3$bias.scale, lse.n200l1s3$rmse.scale, lse.n200l1s3$Dobs, lse.n200l1s3$Dmax),
  wls = c(wlse.n200l1s3$bias.shape,wlse.n200l1s3$rmse.shape, wlse.n200l1s3$bias.scale, wlse.n200l1s3$rmse.scale, wlse.n200l1s3$Dobs, wlse.n200l1s3$Dmax),
  pce = c(pce.n200l1s3$bias.shape,pce.n200l1s3$rmse.shape, pce.n200l1s3$bias.scale, pce.n200l1s3$rmse.scale, pce.n200l1s3$Dobs, pce.n200l1s3$Dmax),
  mps = c(mps.n200l1s3$bias.shape,mps.n200l1s3$rmse.shape, mps.n200l1s3$bias.scale, mps.n200l1s3$rmse.scale, mps.n200l1s3$Dobs, mps.n200l1s3$Dmax),
  cvm = c(cvm.n200l1s3$bias.shape,cvm.n200l1s3$rmse.shape, cvm.n200l1s3$bias.scale, cvm.n200l1s3$rmse.scale, cvm.n200l1s3$Dobs, cvm.n200l1s3$Dmax),
  ad = c(ad.n200l1s3$bias.shape,ad.n200l1s3$rmse.shape, ad.n200l1s3$bias.scale, ad.n200l1s3$rmse.scale, ad.n200l1s3$Dobs, ad.n200l1s3$Dmax),
  rad = c(rtad.n200l1s3$bias.shape,rtad.n200l1s3$rmse.shape, rtad.n200l1s3$bias.scale, rtad.n200l1s3$rmse.scale, rtad.n200l1s3$Dobs, rtad.n200l1s3$Dmax)
)

n200 <- n200l1s3
round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))

sink()


# End#

##################################################
####### Additional Simulation ####################

# (lambda = 50, sigma = 20) DONE
# (lamdda = 100, sigma = 20) DONE
# (lambda = 50, sigma = 50)
# (lambda = 100, sigma = 50)

#################################################
############ RUNNING SIMULATION #################
############ lambda = 50, sigma=20 ##############
#################################################

ini = c(40, 15)
par=c(50, 20)

# MLE (Done for EG)
mle.n20l50s20 <-run.mle(n = 20, R= 5000, param=par, start=ini)
mle.n50l50s20  <- run.mle(n = 50, R= 5000, param=par, start=ini)
mle.n100l50s20  <- run.mle(n = 100, R= 5000, param=par, start=ini)
mle.n200l50s20  <- run.mle(n = 200, R= 5000, param=par, start=ini)


# MME 
mme.n20l50s20<-run.mme(n=20,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n50l50s20<-run.mme(n=50,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n100l50s20<-run.mme(n=100,R=5000, param=par, start=ini,low=1, up=Inf)
mme.n200l50s20<-run.mme(n=200,R=5000, param=par, start=ini,low=1, up=Inf)

# LME (Does not work for EG)
# lme.n20l1s3  <- run.lme(n = 20,R=5000, param=par, start=ini,low=1, up=Inf)
# lme.n50l1s3  <- run.lme(n = 50, R= 5000, param=par, start=ini)
# lme.n100l1s3  <- run.lme(n = 100, R= 5000, param=par, start=ini)
# lme.n200l1s3  <- run.lme(n = 200, R= 5000, param=par, start=ini)

# LSE 
lse.n20l50s20  <- run.lse(n = 20, R= 5000, param=par, start=ini)
lse.n50l50s20  <- run.lse(n = 50, R= 5000, param=par, start=ini)
lse.n100l50s20  <- run.lse(n = 100, R= 5000, param=par, start=ini)
lse.n200l50s20  <- run.lse(n = 200, R= 5000, param=par, start=ini)

# WLSE
wlse.n20l50s20  <- run.wlse(n = 20, R= 5000, param=par, start=ini)
wlse.n50l50s20  <- run.wlse(n = 50, R= 5000, param=par, start=ini)
wlse.n100l50s20  <- run.wlse(n = 100, R= 5000, param=par, start=ini)
wlse.n200l50s20  <- run.wlse(n = 200, R= 5000, param=par, start=ini)

# PCE
pce.n20l50s20  <- run.pce(n = 20, R= 5000, param=par, start=ini)
pce.n50l50s20  <- run.pce(n = 50, R= 5000, param=par, start= ini)
pce.n100l50s20  <- run.pce(n = 100, R= 5000, param=par, start=ini)
pce.n200l50s20  <- run.pce(n = 200, R= 5000, param=par, start=ini)

# MPS
mps.n20l50s20 <- run.mps(n = 20, R= 5000, param=par, start=ini)
mps.n50l50s20 <- run.mps(n = 50, R= 5000, param=par, start=ini)
mps.n100l50s20 <- run.mps(n = 100, R= 5000, param=par, start=ini)
mps.n200l50s20 <- run.mps(n = 200, R= 5000, param=par, start=ini)

# CVM
cvm.n20l50s20  <- run.cvm(n = 20, R= 5000, param=par, start=ini)
cvm.n50l50s20 <- run.cvm(n = 50, R= 5000, param=par, start=ini)
cvm.n100l50s20 <- run.cvm(n = 100, R= 5000, param=par, start=ini)
cvm.n200l50s20 <- run.cvm(n = 200, R= 5000, param=par, start=ini)

# AD
ad.n20l50s20  <- run.ad(n = 20, R= 5000, param=par, start=ini)
ad.n50l50s20  <- run.ad(n = 50, R= 5000, param=par, start=ini)
ad.n100l50s20 <- run.ad(n = 100, R= 5000, param=par, start=c(45, 18))
ad.n200l50s20 <- run.ad(n = 200, R= 5000, param=par, start=c(45, 18))

# RTAD
rtad.n20l50s20 <- run.rtad(n = 20, R= 5000, param=par, start = ini)
rtad.n50l50s20  <- run.rtad(n = 50, R= 5000, param=par, start = ini)
rtad.n100l50s20 <- run.rtad(n = 100, R= 5000, param=par, start = ini)
rtad.n200l50s20 <- run.rtad(n = 200, R= 5000, param=par, start = ini)


save.image("Dec 23 2015.RData")

sink("output_l50s20.txt", type = c("output", "message"))

# Combining the estimates n = 20
n20l50s20 <- data.frame(
  mle = c(mle.n20l50s20$bias.shape,mle.n20l50s20$rmse.shape, mle.n20l50s20$bias.scale, mle.n20l50s20$rmse.scale, mle.n20l50s20$Dobs, mle.n20l50s20$Dmax),
  mme = c(mme.n20l50s20$bias.shape,mme.n20l50s20$rmse.shape, mme.n20l50s20$bias.scale, mme.n20l50s20$rmse.scale, mme.n20l50s20$Dobs, mme.n20l50s20$Dmax),
  lse = c(lse.n20l50s20$bias.shape,lse.n20l50s20$rmse.shape, lse.n20l50s20$bias.scale, lse.n20l50s20$rmse.scale, lse.n20l50s20$Dobs, lse.n20l50s20$Dmax),
  wls = c(wlse.n20l50s20$bias.shape,wlse.n20l50s20$rmse.shape, wlse.n20l50s20$bias.scale, wlse.n20l50s20$rmse.scale, wlse.n20l50s20$Dobs, wlse.n20l50s20$Dmax),
  pce = c(pce.n20l50s20$bias.shape,pce.n20l50s20$rmse.shape, pce.n20l50s20$bias.scale, pce.n20l50s20$rmse.scale, pce.n20l50s20$Dobs, pce.n20l50s20$Dmax),
  mps = c(mps.n20l50s20$bias.shape,mps.n20l50s20$rmse.shape, mps.n20l50s20$bias.scale, mps.n20l50s20$rmse.scale, mps.n20l50s20$Dobs, mps.n20l50s20$Dmax),
  cvm = c(cvm.n20l50s20$bias.shape,cvm.n20l50s20$rmse.shape, cvm.n20l50s20$bias.scale, cvm.n20l50s20$rmse.scale, cvm.n20l50s20$Dobs, cvm.n20l50s20$Dmax),
  ad = c(ad.n20l50s20$bias.shape,ad.n20l50s20$rmse.shape, ad.n20l50s20$bias.scale, ad.n20l50s20$rmse.scale, ad.n20l50s20$Dobs, ad.n20l50s20$Dmax),
  rad = c(rtad.n20l50s20$bias.shape,rtad.n20l50s20$rmse.shape, rtad.n20l50s20$bias.scale, rtad.n20l50s20$rmse.scale, rtad.n20l50s20$Dobs, rtad.n20l50s20$Dmax)
)

n20 <- n20l50s20
round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))


# Combining the estimates n = 50
n50l50s20 <- data.frame(
  mle = c(mle.n50l50s20$bias.shape,mle.n50l50s20$rmse.shape, mle.n50l50s20$bias.scale, mle.n50l50s20$rmse.scale, mle.n50l50s20$Dobs, mle.n50l50s20$Dmax),
  mme = c(mme.n50l50s20$bias.shape,mme.n50l50s20$rmse.shape, mme.n50l50s20$bias.scale, mme.n50l50s20$rmse.scale, mme.n50l50s20$Dobs, mme.n50l50s20$Dmax),
  lse = c(lse.n50l50s20$bias.shape,lse.n50l50s20$rmse.shape, lse.n50l50s20$bias.scale, lse.n50l50s20$rmse.scale, lse.n50l50s20$Dobs, lse.n50l50s20$Dmax),
  wls = c(wlse.n50l50s20$bias.shape,wlse.n50l50s20$rmse.shape, wlse.n50l50s20$bias.scale, wlse.n50l50s20$rmse.scale, wlse.n50l50s20$Dobs, wlse.n50l50s20$Dmax),
  pce = c(pce.n50l50s20$bias.shape,pce.n50l50s20$rmse.shape, pce.n50l50s20$bias.scale, pce.n50l50s20$rmse.scale, pce.n50l50s20$Dobs, pce.n50l50s20$Dmax),
  mps = c(mps.n50l50s20$bias.shape,mps.n50l50s20$rmse.shape, mps.n50l50s20$bias.scale, mps.n50l50s20$rmse.scale, mps.n50l50s20$Dobs, mps.n50l50s20$Dmax),
  cvm = c(cvm.n50l50s20$bias.shape,cvm.n50l50s20$rmse.shape, cvm.n50l50s20$bias.scale, cvm.n50l50s20$rmse.scale, cvm.n50l50s20$Dobs, cvm.n50l50s20$Dmax),
  ad = c(ad.n50l50s20$bias.shape,ad.n50l50s20$rmse.shape, ad.n50l50s20$bias.scale, ad.n50l50s20$rmse.scale, ad.n50l50s20$Dobs, ad.n50l50s20$Dmax),
  rad = c(rtad.n50l50s20$bias.shape,rtad.n50l50s20$rmse.shape, rtad.n50l50s20$bias.scale, rtad.n50l50s20$rmse.scale, rtad.n50l50s20$Dobs, rtad.n50l50s20$Dmax)
)

n50 <- n50l50s20
round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))


# Combining the estimates n = 100
n100l50s20 <- data.frame(
  mle = c(mle.n100l50s20$bias.shape,mle.n100l50s20$rmse.shape, mle.n100l50s20$bias.scale, mle.n100l50s20$rmse.scale, mle.n100l50s20$Dobs, mle.n100l50s20$Dmax),
  mme = c(mme.n100l50s20$bias.shape,mme.n100l50s20$rmse.shape, mme.n100l50s20$bias.scale, mme.n100l50s20$rmse.scale, mme.n100l50s20$Dobs, mme.n100l50s20$Dmax),
  lse = c(lse.n100l50s20$bias.shape,lse.n100l50s20$rmse.shape, lse.n100l50s20$bias.scale, lse.n100l50s20$rmse.scale, lse.n100l50s20$Dobs, lse.n100l50s20$Dmax),
  wls = c(wlse.n100l50s20$bias.shape,wlse.n100l50s20$rmse.shape, wlse.n100l50s20$bias.scale, wlse.n100l50s20$rmse.scale, wlse.n100l50s20$Dobs, wlse.n100l50s20$Dmax),
  pce = c(pce.n100l50s20$bias.shape,pce.n100l50s20$rmse.shape, pce.n100l50s20$bias.scale, pce.n100l50s20$rmse.scale, pce.n100l50s20$Dobs, pce.n100l50s20$Dmax),
  mps = c(mps.n100l50s20$bias.shape,mps.n100l50s20$rmse.shape, mps.n100l50s20$bias.scale, mps.n100l50s20$rmse.scale, mps.n100l50s20$Dobs, mps.n100l50s20$Dmax),
  cvm = c(cvm.n100l50s20$bias.shape,cvm.n100l50s20$rmse.shape, cvm.n100l50s20$bias.scale, cvm.n100l50s20$rmse.scale, cvm.n100l50s20$Dobs, cvm.n100l50s20$Dmax),
  ad = c(ad.n100l50s20$bias.shape,ad.n100l50s20$rmse.shape, ad.n100l50s20$bias.scale, ad.n100l50s20$rmse.scale, ad.n100l50s20$Dobs, ad.n100l50s20$Dmax),
  rad = c(rtad.n100l50s20$bias.shape,rtad.n100l50s20$rmse.shape, rtad.n100l50s20$bias.scale, rtad.n100l50s20$rmse.scale, rtad.n100l50s20$Dobs, rtad.n100l50s20$Dmax)
)

n100 <- n100l50s20
round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))

# Combining the estimates n = 200
n200l50s20 <- data.frame(
  mle = c(mle.n200l50s20$bias.shape,mle.n200l50s20$rmse.shape, mle.n200l50s20$bias.scale, mle.n200l50s20$rmse.scale, mle.n200l50s20$Dobs, mle.n200l50s20$Dmax),
  mme = c(mme.n200l50s20$bias.shape,mme.n200l50s20$rmse.shape, mme.n200l50s20$bias.scale, mme.n200l50s20$rmse.scale, mme.n200l50s20$Dobs, mme.n200l50s20$Dmax),
  lse = c(lse.n200l50s20$bias.shape,lse.n200l50s20$rmse.shape, lse.n200l50s20$bias.scale, lse.n200l50s20$rmse.scale, lse.n200l50s20$Dobs, lse.n200l50s20$Dmax),
  wls = c(wlse.n200l50s20$bias.shape,wlse.n200l50s20$rmse.shape, wlse.n200l50s20$bias.scale, wlse.n200l50s20$rmse.scale, wlse.n200l50s20$Dobs, wlse.n200l50s20$Dmax),
  pce = c(pce.n200l50s20$bias.shape,pce.n200l50s20$rmse.shape, pce.n200l50s20$bias.scale, pce.n200l50s20$rmse.scale, pce.n200l50s20$Dobs, pce.n200l50s20$Dmax),
  mps = c(mps.n200l50s20$bias.shape,mps.n200l50s20$rmse.shape, mps.n200l50s20$bias.scale, mps.n200l50s20$rmse.scale, mps.n200l50s20$Dobs, mps.n200l50s20$Dmax),
  cvm = c(cvm.n200l50s20$bias.shape,cvm.n200l50s20$rmse.shape, cvm.n200l50s20$bias.scale, cvm.n200l50s20$rmse.scale, cvm.n200l50s20$Dobs, cvm.n200l50s20$Dmax),
  ad = c(ad.n200l50s20$bias.shape,ad.n200l50s20$rmse.shape, ad.n200l50s20$bias.scale, ad.n200l50s20$rmse.scale, ad.n200l50s20$Dobs, ad.n200l50s20$Dmax),
  rad = c(rtad.n200l50s20$bias.shape,rtad.n200l50s20$rmse.shape, rtad.n200l50s20$bias.scale, rtad.n200l50s20$rmse.scale, rtad.n200l50s20$Dobs, rtad.n200l50s20$Dmax)
)

n200 <- n200l50s20
round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))

sink()


# End#



#################################################
############ RUNNING SIMULATION TO RUN #################
############ lambda = 100, sigma=20 ##############
#################################################

ini = c(80, 15)
par=c(100, 20)

# MLE (Done for EG)
mle.n20l100s20 <-run.mle(n = 20, R= 5000, param=par, start=ini)
mle.n50l100s20  <- run.mle(n = 50, R= 5000, param=par, start=ini)
mle.n100l100s20  <- run.mle(n = 100, R= 5000, param=par, start=ini)
mle.n200l100s20  <- run.mle(n = 200, R= 5000, param=par, start=ini)


# MME 
mme.n20l100s20<-run.mme(n=20,R=5000, param=par, start=ini,low=1, up=500)
mme.n50l100s20<-run.mme(n=50,R=5000, param=par, start=ini,low=1, up=500)
mme.n100l100s20<-run.mme(n=100,R=5000, param=par, start=ini,low=1, up=500)
mme.n200l100s20<-run.mme(n=200,R=5000, param=par, start=ini,low=1, up=500)

# LME (Does not work for EG)
# lme.n20l1s3  <- run.lme(n = 20,R=5000, param=par, start=ini,low=1, up=Inf)
# lme.n50l1s3  <- run.lme(n = 50, R= 5000, param=par, start=ini)
# lme.n100l1s3  <- run.lme(n = 100, R= 5000, param=par, start=ini)
# lme.n200l1s3  <- run.lme(n = 200, R= 5000, param=par, start=ini)

# LSE 
lse.n20l100s20  <- run.lse(n = 20, R= 5000, param=par, start=ini)
lse.n50l100s20  <- run.lse(n = 50, R= 5000, param=par, start=ini)
lse.n100l100s20  <- run.lse(n = 100, R= 5000, param=par, start=ini)
lse.n200l100s20  <- run.lse(n = 200, R= 5000, param=par, start=ini)

# WLSE
wlse.n20l100s20  <- run.wlse(n = 20, R= 5000, param=par, start=ini)
wlse.n50l100s20  <- run.wlse(n = 50, R= 5000, param=par, start=ini)
wlse.n100l100s20  <- run.wlse(n = 100, R= 5000, param=par, start=ini)
wlse.n200l100s20  <- run.wlse(n = 200, R= 5000, param=par, start=ini)

# PCE
pce.n20l100s20  <- run.pce(n = 20, R= 5000, param=par, start=ini)
pce.n50l100s20  <- run.pce(n = 50, R= 5000, param=par, start= ini)
pce.n100l100s20  <- run.pce(n = 100, R= 5000, param=par, start=ini)
pce.n200l100s20  <- run.pce(n = 200, R= 5000, param=par, start=ini)

# MPS
mps.n20l100s20 <- run.mps(n = 20, R= 5000, param=par, start=ini)
mps.n50l100s20 <- run.mps(n = 50, R= 5000, param=par, start=ini)
mps.n100l100s20 <- run.mps(n = 100, R= 5000, param=par, start=ini)
mps.n200l100s20 <- run.mps(n = 200, R= 5000, param=par, start=ini)

# CVM
cvm.n20l100s20  <- run.cvm(n = 20, R= 5000, param=par, start=ini)
cvm.n50l100s20 <- run.cvm(n = 50, R= 5000, param=par, start=ini)
cvm.n100l100s20 <- run.cvm(n = 100, R= 5000, param=par, start=ini)
cvm.n200l100s20 <- run.cvm(n = 200, R= 5000, param=par, start=ini)

# AD
ini.ad <- c(90, 17)
ad.n20l100s20  <- run.ad(n = 20, R= 5000, param=par, start=ini.ad
ad.n50l100s20  <- run.ad(n = 50, R= 5000, param=par, start=ini.ad)
ad.n100l100s20 <- run.ad(n = 100, R= 5000, param=par, start=ini.ad)
ad.n200l100s20 <- run.ad(n = 200, R= 5000, param=par, start=ini.ad)

# RTAD
rtad.n20l100s20 <- run.rtad(n = 20, R= 5000, param=par, start = ini)
rtad.n50l100s20  <- run.rtad(n = 50, R= 5000, param=par, start = ini)
rtad.n100l100s20 <- run.rtad(n = 100, R= 5000, param=par, start = ini)
rtad.n200l100s20 <- run.rtad(n = 200, R= 5000, param=par, start = ini)


save.image("Dec 23 2015_2.RData")

sink("output_l100s20.txt", type = c("output", "message"))

# Combining the estimates n = 20
n20l100s20 <- data.frame(
  mle = c(mle.n20l100s20$bias.shape,mle.n20l100s20$rmse.shape, mle.n20l100s20$bias.scale, mle.n20l100s20$rmse.scale, mle.n20l100s20$Dobs, mle.n20l100s20$Dmax),
  mme = c(mme.n20l100s20$bias.shape,mme.n20l100s20$rmse.shape, mme.n20l100s20$bias.scale, mme.n20l100s20$rmse.scale, mme.n20l100s20$Dobs, mme.n20l100s20$Dmax),
  lse = c(lse.n20l100s20$bias.shape,lse.n20l100s20$rmse.shape, lse.n20l100s20$bias.scale, lse.n20l100s20$rmse.scale, lse.n20l100s20$Dobs, lse.n20l100s20$Dmax),
  wls = c(wlse.n20l100s20$bias.shape,wlse.n20l100s20$rmse.shape, wlse.n20l100s20$bias.scale, wlse.n20l100s20$rmse.scale, wlse.n20l100s20$Dobs, wlse.n20l100s20$Dmax),
  pce = c(pce.n20l100s20$bias.shape,pce.n20l100s20$rmse.shape, pce.n20l100s20$bias.scale, pce.n20l100s20$rmse.scale, pce.n20l100s20$Dobs, pce.n20l100s20$Dmax),
  mps = c(mps.n20l100s20$bias.shape,mps.n20l100s20$rmse.shape, mps.n20l100s20$bias.scale, mps.n20l100s20$rmse.scale, mps.n20l100s20$Dobs, mps.n20l100s20$Dmax),
  cvm = c(cvm.n20l100s20$bias.shape,cvm.n20l100s20$rmse.shape, cvm.n20l100s20$bias.scale, cvm.n20l100s20$rmse.scale, cvm.n20l100s20$Dobs, cvm.n20l100s20$Dmax),
  ad = c(ad.n20l100s20$bias.shape,ad.n20l100s20$rmse.shape, ad.n20l100s20$bias.scale, ad.n20l100s20$rmse.scale, ad.n20l100s20$Dobs, ad.n20l100s20$Dmax),
  rad = c(rtad.n20l100s20$bias.shape,rtad.n20l100s20$rmse.shape, rtad.n20l100s20$bias.scale, rtad.n20l100s20$rmse.scale, rtad.n20l100s20$Dobs, rtad.n20l100s20$Dmax)
)

n20 <- n20l100s20
round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))


# Combining the estimates n = 50
n50l100s20 <- data.frame(
  mle = c(mle.n50l100s20$bias.shape,mle.n50l100s20$rmse.shape, mle.n50l100s20$bias.scale, mle.n50l100s20$rmse.scale, mle.n50l100s20$Dobs, mle.n50l100s20$Dmax),
  mme = c(mme.n50l100s20$bias.shape,mme.n50l100s20$rmse.shape, mme.n50l100s20$bias.scale, mme.n50l100s20$rmse.scale, mme.n50l100s20$Dobs, mme.n50l100s20$Dmax),
  lse = c(lse.n50l100s20$bias.shape,lse.n50l100s20$rmse.shape, lse.n50l100s20$bias.scale, lse.n50l100s20$rmse.scale, lse.n50l100s20$Dobs, lse.n50l100s20$Dmax),
  wls = c(wlse.n50l100s20$bias.shape,wlse.n50l100s20$rmse.shape, wlse.n50l100s20$bias.scale, wlse.n50l100s20$rmse.scale, wlse.n50l100s20$Dobs, wlse.n50l100s20$Dmax),
  pce = c(pce.n50l100s20$bias.shape,pce.n50l100s20$rmse.shape, pce.n50l100s20$bias.scale, pce.n50l100s20$rmse.scale, pce.n50l100s20$Dobs, pce.n50l100s20$Dmax),
  mps = c(mps.n50l100s20$bias.shape,mps.n50l100s20$rmse.shape, mps.n50l100s20$bias.scale, mps.n50l100s20$rmse.scale, mps.n50l100s20$Dobs, mps.n50l100s20$Dmax),
  cvm = c(cvm.n50l100s20$bias.shape,cvm.n50l100s20$rmse.shape, cvm.n50l100s20$bias.scale, cvm.n50l100s20$rmse.scale, cvm.n50l100s20$Dobs, cvm.n50l100s20$Dmax),
  ad = c(ad.n50l100s20$bias.shape,ad.n50l100s20$rmse.shape, ad.n50l100s20$bias.scale, ad.n50l100s20$rmse.scale, ad.n50l100s20$Dobs, ad.n50l100s20$Dmax),
  rad = c(rtad.n50l100s20$bias.shape,rtad.n50l100s20$rmse.shape, rtad.n50l100s20$bias.scale, rtad.n50l100s20$rmse.scale, rtad.n50l100s20$Dobs, rtad.n50l100s20$Dmax)
)

n50 <- n50l100s20
round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))


# Combining the estimates n = 100
n100l100s20 <- data.frame(
  mle = c(mle.n100l100s20$bias.shape,mle.n100l100s20$rmse.shape, mle.n100l100s20$bias.scale, mle.n100l100s20$rmse.scale, mle.n100l100s20$Dobs, mle.n100l100s20$Dmax),
  mme = c(mme.n100l100s20$bias.shape,mme.n100l100s20$rmse.shape, mme.n100l100s20$bias.scale, mme.n100l100s20$rmse.scale, mme.n100l100s20$Dobs, mme.n100l100s20$Dmax),
  lse = c(lse.n100l100s20$bias.shape,lse.n100l100s20$rmse.shape, lse.n100l100s20$bias.scale, lse.n100l100s20$rmse.scale, lse.n100l100s20$Dobs, lse.n100l100s20$Dmax),
  wls = c(wlse.n100l100s20$bias.shape,wlse.n100l100s20$rmse.shape, wlse.n100l100s20$bias.scale, wlse.n100l100s20$rmse.scale, wlse.n100l100s20$Dobs, wlse.n100l100s20$Dmax),
  pce = c(pce.n100l100s20$bias.shape,pce.n100l100s20$rmse.shape, pce.n100l100s20$bias.scale, pce.n100l100s20$rmse.scale, pce.n100l100s20$Dobs, pce.n100l100s20$Dmax),
  mps = c(mps.n100l100s20$bias.shape,mps.n100l100s20$rmse.shape, mps.n100l100s20$bias.scale, mps.n100l100s20$rmse.scale, mps.n100l100s20$Dobs, mps.n100l100s20$Dmax),
  cvm = c(cvm.n100l100s20$bias.shape,cvm.n100l100s20$rmse.shape, cvm.n100l100s20$bias.scale, cvm.n100l100s20$rmse.scale, cvm.n100l100s20$Dobs, cvm.n100l100s20$Dmax),
  ad = c(ad.n100l100s20$bias.shape,ad.n100l100s20$rmse.shape, ad.n100l100s20$bias.scale, ad.n100l100s20$rmse.scale, ad.n100l100s20$Dobs, ad.n100l100s20$Dmax),
  rad = c(rtad.n100l100s20$bias.shape,rtad.n100l100s20$rmse.shape, rtad.n100l100s20$bias.scale, rtad.n100l100s20$rmse.scale, rtad.n100l100s20$Dobs, rtad.n100l100s20$Dmax)
)

n100 <- n100l100s20
round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))

# Combining the estimates n = 200
n200l100s20 <- data.frame(
  mle = c(mle.n200l100s20$bias.shape,mle.n200l100s20$rmse.shape, mle.n200l100s20$bias.scale, mle.n200l100s20$rmse.scale, mle.n200l100s20$Dobs, mle.n200l100s20$Dmax),
  mme = c(mme.n200l100s20$bias.shape,mme.n200l100s20$rmse.shape, mme.n200l100s20$bias.scale, mme.n200l100s20$rmse.scale, mme.n200l100s20$Dobs, mme.n200l100s20$Dmax),
  lse = c(lse.n200l100s20$bias.shape,lse.n200l100s20$rmse.shape, lse.n200l100s20$bias.scale, lse.n200l100s20$rmse.scale, lse.n200l100s20$Dobs, lse.n200l100s20$Dmax),
  wls = c(wlse.n200l100s20$bias.shape,wlse.n200l100s20$rmse.shape, wlse.n200l100s20$bias.scale, wlse.n200l100s20$rmse.scale, wlse.n200l100s20$Dobs, wlse.n200l100s20$Dmax),
  pce = c(pce.n200l100s20$bias.shape,pce.n200l100s20$rmse.shape, pce.n200l100s20$bias.scale, pce.n200l100s20$rmse.scale, pce.n200l100s20$Dobs, pce.n200l100s20$Dmax),
  mps = c(mps.n200l100s20$bias.shape,mps.n200l100s20$rmse.shape, mps.n200l100s20$bias.scale, mps.n200l100s20$rmse.scale, mps.n200l100s20$Dobs, mps.n200l100s20$Dmax),
  cvm = c(cvm.n200l100s20$bias.shape,cvm.n200l100s20$rmse.shape, cvm.n200l100s20$bias.scale, cvm.n200l100s20$rmse.scale, cvm.n200l100s20$Dobs, cvm.n200l100s20$Dmax),
  ad = c(ad.n200l100s20$bias.shape,ad.n200l100s20$rmse.shape, ad.n200l100s20$bias.scale, ad.n200l100s20$rmse.scale, ad.n200l100s20$Dobs, ad.n200l100s20$Dmax),
  rad = c(rtad.n200l100s20$bias.shape,rtad.n200l100s20$rmse.shape, rtad.n200l100s20$bias.scale, rtad.n200l100s20$rmse.scale, rtad.n200l100s20$Dobs, rtad.n200l100s20$Dmax)
)

n200 <- n200l100s20
round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))

sink()


# End#

############ Simulation to confirm real data results ##########
# The real data has lambda = 1.65, sigma = 10

#################################################
############ RUNNING SIMULATION for real data ########
############ lambda = 100, sigma=20 ##############
#################################################


ini = c(1, 8)
par=c(1.65, 10)

# MLE (Done for EG)
mle.n20l1.6s10 <-run.mle(n = 20, R=1000, param=par, start=ini)
mle.n50l1.6s10  <- run.mle(n = 50, R= 1000, param=par, start=ini)
mle.n100l1.6s10  <- run.mle(n = 100, R= 1000, param=par, start=ini)
mle.n200l1.6s10  <- run.mle(n = 200, R= 1000, param=par, start=ini)


# MME 
mme.n20l1.6s10<-run.mme(n=20,R=1000, param=par, start=ini,low=1, up=500)
mme.n50l1.6s10<-run.mme(n=50,R=1000, param=par, start=ini,low=1, up=500)
mme.n100l1.6s10<-run.mme(n=100,R=1000, param=par, start=ini,low=1, up=500)
mme.n200l1.6s10<-run.mme(n=200,R=1000, param=par, start=ini,low=1, up=500)

# LME (Does not work for EG)
# lme.n20l1s3  <- run.lme(n = 20,R=5000, param=par, start=ini,low=1, up=Inf)
# lme.n50l1s3  <- run.lme(n = 50, R= 5000, param=par, start=ini)
# lme.n100l1s3  <- run.lme(n = 100, R= 5000, param=par, start=ini)
# lme.n200l1s3  <- run.lme(n = 200, R= 5000, param=par, start=ini)

# LSE 
lse.n20l1.6s10  <- run.lse(n = 20, R= 1000, param=par, start=ini)
lse.n50l1.6s10  <- run.lse(n = 50, R= 1000, param=par, start=ini)
lse.n100l1.6s10  <- run.lse(n = 100, R= 1000, param=par, start=ini)
lse.n200l1.6s10  <- run.lse(n = 200, R= 1000, param=par, start=ini)

# WLSE
wlse.n20l1.6s10  <- run.wlse(n = 20, R= 1000, param=par, start=ini)
wlse.n50l1.6s10  <- run.wlse(n = 50, R= 1000, param=par, start=ini)
wlse.n100l1.6s10  <- run.wlse(n = 100, R= 1000, param=par, start=ini)
wlse.n200l1.6s10  <- run.wlse(n = 200, R= 1000, param=par, start=ini)

# PCE
pce.n20l1.6s10  <- run.pce(n = 20, R= 1000, param=par, start=ini)
pce.n50l1.6s10  <- run.pce(n = 50, R= 1000, param=par, start= ini)
pce.n100l1.6s10  <- run.pce(n = 100, R= 1000, param=par, start=ini)
pce.n200l1.6s10  <- run.pce(n = 200, R= 1000, param=par, start=ini)

# MPS
mps.n20l1.6s10 <- run.mps(n = 20, R= 1000, param=par, start=ini)
mps.n50l1.6s10 <- run.mps(n = 50, R= 1000, param=par, start=ini)
mps.n100l1.6s10 <- run.mps(n = 100, R= 1000, param=par, start=ini)
mps.n200l1.6s10 <- run.mps(n = 200, R= 1000, param=par, start=ini)

# CVM
cvm.n20l1.6s10  <- run.cvm(n = 20, R=1000, param=par, start=ini)
cvm.n50l1.6s10 <- run.cvm(n = 50, R= 1000, param=par, start=ini)
cvm.n100l1.6s10 <- run.cvm(n = 100, R= 1000, param=par, start=ini)
cvm.n200l1.6s10 <- run.cvm(n = 200, R= 1000, param=par, start=ini)

# AD
ini.ad <- c(1, 8)
ad.n20l1.6s10  <- run.ad(n = 20, R= 1000, param=par, start=ini.ad)
ad.n50l1.6s10  <- run.ad(n = 50, R= 1000, param=par, start=ini.ad)
ad.n100l1.6s10 <- run.ad(n = 100, R= 1000, param=par, start=ini.ad)
ad.n200l1.6s10 <- run.ad(n = 200, R= 1000, param=par, start=ini.ad)
                         
# RTAD
rtad.n20l1.6s10 <- run.rtad(n = 20, R= 1000, param=par, start = ini)
rtad.n50l1.6s10  <- run.rtad(n = 50, R= 1000, param=par, start = ini)
rtad.n100l1.6s10 <- run.rtad(n = 100, R= 1000, param=par, start = ini)
rtad.n200l1.6s10 <- run.rtad(n = 200, R= 1000, param=par, start = ini)

save.image("Dec 24 2015_realData.RData")

sink("output_l1.6s10.txt", type = c("output", "message"))
                         
# Combining the estimates n = 20
n20l1.6s10 <- data.frame(
mle = c(mle.n20l1.6s10$bias.shape,mle.n20l1.6s10$rmse.shape, mle.n20l1.6s10$bias.scale, mle.n20l1.6s10$rmse.scale, mle.n20l1.6s10$Dobs, mle.n20l1.6s10$Dmax),
mme = c(mme.n20l1.6s10$bias.shape,mme.n20l1.6s10$rmse.shape, mme.n20l1.6s10$bias.scale, mme.n20l1.6s10$rmse.scale, mme.n20l1.6s10$Dobs, mme.n20l1.6s10$Dmax),
lse = c(lse.n20l1.6s10$bias.shape,lse.n20l1.6s10$rmse.shape, lse.n20l1.6s10$bias.scale, lse.n20l1.6s10$rmse.scale, lse.n20l1.6s10$Dobs, lse.n20l1.6s10$Dmax),
wls = c(wlse.n20l1.6s10$bias.shape,wlse.n20l1.6s10$rmse.shape, wlse.n20l1.6s10$bias.scale, wlse.n20l1.6s10$rmse.scale, wlse.n20l1.6s10$Dobs, wlse.n20l1.6s10$Dmax),
pce = c(pce.n20l1.6s10$bias.shape,pce.n20l1.6s10$rmse.shape, pce.n20l1.6s10$bias.scale, pce.n20l1.6s10$rmse.scale, pce.n20l1.6s10$Dobs, pce.n20l1.6s10$Dmax),
mps = c(mps.n20l1.6s10$bias.shape,mps.n20l1.6s10$rmse.shape, mps.n20l1.6s10$bias.scale, mps.n20l1.6s10$rmse.scale, mps.n20l1.6s10$Dobs, mps.n20l1.6s10$Dmax),
cvm = c(cvm.n20l1.6s10$bias.shape,cvm.n20l1.6s10$rmse.shape, cvm.n20l1.6s10$bias.scale, cvm.n20l1.6s10$rmse.scale, cvm.n20l1.6s10$Dobs, cvm.n20l1.6s10$Dmax),
ad = c(ad.n20l1.6s10$bias.shape,ad.n20l1.6s10$rmse.shape, ad.n20l1.6s10$bias.scale, ad.n20l1.6s10$rmse.scale, ad.n20l1.6s10$Dobs, ad.n20l1.6s10$Dmax),
rad = c(rtad.n20l1.6s10$bias.shape,rtad.n20l1.6s10$rmse.shape, rtad.n20l1.6s10$bias.scale, rtad.n20l1.6s10$rmse.scale, rtad.n20l1.6s10$Dobs, rtad.n20l1.6s10$Dmax)
                         )
                         
n20 <- n20l1.6s10
round(n20, 3)
t(apply(abs(n20), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n20), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n20), 1, rank)), 2, sum))
                         
                         
# Combining the estimates n = 50
n50l1.6s10 <- data.frame(
mle = c(mle.n50l1.6s10$bias.shape,mle.n50l1.6s10$rmse.shape, mle.n50l1.6s10$bias.scale, mle.n50l1.6s10$rmse.scale, mle.n50l1.6s10$Dobs, mle.n50l1.6s10$Dmax),
mme = c(mme.n50l1.6s10$bias.shape,mme.n50l1.6s10$rmse.shape, mme.n50l1.6s10$bias.scale, mme.n50l1.6s10$rmse.scale, mme.n50l1.6s10$Dobs, mme.n50l1.6s10$Dmax),
lse = c(lse.n50l1.6s10$bias.shape,lse.n50l1.6s10$rmse.shape, lse.n50l1.6s10$bias.scale, lse.n50l1.6s10$rmse.scale, lse.n50l1.6s10$Dobs, lse.n50l1.6s10$Dmax),
wls = c(wlse.n50l1.6s10$bias.shape,wlse.n50l1.6s10$rmse.shape, wlse.n50l1.6s10$bias.scale, wlse.n50l1.6s10$rmse.scale, wlse.n50l1.6s10$Dobs, wlse.n50l1.6s10$Dmax),
pce = c(pce.n50l1.6s10$bias.shape,pce.n50l1.6s10$rmse.shape, pce.n50l1.6s10$bias.scale, pce.n50l1.6s10$rmse.scale, pce.n50l1.6s10$Dobs, pce.n50l1.6s10$Dmax),
mps = c(mps.n50l1.6s10$bias.shape,mps.n50l1.6s10$rmse.shape, mps.n50l1.6s10$bias.scale, mps.n50l1.6s10$rmse.scale, mps.n50l1.6s10$Dobs, mps.n50l1.6s10$Dmax),
cvm = c(cvm.n50l1.6s10$bias.shape,cvm.n50l1.6s10$rmse.shape, cvm.n50l1.6s10$bias.scale, cvm.n50l1.6s10$rmse.scale, cvm.n50l1.6s10$Dobs, cvm.n50l1.6s10$Dmax),
ad = c(ad.n50l1.6s10$bias.shape,ad.n50l1.6s10$rmse.shape, ad.n50l1.6s10$bias.scale, ad.n50l1.6s10$rmse.scale, ad.n50l1.6s10$Dobs, ad.n50l1.6s10$Dmax),
rad = c(rtad.n50l1.6s10$bias.shape,rtad.n50l1.6s10$rmse.shape, rtad.n50l1.6s10$bias.scale, rtad.n50l1.6s10$rmse.scale, rtad.n50l1.6s10$Dobs, rtad.n50l1.6s10$Dmax)
)
                         
n50 <- n50l1.6s10
round(n50, 3)
t(apply(abs(n50), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n50), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n50), 1, rank)), 2, sum))

# Combining the estimates n = 100
n100l1.6s10 <- data.frame(
mle = c(mle.n100l1.6s10$bias.shape,mle.n100l1.6s10$rmse.shape, mle.n100l1.6s10$bias.scale, mle.n100l1.6s10$rmse.scale, mle.n100l1.6s10$Dobs, mle.n100l1.6s10$Dmax),
mme = c(mme.n100l1.6s10$bias.shape,mme.n100l1.6s10$rmse.shape, mme.n100l1.6s10$bias.scale, mme.n100l1.6s10$rmse.scale, mme.n100l1.6s10$Dobs, mme.n100l1.6s10$Dmax),
lse = c(lse.n100l1.6s10$bias.shape,lse.n100l1.6s10$rmse.shape, lse.n100l1.6s10$bias.scale, lse.n100l1.6s10$rmse.scale, lse.n100l1.6s10$Dobs, lse.n100l1.6s10$Dmax),
wls = c(wlse.n100l1.6s10$bias.shape,wlse.n100l1.6s10$rmse.shape, wlse.n100l1.6s10$bias.scale, wlse.n100l1.6s10$rmse.scale, wlse.n100l1.6s10$Dobs, wlse.n100l1.6s10$Dmax),
pce = c(pce.n100l1.6s10$bias.shape,pce.n100l1.6s10$rmse.shape, pce.n100l1.6s10$bias.scale, pce.n100l1.6s10$rmse.scale, pce.n100l1.6s10$Dobs, pce.n100l1.6s10$Dmax),
mps = c(mps.n100l1.6s10$bias.shape,mps.n100l1.6s10$rmse.shape, mps.n100l1.6s10$bias.scale, mps.n100l1.6s10$rmse.scale, mps.n100l1.6s10$Dobs, mps.n100l1.6s10$Dmax),
cvm = c(cvm.n100l1.6s10$bias.shape,cvm.n100l1.6s10$rmse.shape, cvm.n100l1.6s10$bias.scale, cvm.n100l1.6s10$rmse.scale, cvm.n100l1.6s10$Dobs, cvm.n100l1.6s10$Dmax),
ad = c(ad.n100l1.6s10$bias.shape,ad.n100l1.6s10$rmse.shape, ad.n100l1.6s10$bias.scale, ad.n100l1.6s10$rmse.scale, ad.n100l1.6s10$Dobs, ad.n100l1.6s10$Dmax),
rad = c(rtad.n100l1.6s10$bias.shape,rtad.n100l1.6s10$rmse.shape, rtad.n100l1.6s10$bias.scale, rtad.n100l1.6s10$rmse.scale, rtad.n100l1.6s10$Dobs, rtad.n100l1.6s10$Dmax)
)

n100 <- n100l1.6s10
round(n100, 3)
t(apply(abs(n100), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n100), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n100), 1, rank)), 2, sum))
                         
# Combining the estimates n = 200
n200l1.6s10 <- data.frame(
mle = c(mle.n200l1.6s10$bias.shape,mle.n200l1.6s10$rmse.shape, mle.n200l1.6s10$bias.scale, mle.n200l1.6s10$rmse.scale, mle.n200l1.6s10$Dobs, mle.n200l1.6s10$Dmax),
mme = c(mme.n200l1.6s10$bias.shape,mme.n200l1.6s10$rmse.shape, mme.n200l1.6s10$bias.scale, mme.n200l1.6s10$rmse.scale, mme.n200l1.6s10$Dobs, mme.n200l1.6s10$Dmax),
lse = c(lse.n200l1.6s10$bias.shape,lse.n200l1.6s10$rmse.shape, lse.n200l1.6s10$bias.scale, lse.n200l1.6s10$rmse.scale, lse.n200l1.6s10$Dobs, lse.n200l1.6s10$Dmax),
wls = c(wlse.n200l1.6s10$bias.shape,wlse.n200l1.6s10$rmse.shape, wlse.n200l1.6s10$bias.scale, wlse.n200l1.6s10$rmse.scale, wlse.n200l1.6s10$Dobs, wlse.n200l1.6s10$Dmax),
pce = c(pce.n200l1.6s10$bias.shape,pce.n200l1.6s10$rmse.shape, pce.n200l1.6s10$bias.scale, pce.n200l1.6s10$rmse.scale, pce.n200l1.6s10$Dobs, pce.n200l1.6s10$Dmax),
mps = c(mps.n200l1.6s10$bias.shape,mps.n200l1.6s10$rmse.shape, mps.n200l1.6s10$bias.scale, mps.n200l1.6s10$rmse.scale, mps.n200l1.6s10$Dobs, mps.n200l1.6s10$Dmax),
cvm = c(cvm.n200l1.6s10$bias.shape,cvm.n200l1.6s10$rmse.shape, cvm.n200l1.6s10$bias.scale, cvm.n200l1.6s10$rmse.scale, cvm.n200l1.6s10$Dobs, cvm.n200l1.6s10$Dmax),
ad = c(ad.n200l1.6s10$bias.shape,ad.n200l1.6s10$rmse.shape, ad.n200l1.6s10$bias.scale, ad.n200l1.6s10$rmse.scale, ad.n200l1.6s10$Dobs, ad.n200l1.6s10$Dmax),
rad = c(rtad.n200l1.6s10$bias.shape,rtad.n200l1.6s10$rmse.shape, rtad.n200l1.6s10$bias.scale, rtad.n200l1.6s10$rmse.scale, rtad.n200l1.6s10$Dobs, rtad.n200l1.6s10$Dmax)
)
                         
n200 <- n200l1.6s10
round(n200, 3)
t(apply(abs(n200), 1, rank))
# Sum of the ranks
apply(t(apply(abs(n200), 1, rank)), 2, sum)
rank(apply(t(apply(abs(n200), 1, rank)), 2, sum))

sink()
# End of simulation for param of real data

