# logistic growth model
# random coefficients
# reparameterization  of model coefficients to be on the entire real line
# and reparameterized to Ymin, Ymax, and Ghalf
model{
  for(i in 1:n){
    # likelihood
    height[i] ~ dnorm(mu[i], tau)
    height.rep[i] ~ dnorm(mu[i], tau)
    
    # logistic expression
    mu[i] = c[block[i]] / (1 + a[block[i]] * exp(b[block[i]] * gdd[i]))
    
    # calculation for posterior predictive loss
    Sqdiff[i] <- pow(height.rep[i] - height[i], 2)
  }
  
  for(j in 1:nblocks){
    # reparameterization to real line
    c[j] <- exp(theta.c[j])
    a[j] <- exp(theta.a[j]) - 1
    b[j] <- -1 * exp(theta.b[j])
    
    # reparameterization to meaningful parameters
    ymax[j] <- c[j]
    ymin[j] <- c[j] / (1 + a[j])
    ghalf[j] <- -1 * (c[j] * b[j]) / 4
    
    theta.a[j] ~ dnorm(mu.theta.a, tau.a)
    theta.b[j] ~ dnorm(mu.theta.b, tau.b)
    theta.c[j] ~ dnorm(mu.theta.c, tau.c)
  }
  
  # root node priors - population means
  mu.theta.a ~ dnorm(5, 0.1)
  mu.theta.b ~ dnorm(-5, 0.1)
  mu.theta.c ~ dnorm(5, 0.1)
  
  # compute population-level values for phi's
  phi.c <- exp(mu.theta.c)
  phi.a <- exp(mu.theta.a) - 1
  phi.b <- -1 * exp(mu.theta.b)
  
  # compute population-level values for reparameterized Ymax, Ymin, and Ghalf
  Ymax <- phi.c
  Ymin <- phi.c / (1 + phi.a)
  Ghalf <-  -1 * (phi.c * phi.b) / 4
  
  # root node priors - population precisions
  # folded t for small group size
  tau.a.eps ~ dt(0, Pa, 2)
  sig.a.eps <- abs(tau.a.eps)
  tau.a <- pow(sig.a.eps, -2)

  tau.b.eps ~ dt(0, Pb, 2)
  sig.b.eps <- abs(tau.b.eps)
  tau.b <- pow(sig.b.eps, -2)

  tau.c.eps ~ dt(0, Pc, 2)
  sig.c.eps <- abs(tau.c.eps)
  tau.c <- pow(sig.c.eps, -2)

  # set as data
  Pa <- 1/stda*stda
  Pb <- 1/stdb*stdb
  Pc <- 1/stdc*stdc
  
  # root node prior - global precision
  tau ~ dgamma(0.1, 0.1)
  sig <- pow(tau, -0.5)
  
  # sigs to monitor
  sigs[1] <- sig
  sigs[2] <- sig.a.eps
  sigs[3] <- sig.b.eps
  sigs[4] <- sig.c.eps

  # posterior predictive loss is the posterior mean of Dsum
  Dsum <- sum(Sqdiff[])
}