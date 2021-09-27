# Function to run the hierarchical RE or simple version of the JAGS logistic regression model
# Input dataframe should be for a single cultivar

fit_logistic_growth <- function(data, type = "simple", outdir) {
  
  require(rjags)
  load.module('dic')
  require(coda)
  require(postjags)
  require(mcmcplots)
  require(ggplot2)
  
  if (length(unique(data$cultivar)) != 1){
    stop("Only runs with a single cultivar.")
  }
  
  cultivar <- unique(data$cultivar)
  
  if (type == "RE" & length(unique(data$site)) == 1){
    type <- "simple"
    return("Number of sites < 1; running simple model without random effects.")
  }
  
  # Create outdir if does not already exist
  if(!dir.exists(paste0(outdir))){
    dir.create(paste0(outdir), recursive = T)
  }
  
  # Create cultivar folder if does not already exist
  if(!dir.exists(paste0(outdir, "/", cultivar))){
    dir.create(paste0(outdir, "/", cultivar), recursive = T)
  }
  
  # Create list of data for model
  if (type == "simple") {
    datlist <- list(height = data$canopy_height,
                    gdd = data$gdd,
                    n = nrow(data))
  } else {
    datlist <- list(height = data$canopy_height,
                    gdd = data$gdd,
                    block = data$site,
                    n = nrow(data),
                    stdc = 10, stda = 10, stdb = 10,
                    nblocks = length(unique(data$site)))
  }

  
  # Generate initials list
  if (type == "simple") {
    inits <- function(){list(theta.c = rnorm(1, 0, 10), 
                             theta.a = rnorm(1, 0, 10),
                             theta.b = rnorm(1, 0, 10),
                             tau = runif(1, 0, 1))}
  } else {
    inits <- function(){list(mu.theta.c = rnorm(1, 5, 1), 
                            mu.theta.a = rnorm(1, -5, 1),
                            mu.theta.b = rnorm(1, 5, 1),
                            tau.c.eps = runif(1, 0, 1),
                            tau.a.eps = runif(1, 0, 1),
                            tau.b.eps = runif(1, 0, 1),
                            tau = runif(1, 0, 1))}
  }
  
  initslist <- list(inits(), inits(), inits())
  
  # Initialize model
  if (type == "simple") {
    jm <- jags.model(file = "source/jags_simple.jags", 
                     data = datlist, 
                     inits = initslist,
                     n.chains = 3)
  } else {
    jm <- jags.model(file = "source/jags_hierarchical.jags", 
                     data = datlist, 
                     inits = initslist,
                     n.chains = 3)
    
  }
  
  # Print cultivar ID for initialized model
  print(paste0("Cultivar ", cultivar, "model initialized"))
  
  # Update and monitor samples for effective number of parameters
  update(jm, n.iter = 5000)
  # dic <- dic.samples(jm, n.iter = 5000)
  
  # Set parameters to monitor
  if (type == "simple") {
    params <- c("deviance", "Dsum", 
                "theta.a", "theta.b", "theta.c",
                "tau", "sig",
                "Ymax", "Ymin", "Ghalf",
                "height.rep")
  } else {
    params <- c("deviance", "Dsum", 
                "mu.theta.a", "mu.theta.b", "mu.theta.c", 
                "tau.a.eps", "tau.b.eps", "tau.c.eps",
                "tau", "sig", "sigs",
                "Ymax", "Ymin", "Ghalf",
                "ymax", "ymin", "ghalf",
                "height.rep")
  }
  # Parameters to check for convergence
  conv.params <- c("deviance", "Dsum", 
                   "tau", "sig",
                   "Ymax", "Ymin", "Ghalf")

  # Monitor coda samples
  jm_coda <- coda.samples(model = jm,
                          variable.names = params,
                          n.iter = 10000,
                          thin = 10)
  
  # First convergence check of key parameters: update and re-run if not converged
  gel <- gelman.diag(jm_coda, multivariate = F)
  gel2 <- gel$psrf[match(conv.params, row.names(gel$psrf)),1]
  n.update <- 0 
  if (max(gel2) > 1.3) {
    warning("model did not converge; restarting with saved state")
    n.update <- 1
    saved.state <- initfind(jm_coda)
    if (type == "simple") {
      new_inits <- removevars(saved.state, variables = c(1:6))[[2]]
      jm <- jags.model(file = "source/jags_simple.jags", 
                       data = datlist, 
                       inits = new_inits,
                       n.chains = 3)
    } else {
      new_inits <- removevars(saved.state, variables = c(1:6, 10:11, 16:17))[[2]]
      jm <- jags.model(file = "source/jags_hierarchical.jags", 
                       data = datlist, 
                       inits = new_inits,
                       n.chains = 3)
    }
    update(jm, n.iter = 5000)
    # dic <- dic.samples(jm, n.iter = 5000)
    jm_coda <- coda.samples(model = jm,
                            variable.names = params,
                            n.iter = 10000,
                            thin = 10)
    gel <- gelman.diag(jm_coda, multivariate = F)
    gel2 <- gel$psrf[match(conv.params, row.names(gel$psrf)),1]
  }
  
  # Second convergence check: update with lowest saved state and re-run if not converged
  if (max(gel2) > 1.3) {
    warning("model did not converge; restarting with saved state with lowest Dsum")
    n.update <- 2
    saved.state <- initfind(jm_coda)
    # Find chain with lowest Dsum
    i <- which(names(saved.state$initials[[1]]) == "Dsum")
    Dsums <- c(saved.state$initials[[1]][i],
               saved.state$initials[[2]][i],
               saved.state$initials[[3]][i])
    ind <- which.min(Dsums)
    if (type == "simple") {
      new_inits <- removevars(saved.state, variables = c(1:6))[[2]]
      new_inits2 <- list(new_inits[[ind]],
                         new_inits[[ind]],
                         new_inits[[ind]])
      jm <- jags.model(file = "source/jags_simple.jags", 
                       data = datlist, 
                       inits = new_inits2,
                       n.chains = 3)
    } else {
      new_inits <- removevars(saved.state, variables = c(1:6, 10:11, 16:17))[[2]]
      new_inits2 <- list(new_inits[[ind]],
                         new_inits[[ind]],
                         new_inits[[ind]])
      jm <- jags.model(file = "source/jags_hierarchical.jags", 
                       data = datlist, 
                       inits = new_inits2,
                       n.chains = 3)
    }
    update(jm, n.iter = 5000)
    # dic <- dic.samples(jm, n.iter = 5000)
    jm_coda <- coda.samples(model = jm,
                            variable.names = params,
                            n.iter = 10000,
                            thin = 10)
    gel <- gelman.diag(jm_coda, multivariate = F)
  }
  
  # Outputs and model diagnostics
  
  # Gelman diagnostic (Rhat)
  save(gel, file = file.path(outdir, cultivar,  paste0("gelman_", cultivar, ".Rdata")))
  
  # Coda and traceplots
  mcmc <- jm_coda
  save(mcmc, file = file.path(outdir, cultivar,  paste0("mcmc_", cultivar, ".Rdata")))
  
 
  jpeg(filename = file.path(outdir, cultivar,  paste0("trace_", cultivar, ".jpg")),
       width = 6, 
       height = 4,
       units = "in",
       res = 300)
  traplot(mcmc, parms = c("Ymax", "Ymin", "Ghalf", "sig"))
  dev.off()
  
  
  # Posterior summary
  mcmc_sum <- coda.fast(jm_coda)
  save(mcmc_sum, file = file.path(outdir, cultivar,  paste0("mcmc_sum_", cultivar, ".Rdata")))
  
  # Model fit
  inds <-grep("height.rep", rownames(mcmc_sum))
  pred <- data.frame(data, 
                    h.median = mcmc_sum[inds,2],
                    h.lower = mcmc_sum[inds,4],
                    h.upper = mcmc_sum[inds,5])
  
  fit <- ggplot(pred, aes(x = canopy_height, y = h.median)) +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    geom_pointrange(aes(ymin = h.lower, ymax = h.upper)) +
    scale_x_continuous("Observed height (cm)") +
    scale_y_continuous("Predicted height (cm)") +
    theme_bw(base_size = 14) +
    coord_equal()
  
  # Other fit parameters
  m <- summary(lm(h.median ~ canopy_height, data = pred))
  pred$cov <- ifelse(pred$canopy_height <= pred$h.upper & pred$canopy_height >= pred$h.lower,
                     1, 0)
  ggsave(filename = file.path(outdir, cultivar,  paste0("fit_", cultivar, ".jpg")),
         plot = fit,
         width = 4,
         height = 4, 
         units = "in")
  
  # Return a dataframe of cultivar, 
  # Posterior median and 95% CI of Ymax, Ymin, and Ghalf
  # Posterior mean of deviance and Dsum
  # pD (effective number of parameters) and DIC (deviance + pD) (Leave out for now)
  # Rhat (Gelman diagnostic),
  # R^2, bias, coverage
  out <- data.frame(cultivar = cultivar,
                    type = type, 
                    n.update = n.update, 
                    Ymax.median = mcmc_sum[grep("Ymax", rownames(mcmc_sum)),2],
                    Ymax.lower = mcmc_sum[grep("Ymax", rownames(mcmc_sum)),4],
                    Ymax.upper = mcmc_sum[grep("Ymax", rownames(mcmc_sum)),5],
                    Ymin.median = mcmc_sum[grep("Ymin", rownames(mcmc_sum)),2],
                    Ymin.lower = mcmc_sum[grep("Ymin", rownames(mcmc_sum)),4],
                    Ymin.upper = mcmc_sum[grep("Ymin", rownames(mcmc_sum)),5],
                    Ghalf.median = mcmc_sum[grep("Ghalf", rownames(mcmc_sum)),2],
                    Ghalf.lower = mcmc_sum[grep("Ghalf", rownames(mcmc_sum)),4],
                    Ghalf.upper = mcmc_sum[grep("Ghalf", rownames(mcmc_sum)),5],
                    deviance = mcmc_sum[grep("deviance", rownames(mcmc_sum)),1],
                    Dsum = mcmc_sum[grep("Dsum", rownames(mcmc_sum)),1],
                    # pD = sum(dic$penalty),
                    # DIC = sum(dic$deviance) + sum(dic$penalty),
                    Rhat = max(gel$psrf[,1]),
                    r2 = m$adj.r.squared,
                    bias = m$coefficients[2,1],
                    coverage = mean(pred$cov))
  save(out, file = file.path(outdir, cultivar, "out.Rdata"))
  
  return(out) 

}
