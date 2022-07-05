# Script version to test models
# Updated 7/2022 to account for
# a)  non-random block, was a stratified random that clustered short, medium and
# tall genotypes into distinct blocks
# b)  drought treatment in season 4 that should be a fixed factor

library(dplyr)
library(ggplot2)
library(rjags)
load.module('dic')
load.module('glm')
library(mcmcplots)
library(postjags)
library(cowplot)

#set seed for R's random sample
set.seed(5)
# read in season 6 wide format dataframe
season6 <- na.omit(read.table(file = "data_clean/season6_combined.txt", sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE))

season4 <- na.omit(read.table(file = "data_clean/season4_combined.txt", sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE))

#randomly sample 1 cultivar from the season for testing model
s6_cultivars <- sample(unique(season6$cultivar), size = 1)
s4_cultivars <- sample(unique(season4$cultivar), size = 5)

#subset  dataframe by 10 randomly selected cultivars
s6_subset <- season6 %>%  filter(cultivar %in% s6_cultivars) %>% 
  select(sitename, gdd, canopy_height, cultivar, date) %>% 
  arrange(date)

s4_subset <- season4 %>%  filter(cultivar %in% s4_cultivars) %>% 
  select(sitename, gdd, canopy_height, cultivar, date,
         first_water_deficit_treatment, second_water_deficit_treatment) %>% 
  arrange(date)

ggplot(data = s6_subset, 
       aes(x = gdd, y = canopy_height,
           color = cultivar,
           shape = sitename)) +
  geom_point() +
  geom_smooth() +
  guides(shape = "none")

ggplot(data = s4_subset, 
       aes(x = gdd, y = canopy_height,
           color = cultivar,
           group = sitename)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~second_water_deficit_treatment)

#nls model
c <- 90
a <- 0.1
y <- s6_subset$canopy_height[3]
g <- s6_subset$gdd[3]
b <- ((log((c/y) - 1)) - a)/g
model_single_cultivar <- nls(canopy_height ~ c / (1 + a* exp(b * gdd)), 
                             start = list(c = c, a = a, b = b),
                             data = s6_subset)
summary(model_single_cultivar)
coef(model_single_cultivar)

#nls model reparameterized version
c <- 300
a <- 1.5
b <- 0.5
model_reparam <- nls(canopy_height ~ c / (1 + exp((-4 * b)/c * gdd) * (c/a -1)), 
                             start = list(c = c, a = a, b = b),
                             data = s6_subset)
summary(model_reparam)
coef(model_reparam)

#data list
# s6_subset$block <- as.numeric(as.factor(s6_subset$sitename))
datlist <- list(height = s6_subset$canopy_height,
                gdd = s6_subset$gdd,
                # block = s6_subset$block,
                n = nrow(s6_subset)
                # stdc = 10, stda = 10, stdb = 10,
                # nblocks = length(unique(s6_subset$block))
                )

#initials list
# inits <- function(){list(mu.theta.c = rnorm(1, 0, 10), 
#                          mu.theta.a = rnorm(1, 0, 10),
#                          mu.theta.b = rnorm(1, 0, 10),
#                          tau.c.eps = runif(1, 0, 1),
#                          tau.a.eps = runif(1, 0, 1),
#                          tau.b.eps = runif(1, 0, 1),
#                          tau = runif(1, 0, 1))}

inits <- function(){list(theta.c = rnorm(1, 0, 10), 
                         theta.a = rnorm(1, 0, 10),
                         theta.b = rnorm(1, 0, 10),
                         tau = runif(1, 0, 1))}

initslist <- list(inits(), inits(), inits())

#initialize model
jm <- jags.model(file = "source/jags_simple.jags", 
                 data = datlist, 
                 inits = initslist,
                 n.chains = 3)

#set parameters to monitor
params <- c("deviance", "Dsum", 
            "theta.a", "theta.b", "theta.c",
            "tau", "sig",
            "Ymax", "Ymin", "Ghalf",
            "height.rep")

#update and monitor samples
update(jm, n.iter = 5000)

dic <- dic.samples(jm, n.iter = 5000)

jm_coda <- coda.samples(model = jm,
                        variable.names = params,
                        n.iter = 10000,
                        thin = 10)

#diagnostic plots via mcmcplots
mcmcplot(mcmcout = jm_coda,
         parms = c("deviance", "Dsum", "Ymax", "Ymin", "Ghalf",
                   "sig"))

#summarize posterior chains
post.sum <- coda.fast(jm_coda)

#match to data, plot model fit
pred <- data.frame(s6_subset, 
                   h.lower = post.sum[match("height.rep[1]", row.names(post.sum)):match("height.rep[73]", row.names(post.sum)),4],
                   h.upper = post.sum[match("height.rep[1]", row.names(post.sum)):match("height.rep[73]", row.names(post.sum)),5],
                   h.med = post.sum[match("height.rep[1]", row.names(post.sum)):match("height.rep[73]", row.names(post.sum)),2])

ggplot(pred, aes(x = canopy_height, y = h.med)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_pointrange(aes(ymin = h.lower, ymax = h.upper)) +
  scale_x_continuous("Observed height (cm)") +
  scale_y_continuous("Predicted height (cm)") +
  theme_cowplot() +
  coord_equal()

