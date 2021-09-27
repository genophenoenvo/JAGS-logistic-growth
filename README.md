# *Sorghum bicolor bicolor* parameters from Bayesian logistic regression

The dense timeseries of canopy height from the gantry at Maricopa Agricultural Center capture the trajectory of plant growth, but result in too many variables to be used in other analyses. While simple summary statistics such as mean, max, and min can be easily calculated, they may not be the most biologically relevant and are subject to outliers. Therefore, we developed a cleaning algorithm to QA/QC the MAC canopy height data in `phenophasebbn/bnprocess_mac.R` and a set of scripts and functions to apply a Bayesian logistic regression to obtain height and growth parameters for each cultivar `phenophasebbn/jags/`. 

## Data prep and QA/QC
`phenophasebbn/bnprocess_mac.R` takes MAC season 4 and MAC season 6 data processed by Emily's Jupyter notebooks, extracts canopy_height, and filters by cultivars present in both seasons that have corresponding genomic data (n = 274 cultivars). Additional filtering steps include removing NA values and joining weather data by date. 

`phenophasebbn/data_figs/`contains MAC season 4 (2017-04-13 to 2017-09-21) and season 6 (2018-04-20 to 2018-08-02) timeseries of height vs gdd by cultivar. 

Criteria for inclusion in the Bayesian logistic regression were if n observations >= 35 OR if density (heights per day) >= 0.4. This resulted in 262 cultivars from MAC season 4 and 274 cultivars from MAC season 6, saved separately as `phenophasebbn/season4_combined.txt` and `phenophasebbn/season6_combined.txt`; `phenophasebbn/bn_input.txt` was the combined filtered data for both MAC seasons. 

## Bayesian logistic regression
`jags/jags_canopy_height_model.R` is a test script used during model development. 

`jags/jags_simple.R` is the model script for a simple logistic regression model. The data were modeled with a normal likelihood, with the mean described by a standard logistic model. Modeled variables were reparameterized on the full real line and given relatively uninformative normal priors, which improves convergence and mixing. Monitored variables were reparameterized to indicate minimum height, maximum height, and maximum growth rate. The global precision was given a diffuse gamma prior. Posterior predictive loss and replicated data (prediction interval) were calculated and monitored to assess model fit. 

`jags/jags_hierarchical.R` is the model script for the hierarchical logistic regression model, which is similar to above but accounts for the random block design. Here, the parameters for each block are drawn from a population-level distribution. The population-level means were given semi-informative normal priors, while the population-level precisions were drawn from a diffuse folded-t distribution. Only population-level mean parameters are monitored, to obtain similar output to `jags/jags_simple.R`. 

`jags/fit_logistic_growth.R` is a function that runs either the hierarchical or simple model for each cultivar, depending on whether site number was > or = 1, respectively. Models were initialized with random starting values for 3 chains, 5000 samples discarded as burn-in, and 10000 samples retained. If the Gelman-Rubin (Rhat) diagnostic exceeded 1.3 for any parameter, the model was re-initialized with the ending values of the previous model run. If the Gelman-Rubin diagnostic still exceeded 1.3, the model was re-initialized with the ending values of the chain with the lowest posterior predictive loss (Dsum). Gelman-Rubin, posterior samples, traceplots, model fit, and model summary were saved out. Posterior samples were summarized as the median and central 95% credible interval for the final 10000 samples. 

`jags/benchmarking_logistic_growth.Rmd` tests the efficiency of the code for 10 cultivars. 

`jags/logistic_growth_by_cultivar.Rmd` runs the fit_logistic_growth.R function for each cultivar in season 4 and season6; model sumamries are combined for both seasons into `jags/mac_growth_rate_modeled.csv`. 

## Output
`jags/mac_growth_rate_modeled.csv` contains the posterior median of the three modeled parameters: minimum height (cm), maximum height (cm), and maximum growth rate (cm/gdd) for each season and cultivar. Model fit is described by the $R^2$ of the predicted versus observed canopy heights. Finally, the method and method type (RE or simple) are described. This provides the phenotype data about sorghum cultivars needed for the PEG Bayesian Belief Network. 

