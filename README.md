# *Sorghum bicolor bicolor* parameters from Bayesian logistic regression

Preferred citation:  Guo, Jessica S.; Ryan P. Bartelme; David S. LeBauer.(2021) "Sorghum bicolor bicolor parameters from Bayesian logistic regression".  

Corresponding author: Jessica S. Guo, University of Arizona, jessicaguo@email.arizona.edu

License: CC0

DOI:

## Summary

The dense timeseries of canopy height from the gantry at Maricopa Agricultural Center capture the trajectory of plant growth, but result in too many variables to be used in other analyses. While simple summary statistics such as mean, max, and min can be easily calculated, they may not be the most biologically relevant and are subject to outliers. Therefore, we developed a cleaning algorithm to QA/QC the MAC canopy height data and a function to apply a Bayesian logistic regression to obtain height and growth parameters for each cultivar.  

## Files and Folders

 - `logistic_growth_by_cultivar.Rmd` runs the fit_logistic_growth.R function for each cultivar in season 4 and season6. Cultivar SP1516 in season 6 was grown at 99 sites, including many with sparse measurements. Thus, it was modeled with the simple model. Model summaries are combined for both seasons into `data_clean/mac_growth_rate_modeled.csv`.

 - `logistic_growth_by_cultivar_season6.Rmd` runs the fit_logistic_growth.R function for each cultivar in season6 only. All cultivars were modeled with the simple model, as each both locations were treated as iid. Model summaries are available as `data_clean/mac_growth_rate_modeled_season6.csv`.
 
 - `benchmarking_logistic_growth.Rmd` tests the efficiency of the code for 10 cultivars. 

#### `data_clean/`

 - `mac_growth_rate_modeled.csv` contains the posterior median of the three modeled parameters: minimum height (cm), maximum height (cm), and maximum growth rate (cm/gdd) for each season and cultivar. Model fit is described by the $R^2$ of the predicted versus observed canopy heights. Finally, the method and method type (RE or simple) are described. This provides the phenotype data about sorghum cultivars that can be added to the original dataset from [LeBauer et al 2020](https://datadryad.org/stash/dataset/doi:10.5061/dryad.4b8gtht99).
 
  - `mac_growth_rate_modeled_season.csv` contains the posterior median of the three modeled parameters: minimum height (cm), maximum height (cm), and maximum growth rate (cm/gdd) for each season and cultivar. Model fit is described by the $R^2$ of the predicted versus observed canopy heights. All cultivars were run with the simple model. 

#### `scripts/`

 - *`01-clean-data.R` takes MAC season 4 and MAC season 6 data hosted on CyVerse, extracts canopy_height, joins with weather data, and filters by cultivars with sample size of at least 35 or sampling density of at least 40% (density defined as sample size over sampling duration). This resulted in 336 cultivars from MAC season 4 and 326 cultivars from MAC season 6, located in `data_clean/`. 


#### `source/`
 - `jags_canopy_height_model.R` is a test script used during model development. 

 - `jags_simple.jags` is the JAGS model script for a simple logistic regression model. The data were modeled with a normal likelihood, with the mean described by a standard logistic model. Modeled variables were transformed to the full real line and given relatively uninformative normal priors, which improves convergence and mixing. Monitored variables were reparameterized to indicate minimum height, maximum height, and maximum growth rate. The global precision was given a diffuse gamma prior. Posterior predictive loss and replicated data (prediction interval) were calculated and monitored to assess model fit. 

 - `jags_hierarchical.jags` is the model script for the hierarchical logistic regression model, which is similar to above but accounts for the random location of each cultivar. Here, the parameters for each location are drawn from a population-level distribution. The population-level means were given semi-informative normal priors, while the population-level precisions were drawn from a diffuse folded-t distribution. Only population-level mean parameters are monitored, to obtain similar output to `source/jags_simple.jags`. 

 - `fit_logistic_growth.R` is a function that runs either the hierarchical or simple model for each cultivar, depending on whether site number was > or = 1, respectively. Models were initialized with random starting values for 3 chains, 5000 samples discarded as burn-in, and 10000 samples retained. If the Gelman-Rubin (Rhat) diagnostic exceeded 1.2 for any parameter, the model was re-initialized with the ending values of the previous model run. If the Gelman-Rubin diagnostic still exceeded 1.2, the model was re-initialized with the ending values of the chain with the lowest posterior predictive loss (Dsum). Gelman-Rubin, posterior samples, traceplots, model fit, and model summary were saved out. Posterior samples were summarized as the median and central 95% credible interval for the final 10000 samples. 

## Additional notes

The data and code are available on [GitHub](https://github.com/genophenoenvo/JAGS-logistic-growth). 

**References** 

LeBauer, David et al. (2020), Data From: TERRA-REF, An open reference data set from high resolution genomics, phenomics, and imaging sensors, Dryad, Dataset, https://doi.org/10.5061/dryad.4b8gtht99
