# *Sorghum bicolor bicolor* parameters from Bayesian logistic regression

Preferred citation: Guo, Jessica S.; Ryan P. Bartelme; David S. LeBauer.(2021) "Sorghum bicolor bicolor parameters from Bayesian logistic regression".

Corresponding author: Jessica S. Guo, University of Arizona, [jessicaguo\@email.arizona.edu](mailto:jessicaguo@email.arizona.edu){.email}

License: CC0

[DOI](https://doi.org/10.5281/zenodo.7114676)

## Summary

The dense timeseries of canopy height from the gantry at Maricopa Agricultural Center capture the trajectory of plant growth, but result in too many variables to be used in other analyses. While simple summary statistics such as mean, max, and min can be easily calculated, they may not be the most biologically relevant and are subject to outliers. Therefore, we developed a cleaning algorithm to QA/QC the MAC canopy height data and a function to apply a Bayesian logistic regression to obtain height and growth parameters for each cultivar.

Season 4 and season 6 data were both available, but season 4 cultivars were subject to two late-season drought treatments. Due to lack of controls and unbalanced design, season 4 and season 6 parameters were not comparable.

## Model description

Automated measurements of plant height over the growing season were fit with a logistic growth curve to obtain the traits of maximum plant height and maximum growth rate relative to growing degree days (gdd, [McMaster & Wilhelm 1997](https://www.google.com/url?q=https://doi.org/10.1016/S0168-1923(97)00027-0&sa=D&source=docs&ust=1664228124277588&usg=AOvVaw2-5UHgL2m_Z-xosjdx37RY)). This approach can yield phenotypic traits that are less sensitive to outliers, more comparable across environmental conditions, and permit greater biological inference. First, we applied a cleaning algorithm which excluded cultivars with N \< 35 or density or had measurements on fewer than 40% of days which yielded a total of 326 cultivars.

Plant height and gdd for each cultivar were modeled in a Bayesian framework with the likelihood: 

$height_i \sim Normal(\mu_i, \sigma^2)$

where $\mu$ is the expected value of each height, $\sigma^2$ is the measurement error variance, and $i$ indexes each observation. The logistic model was defined as:

$\mu_i = \frac{c}{1 + a \cdot e^{b \cdot gdd_i}}$

With logistic parameters $a$, $b$, and $c$ and the covariate gdd associated with each observation. However, we reparameterized the model to obtain more biologically meaningful parameters by defining:

$Y_{max} = c$

$Y_{min} = \frac{c}{1+a}$

$R_{half} = - \frac{c \cdot b}{4}$

Logistic parameters were further transformed to the whole real line: 

$c = e^{\theta_c}$

$a = e^{\theta_a} - 1$

$b = - e^{\theta_b}$

All root nodes were given wide, relatively non-informative standard priors, including Normal(0, 1000) for all transformed logistic parameters ($\theta_a$, $\theta_b$, $\theta_c$) and Gamma(0.1, 0.1) for the measurement error precision ($1/\sigma^2$).

The above model was implemented in JAGS 4.3.0 ([Plummer, 2003](http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Drafts/Plummer.pdf) via R/rjags. Three parallel Markov chain Monte Carlo (MCMC) sequences were assigned dispersed starting values. Models were initially run until convergence was achieved at a Rubin and Gelman ([1992](https://www.jstor.org/stable/2246093)) diagnostic \< 1.2. An additional run of 10000 iterations per chain thinned by 10 yielded a total of 3000 relatively independent posterior samples for each parameter of interest. The posterior median of $Y_{max}$ (maximum height, cm) and $R_{half}$ (maximum growth rate, cm/gdd) of each cultivar are summarized [here](https://github.com/genophenoenvo/JAGS-logistic-growth/blob/main/data_clean/mac_growth_rate_modeled_season6.csv).

## Files and Folders

-   `logistic_growth_by_cultivar.Rmd` runs the fit_logistic_growth.R function for each cultivar in season 4 and season6. Cultivar SP1516 in season 6 was grown at 99 sites, including many with sparse measurements. Thus, it was modeled with the simple model. Model summaries are combined for both seasons into `data_clean/mac_growth_rate_modeled.csv`.

-   `logistic_growth_by_cultivar_season6.Rmd` runs the fit_logistic_growth.R function for each cultivar in season6 only. All cultivars were modeled with the simple model, as each both locations were treated as iid. Model summaries are available as `data_clean/mac_growth_rate_modeled_season6.csv`.

-   `benchmarking_logistic_growth.Rmd` tests the efficiency of the code for 10 cultivars.

#### `data_clean/`

-   `mac_growth_rate_modeled.csv` contains the posterior median of the three modeled parameters: minimum height (cm), maximum height (cm), and maximum growth rate (cm/gdd) for each season and cultivar. Model fit is described by the $R^2$ of the predicted versus observed canopy heights. Finally, the method and method type (RE or simple) are described. This provides the phenotype data about sorghum cultivars that can be added to the original dataset from [LeBauer et al 2020](https://datadryad.org/stash/dataset/doi:10.5061/dryad.4b8gtht99).

-   `mac_growth_rate_modeled_season.csv` contains the posterior median of the three modeled parameters: minimum height (cm), maximum height (cm), and maximum growth rate (cm/gdd) for each season and cultivar. Model fit is described by the $R^2$ of the predicted versus observed canopy heights. All cultivars were run with the simple model.

#### `scripts/`

-   \*`01-clean-data.R` takes MAC season 4 and MAC season 6 data hosted on CyVerse, extracts canopy_height, joins with weather data, and filters by cultivars with sample size of at least 35 or sampling density of at least 40% (density defined as sample size over sampling duration). This resulted in 336 cultivars from MAC season 4 and 326 cultivars from MAC season 6, located in `data_clean/`.

#### `source/`

-   `jags_canopy_height_model.R` is a test script used during model development.

-   `jags_simple.jags` is the JAGS model script for a simple logistic regression model. The data were modeled with a normal likelihood, with the mean described by a standard logistic model. Modeled variables were transformed to the full real line and given relatively uninformative normal priors, which improves convergence and mixing. Monitored variables were reparameterized to indicate minimum height, maximum height, and maximum growth rate. The global precision was given a diffuse gamma prior. Posterior predictive loss and replicated data (prediction interval) were calculated and monitored to assess model fit.

-   `jags_hierarchical.jags` is the model script for the hierarchical logistic regression model, which is similar to above but accounts for the random location of each cultivar. Here, the parameters for each location are drawn from a population-level distribution. The population-level means were given semi-informative normal priors, while the population-level precisions were drawn from a diffuse folded-t distribution. Only population-level mean parameters are monitored, to obtain similar output to `source/jags_simple.jags`.

-   `fit_logistic_growth.R` is a function that runs either the hierarchical or simple model for each cultivar, depending on whether site number was \> or = 1, respectively. Models were initialized with random starting values for 3 chains, 5000 samples discarded as burn-in, and 10000 samples retained. If the Gelman-Rubin (Rhat) diagnostic exceeded 1.2 for any parameter, the model was re-initialized with the ending values of the previous model run. If the Gelman-Rubin diagnostic still exceeded 1.2, the model was re-initialized with the ending values of the chain with the lowest posterior predictive loss (Dsum). Gelman-Rubin, posterior samples, traceplots, model fit, and model summary were saved out. Posterior samples were summarized as the median and central 95% credible interval for the final 10000 samples.

## Additional notes

The data and code are available on [GitHub](https://github.com/genophenoenvo/JAGS-logistic-growth).

**References**

LeBauer, David et al. (2020), Data From: TERRA-REF, An open reference data set from high resolution genomics, phenomics, and imaging sensors, Dryad, Dataset, <https://doi.org/10.5061/dryad.4b8gtht99>
