# Preparing height data for logistic growth model

# Load libraries
library(tidyverse)
library(lubridate)

# ================================================================
# 1) Trait data processing from MAC season 4 and season 6
# ================================================================

# First time only, use wget to obtain tall format trait data from Cyverse
#season four
# system('wget https://de.cyverse.org/dl/d/B3ADF887-BDE3-435B-9301-4C3FCB4F56F1/tall_season_four.csv -P data_raw/')
#season six (no change from raw data)
# system('wget https://de.cyverse.org/dl/d/FD84112F-FCEA-4089-8486-B1D19D71300B/mac_season_six_2020-04-22.csv -P data_raw/')

# Read in CSV filepaths as list
data_path <- list("data_raw/tall_season_four.csv", "data_raw/mac_season_six_2020-04-22.csv")

# Use lapply to read in csv as elements of list
  # note: function wrapper necessary to provide input argument to read.csv
raw_data <- lapply(data_path, FUN = function(i){
  read.csv(i, header=TRUE, stringsAsFactors = FALSE)})

# Name each dataframe in the list based on data_path list order
experiments <- c("mac_season_4", "mac_season_6")
names(raw_data) <- experiments

# Convert to tibble
raw_data <- map(.x = raw_data, .f = function(i){as_tibble(i)})

# wide format data function 
format_wide <- function(i){
  j <- i %>%
    mutate(row = row_number()) %>%
    pivot_wider(id_cols = c(row, lat, lon, date, sitename,
                            cultivar, treatment),
                names_from = trait, values_from = mean)   %>% 
    select(-row)
  return(j)
}

# Make a list of wide format tibbles
wide_trait_data <- vector(mode = "list", length = length(raw_data))
wide_trait_data <- lapply(raw_data, FUN = function(i){format_wide(i)})
names(wide_trait_data) <- experiments

# ================================================================
# 2) Select traits 
# ================================================================

# Create a vector of desired column names
data2use <- c("sitename", "date", "cultivar", "treatment", "canopy_height")

# Write function to select columns
select_data <- function(df){
  j <- as.data.frame(df[, (colnames(df) %in% data2use)])
  return(j)
}

# Eliminate extraneous data from datasets
filtered_trait_data <- vector(mode = "list", length = length(wide_trait_data))
filtered_trait_data <- map(.x = wide_trait_data, .f = function(df){select_data(df)})
names(filtered_trait_data) <- experiments

# ================================================================
# 3) filter by cultivars in all data sets (including genomic)
# ================================================================
# Read in cultivar lookup table
all_cult <- read.csv(file = "data_raw/cultivar_lookup_table.csv", header = TRUE,
                     stringsAsFactors = FALSE)

# First column is a character vector of all cultivars present across all seasons
# Make character vector of all cultivars in either mac season
all_cult$mac_total <- ifelse(all_cult$season_4 == 1 |all_cult$season_6 == 1, 1, 0)
cultivars4net <- as.vector(all_cult[all_cult$mac_total == 1, 1]) 

# 274 cultivars at 2 mac seasons with genomic data
# 353 cultivars from either mac season 

# Define filter cultivar function
filter_cultivar <- function(df){
  j <- as.data.frame(df[df$cultivar %in% cultivars4net, ])
  return(j)
}

# Filter by cultivars
trait_data <- vector(mode = "list", length = length(filtered_trait_data))
trait_data <- map(.x = filtered_trait_data, .f = function(df){filter_cultivar(df)})

# Remove all NA canopy heights
fix_canopy_height <- function(df){
  j <- as.data.frame(df[!is.na(df$canopy_height), ])
  return(j)
}
fixed_trait_data <- vector(mode = "list", length = length(trait_data))
fixed_trait_data <- map(.x = trait_data, .f = function(i){fix_canopy_height(i)})

# Convert data frames in list to tibbles
trait_tibbs <- vector(mode = "list", length = length(fixed_trait_data))
trait_tibbs <- map(.x = fixed_trait_data, .f = function(i){as_tibble(i)})

# Fix dates in datasets
for(i in 1:length(trait_tibbs)){
  trait_tibbs[[i]]$date <- as_date(trait_tibbs[[i]]$date)
}

# ================================================================
# 4) Join with weather data
# ================================================================

# First time only, use wget to obtain weather data from Cyverse
# system('wget https://de.cyverse.org/dl/d/6D959379-0442-41FE-8BEE-890866ACF037/mac_season_4_weather.csv -P data_raw/')
# system('wget https://de.cyverse.org/dl/d/C6219045-8114-4068-B924-8C2CD54AB9FD/mac_season_6_weather.csv -P data_raw/')

# List names of weather files
weather_raw <- list("data_raw/mac_season_4_weather.csv", "data_raw/mac_season_6_weather.csv")

# Use map from purrr to read in csv files
raw_weather_data <- map(.x = weather_raw, .f = function(i){read.csv(i, header = TRUE, stringsAsFactors = FALSE)})

#assign names to list of dataframes
names(raw_weather_data) <- c("mac_season_4_weather", "mac_season_6_weather")

# Check column names of weather data - why does it include deficit?
colnames(raw_weather_data$mac_season_4_weather)
colnames(raw_weather_data$mac_season_6_weather)

# Convert weather data list of df's to tibbles
weather_tibbs <- map(.x = raw_weather_data, .f = function(i){as_tibble(i)})

# Convert dates in list to date object for join
# Remove columns of water deficit 
for(i in 1:length(weather_tibbs)){
  weather_tibbs[[i]] <- weather_tibbs[[i]] %>%
    mutate(date = as_date(date)) %>%
    select(-first_water_deficit_treatment, -second_water_deficit_treatment)
}

# Join subset weather data and trait data, changed function to all weather data

combined_tibbs <- vector(mode = "list", length = length(trait_tibbs))

# Join traits and weather for both list objects
for(i in 1:length(trait_tibbs)){
  combined_tibbs[[i]] <- as.data.frame(left_join(trait_tibbs[[i]],
                            weather_tibbs[[i]], by = "date"), stringsasfactors = FALSE)
}
names(combined_tibbs) <- experiments

# Retain only unique sitename and date combinations
combined_tibbs_unq <- vector(mode = "list", length = length(combined_tibbs))

for(i in 1:length(trait_tibbs)){
  combined_tibbs_unq[[i]] <- combined_tibbs[[i]] %>% distinct(date, sitename, .keep_all = TRUE)
}
names(combined_tibbs_unq) <- experiments

# Check dimensions, appears that the size is cut in half for each dataframe
dim(combined_tibbs_unq[[2]])
dim(combined_tibbs[[2]])

# ================================================================
# 4) Filter cultivars by sample size
# ================================================================

# Summarize cultivar statistics by treatment
cult_sum <- vector(mode = "list", length = length(combined_tibbs_unq))
for(i in 1:2){
  cult_sum[[i]] <- combined_tibbs_unq[[i]] %>% 
    group_by(cultivar, treatment) %>%
    summarize(count = length(canopy_height),
              min_height = min(canopy_height),
              max_height = max(canopy_height),
              min_date = min(date),
              max_date = max(date),
              dur_days = as.numeric(difftime(max(date), min(date), "days")),
              min_gdd = min(gdd),
              max_gdd = max(gdd),
              dens = count / dur_days) %>%
    arrange(dens)
}

#mac season 4 range 2017-04-13 to 2017-09-21
range(combined_tibbs_unq[[1]]$date)

#mac season 6 range 2018-04-20 to 2018-08-02
range(combined_tibbs_unq[[2]]$date)

#combine into single dataframe and plot by cultivar
mac <- do.call(rbind.data.frame, combined_tibbs_unq) %>%
  mutate(season = case_when(date <= as.Date("2017-12-31") ~ "season 4",
                            date >= as.Date("2018-01-01") ~ "season 6"),
         cultivar = factor(cultivar, levels = unique(cult_sum[[1]]$cultivar)),
         trt = case_when(treatment == "MAC Season 4: BAP water-deficit stress Aug 1-14" ~ "early drought",
                         treatment == "MAC Season 4: BAP water-deficit stress Aug 15-30" ~ "late drought",
                         treatment == "MAC Season 6: Sorghum" ~ "no drought")) %>%
  arrange(desc(season))

# Split by cultivar, then plot
vars <- unique(cult_sum[[1]]$cultivar)
varslist <- split(vars, ceiling(seq_along(vars)/8))

for(i in 1:length(varslist)){
  sub <- subset(mac, cultivar %in% varslist[[i]])
  csum4 <- subset(cult_sum[[1]], cultivar %in% varslist[[i]]) %>%
    mutate(cultivar = factor(cultivar, levels = unique(cult_sum[[1]]$cultivar)),
           season = "season 4")
  csum6 <- subset(cult_sum[[2]], cultivar %in% varslist[[i]]) %>%
    mutate(cultivar = factor(cultivar, levels = unique(cult_sum[[1]]$cultivar)),
           season = "season 6")
  comb <- rbind(csum4, csum6) %>%
    mutate(trt = case_when(treatment == "MAC Season 4: BAP water-deficit stress Aug 1-14" ~ "early drought",
                           treatment == "MAC Season 4: BAP water-deficit stress Aug 15-30" ~ "late drought",
                           treatment == "MAC Season 6: Sorghum" ~ "no drought"),
           lat = case_when(trt %in% c("early drought", "no drought") ~ 0,
                           trt == "late drought" ~ 2000),
           lon = case_when(trt %in% c("early drought", "no drought") ~ 300,
                           trt == "late drought" ~ 50))
  fig <- ggplot() +
    geom_point(data = sub, aes(x = gdd, y = canopy_height,
                               color = trt), alpha = 0.5) +
    geom_text(data = comb, aes(x = lat, y = lon, label = round(dens,3),
                               color = trt),
              size = 3, hjust = 0, vjust = 0) +
    facet_grid(rows = vars(cultivar),
               cols = vars(season)) +
    theme_bw()
  fig
  
  jpeg(filename = paste0("data_figs/height_vs_gdd_", i, ".jpg"), 
       height = 5, width = 4, units = "in", res = 300)
  print(fig)
  dev.off()
}

# Use the filtering criteria of at least n = 35  OR dens > = 0.4
sel_cults <- vector(mode = "list", length = length(cult_sum))
for(i in 1:2){
  sel_cults[[i]] <- cult_sum[[i]] %>% 
    filter(count >= 35 | dens >= 0.4) %>%
    select(cultivar)
}
lapply(sel_cults, nrow)

# Results in 336 cultivars from MAC season 4 and 326 cultivars from MAC season 6

# Final selection of only cultivars that meet criteria
final_tibbs <- vector(mode = "list", length = length(combined_tibbs_unq))
for(i in 1:2){
  final_tibbs[[i]] <- combined_tibbs_unq[[i]] %>% 
    filter(cultivar %in% sel_cults[[i]]$cultivar)
}
names(final_tibbs) <- experiments
length(unique(final_tibbs[[1]]$cultivar))
length(unique(final_tibbs[[2]]$cultivar))
# Compare rows

lapply(final_tibbs, nrow)
lapply(combined_tibbs_unq, nrow)
# Same number of rows for season 6, fewer rows in season 4,

#write out season6 for growth curves
write.table(final_tibbs$mac_season_6, file = "data_clean/season6_combined.txt",
            quote = FALSE, sep = "\t")
write.table(final_tibbs$mac_season_4, file = "data_clean/season4_combined.txt",
            quote = FALSE, sep = "\t")
