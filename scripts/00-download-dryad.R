# Download original traits data from Dryad

# Download once
# download.file("https://datadryad.org/stash/downloads/file_stream/396628",
#               "data_raw/trait_data.zip",
#               mode = "wb")

# download.file("https://datadryad.org/stash/downloads/file_stream/396632",
#               "data_raw/sensor_data_catalogs.zip",
#               mode = "wb")

# List files
trait_files <- unzip("data_raw/trait_data.zip",
      list = TRUE)
trait_files$Name

sensor_files <- unzip("data_raw/sensor_data_catalogs.zip",
                     list = TRUE)
sensor_files$Name

# Unzip two files
s4 <- read.csv((unzip("data_raw/trait_data.zip",
    files = "traits/season_4_traits/season_4_canopy_height_sensor.csv",
    junkpaths = TRUE,
    exdir = "data_raw")))
                
s6 <- read.csv((unzip("data_raw/trait_data.zip",
                files = "traits/season_6_traits/season_6_canopy_height_sensor.csv",
                junkpaths = TRUE,
                exdir = "data_raw")))

# Check distribution of s4 genotypes among treatments
sm <- s4 %>%
  group_by(genotype, treatment) %>%
  summarize(n = n()) %>%
  tidyr::pivot_wider(names_from = treatment,
                     values_from = n) %>%
  mutate(trt = case_when(!is.na(`MAC Season 4: BAP water-deficit stress Aug 1-14`) & 
                           !is.na(`MAC Season 4: BAP water-deficit stress Aug 15-30`) ~ "both",
                         !is.na(`MAC Season 4: BAP water-deficit stress Aug 1-14`) ~ "early",
                         !is.na(`MAC Season 4: BAP water-deficit stress Aug 15-30`) ~ "late"))
sm %>%
  group_by(trt) %>%
  count()

