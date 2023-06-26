# Plot single genotype
# Demonstration of logistic growth curve model

library(tidyverse)

# Raw data
season6 <- read.table(file = "data_clean/season6_combined.txt", 
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  filter(cultivar %in% c("PI655983", "PI213900"))

# Modeled parameters
s6 <- read_csv("data_clean/mac_growth_rate_modeled_season6.csv") %>%
  filter(genotype %in% c("PI655983", "PI213900")) %>%
  rename(cultivar = genotype) %>%
  mutate(label1 = paste0("Y[max]==", round(max_height_cm, 2)),
         label2 = paste0("R[half]==", round(max_growth_cm_gdd, 3)))

season6 %>%
  ggplot(aes(x = gdd, y = canopy_height)) +
  geom_point(size = 2,
             shape = 20) +
  geom_hline(data = s6, aes(yintercept = max_height_cm),
             lty = 2) +
  geom_abline(data = s6, aes(slope = max_growth_cm_gdd, 
                             intercept =  - (max_height_cm - min_height_cm)/2),
              lty = 2) +
  geom_text(data = s6, aes(label = label1,
                           x = 0, y = max_height_cm),
            parse = TRUE,
            hjust = 0,
            vjust = 1.5) +
  geom_text(data = s6, aes(label = label2,
                           x = 0, y = (max_height_cm - min_height_cm)/2),
            parse = TRUE,
            hjust = 0) +
  facet_wrap(~cultivar) +
  scale_x_continuous("Growing degree days",
                     breaks = c(0, 750, 1500)) +
  scale_y_continuous("Canopy height (cm)") +
  theme_bw(base_size = 12)
         

ggsave(filename = "data_figs/Fig2_growth.png",
       height = 3,
       width = 6,
       units = "in")
