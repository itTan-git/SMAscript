# summary of region and isolation source used in this study

# library ====
library(tidyverse)
library(ggnewscale)

# graphic parameter
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

# data import ====
sample_data <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(strainID = "Strain ID in this paper", 
         region = "Region", 
         clinical = "clinical or hospital associated or nonclinical",
         isolation = "isolation source") %>% 
  mutate(region = ifelse(is.na(region), "missing", region)) %>%
  mutate(region = ifelse(region == "Europe w/o UK", "Europe\nw/o UK", region))

# region info ====
## data manipulation ====
### all region data ====
region_summary <- sample_data %>% 
  select(region) %>% 
  group_by(region) %>%
  summarise(N = n()) %>%
  arrange(-N) %>%
  mutate(region = factor(region, levels = unique(region)),
         label = str_c(region, "\n(", N, ")", sep = ""),
         fraction = N/sum(N),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n = -1)),
         label_position = 1 - (ymax + ymin)/2)

### our data info ====
labo_strain_summary <- sample_data %>% 
  mutate(strain_from = case_when(str_starts(string = strainID, pattern = "DL") ~ "database",
                                 TRUE ~ "labo")) %>% 
  group_by(region, strain_from) %>% 
  summarise(N = n()) %>%
  left_join(select(region_summary, region, N), by = "region", suffix = c("_strain_from", "_region")) %>% 
  arrange(-N_region, strain_from) %>% 
  ungroup() %>% 
  mutate(region = factor(region, levels = unique(region)),
         region_from = str_c(region, strain_from, sep = "_"),
         region_from = factor(region_from,
                              levels = c("North America_database",
                                         "UK_database",
                                         "Japan_database",
                                         "Japan_labo",
                                         "Europe\nw/o UK_labo",
                                         "Europe\nw/o UK_database",
                                         "Asia/Oceania_database", 
                                         "South America_database",
                                         "Africa_database",
                                         "missing_database")),
         fraction = N_strain_from/sum(N_strain_from),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n = -1)))

## graphic parameter ==== 
### region
region <- levels(region_summary$region)
region_color <- c("#d8f255", "#ff4b00", "#005aff",
                  "#f6aa00", "#990099", "#03af7a",
                  "#fff100", "#83919e")

### our strain information 
strain_region <- levels(labo_strain_summary$region_from)
strain_region_color <- case_when(str_detect(string = strain_region, pattern = "database") ~ "#ffffff", TRUE ~ "#000000")
  

## graphic ====
c1_width = 2
c2_width = c1_width
c1_center = 1
space_c1c2 = .1
c2_center = c1_center + (c1_width / 2) + space_c1c2 + (c2_width / 2)
seg_position <- c(c1_center - (c1_width / 2),  c1_center + (c1_width / 2), 
                  c1_center + (c1_width / 2) + space_c1c2, 
                  c1_center + (c1_width / 2) + space_c1c2 + c2_width)

gg_region_info <- ggplot() +
  geom_bar(data = labo_strain_summary, 
           mapping = aes(x = c1_center, y = fraction, color = region, fill = region_from), 
           stat = "identity", width = c1_width, show.legend = FALSE) +
  scale_color_manual(values = rep(NA, length(region)), breaks = region, na.value = NA) +
  scale_fill_manual(values = strain_region_color, breaks = strain_region) +
  new_scale_fill() +
  new_scale_color() +
  geom_bar(data = region_summary,
           mapping = aes(x = c2_center, y = fraction, fill = region),
           stat = "identity", width = c2_width, color = "#ffffff", 
           size = .2, show.legend = FALSE) +
  scale_fill_manual(values = region_color, breaks = region) +
  geom_segment(mapping = aes(x = seg_position, y = 0, xend = seg_position, yend = 1), 
               color = "#000000", size = .1) +
  coord_polar(theta = "y", start = 0, direction = -1) +
  theme_void()

ggsave(filename = "F02a_region_summary.pdf", plot = gg_region_info, device = cairo_pdf, 
       path = "RESULTS/sample_summary", width = fig_width * 1.5, height = fig_height / 2.5, units = "mm", dpi = 1000)

# each region label are add with Illustrator

# clinical info ====
## data manipulation ====
### base data ====
ci_data <- sample_data %>% 
  mutate(clinical = ifelse(clinical == "n.a.", "nonclinical", clinical),
         isolation = ifelse(isolation == "clinical n.a." | isolation == "n.a.", "missing", isolation)) %>%
  group_by(clinical, isolation) %>% 
  mutate(isolation_display = ifelse(n() < 10 & isolation != "missing", "others", isolation)) %>% 
  ungroup()

### clinical data ====
source_type <- c("clinical", "hospital associated", "nonclinical")
clinical_summary <- ci_data %>% 
  group_by(clinical) %>% 
  summarise(N = n()) %>% 
  mutate(clinical = factor(clinical, levels = source_type),
         fraction = N/sum(N))

### isolation ====
isolation_summary <- ci_data %>% 
  group_by(clinical, isolation_display) %>% 
  summarise(N = n()) %>% 
  mutate(category = ifelse(isolation_display == "others", 2,
                           ifelse(isolation_display == "missing", 3, 1))) %>% 
  arrange(clinical, category, -N) %>% 
  ungroup() %>% 
  mutate(isolation_fake = str_c(clinical, isolation_display, sep = "_")) %>% 
  mutate(isolation_fake = factor(isolation_fake, levels = unique(isolation_fake)),
         fraction = N/sum(N)
  )

### labo strain ====
labo_clinical_summary <- ci_data %>% 
  mutate(strain_from = case_when(str_starts(string = strainID, pattern = "DL") ~ "database",
                                 TRUE ~ "labo")) %>% 
  group_by(clinical, strain_from) %>% 
  summarise(N = n()) %>% 
  ungroup() %>% 
  mutate(clinical_from = str_c(clinical, strain_from, sep = "_"),
         clinical_from = factor(clinical_from, levels = c("clinical_labo",
                                                          "clinical_database",
                                                          "hospital associated_database",
                                                          "nonclinical_database",
                                                          "nonclinical_labo")),
         fraction = N/sum(N))

## graphic parameter ====
source_color <- c("#ff4b00", "#ff8082", "#005aff")

isolation_type <- unique(isolation_summary$isolation_fake)
isolation_color <- c(rep("#ff4b00", 7), "#83919e", rep(c("#ff8082", "#005aff"), c(2,5)), "#83919e")

strain_clinical <- labo_clinical_summary$clinical_from
strain_clinical_color <- case_when(str_detect(string = strain_clinical, pattern = "database") ~ "#ffffff", TRUE ~ "#000000")

### pie chart of clinical info 
c1_width = 2
c2_width = c1_width / 4 * 3
c3_width = c1_width / 4
c1_center = 1
space_c1c2 = .1
c2_center = c1_center + (c1_width / 2) + space_c1c2 + (c2_width / 2)
c3_center = c1_center + (c1_width / 2) + space_c1c2 + (c2_width) + (c3_width / 2)
seg_position <- c(c1_center - (c1_width / 2),  c1_center + (c1_width / 2), 
                  c1_center + (c1_width / 2) + space_c1c2, 
                  c1_center + (c1_width / 2) + space_c1c2 + c2_width + c3_width)

gg_clinical_info <- ggplot() +
  geom_bar(data = labo_clinical_summary, mapping = aes(x = c1_center, y = fraction, color = clinical, fill = clinical_from),
           stat = "identity", width = c1_width, size = .2, 
           color = NA, position = "fill", show.legend = FALSE) + 
  scale_fill_manual(values = strain_clinical_color, breaks = strain_clinical) +
  new_scale_fill() +
  geom_bar(data = clinical_summary, mapping = aes(x = c2_center, y= fraction, fill = clinical),
           stat = "identity", width = c2_width, size = .2, 
           color = "#ffffff", position = "fill", show.legend = FALSE) +
  scale_fill_manual(values = source_color, breaks = source_type) +
  new_scale_fill() +
  geom_bar(data = isolation_summary, 
           mapping = aes(x = c3_center, y = fraction, fill = isolation_fake),
           stat = "identity", width = c3_width, size = .2, 
           color = "#ffffff", position = "fill", show.legend = FALSE) +
  scale_fill_manual(values = isolation_color, breaks = isolation_type) +
  geom_segment(mapping = aes(x = seg_position, y = 0, xend = seg_position, yend = 1), 
               color = "#000000", size = .1) +
  coord_polar(theta = "y", start = 0, direction = -1) +
  theme_void()

ggsave(filename = "F02b_isolation_summary.pdf", plot = gg_clinical_info, device = cairo_pdf, 
       path = "RESULTS/sample_summary", width = fig_width * 1.5, height = fig_height / 2.5, units = "mm", dpi = 1000)

# each region label are add with Illustrator