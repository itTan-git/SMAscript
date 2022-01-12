# interclade comparison of the ratio isolated from clinical/hospital associated 

# script option ====
## library import ====
library(tidyverse)
library(data.table)
library(ggbeeswarm)
library(patchwork)

## graphic parameter ====
fontsize = 2
fontsize_theme = 2 * ggplot2::.pt
fontsize_legend = fontsize * 4
fig_width = 85
fig_height = 225

## function ====

# data import and convert ====
## sample data ====
df <- read_tsv("DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper",
         clade = "ANI clade",
         cluster = "SNP10 cluster",
         picked = "Representative strain in each SNP10 cluster",
         source = "clinical or hospital associated or nonclinical") %>% 
  mutate(source = ifelse(source == "hospital associated", "clinical", 
                         ifelse(source == "n.a.", "nonclinical", 
                                source))) %>% 
  arrange(clade)

### calculate clinical ratio in all strain ====
df_summary <- df %>% 
  group_by(clade) %>% 
  nest() %>% 
  mutate(clade = factor(clade),
         sample_num = map_dbl(data, nrow),
         clinical = map_dbl(data, function(x){nrow(filter(x, source == "clinical"))}),
         clinical_ratio = map_dbl(data, function(x){nrow(filter(x, source == "clinical")) / nrow(x) * 100}),
         nonclinical = map_dbl(data, function(x){nrow(filter(x, source == "nonclinical"))}),
         nonclinical_ratio = map_dbl(data, function(x){nrow(filter(x, source == "nonclinical")) / nrow(x) * 100})) %>% 
  select(-data)

### calculate clinical ratio in SNP10 cluster  ====
df_summary_SNP10 <- df %>% 
  filter(picked == "Rep") %>% 
  group_by(clade) %>% 
  nest() %>% 
  mutate(clade = factor(clade),
         sample_num = map_dbl(data, function(x){nrow(na.omit(x))}), 
         clinical = map_dbl(data, function(x){nrow(filter(x, source == "clinical"))}),
         clinical_ratio = map_dbl(data, function(x){nrow(filter(x, source == "clinical")) / nrow(x) * 100}),
         nonclinical = map_dbl(data, function(x){nrow(filter(x, source == "nonclinical"))}),
         nonclinical_ratio = map_dbl(data, function(x){nrow(filter(x, source == "nonclinical")) / nrow(x) * 100})) %>% 
  select(-data)

#### combine all data and SNP10 data ====
df_all <- df %>% 
  filter(picked == "Rep") %>% 
  summarise(clade = "all", 
            sample_num = nrow(.),
            clinical = nrow(filter(., source == "clinical")),
            clinical_ratio = nrow(filter(., source == "clinical")) / nrow(.) * 100, 
            nonclinical = nrow(filter(., source == "nonclinical")),
            nonclinical_ratio = nrow(filter(., source == "nonclinical")) / nrow(.) * 100) %>%  
  bind_rows(df_summary_SNP10, .) %>% 
  ungroup() %>% 
  mutate(clade = factor(x = clade, levels = c(1:14, "all")))

### statistical test to all clinical ratio ====
stat_data <- as.matrix(select(df_all, clinical, nonclinical))
fm_test <- RVAideMemoire::fisher.multcomp(stat_data, p.method = "BH")
fm_p_value_data <- as.data.frame(fm_test$p.value)

clinical_stat_result <- fm_p_value_data %>% 
  tail(1) %>% 
  data.table::transpose() %>% 
  mutate(clade = factor(x = 1:14)) %>% 
  select(clade, Pvalue = V1) %>% 
  mutate(show = case_when(Pvalue < 0.01 ~ "**",
                          Pvalue < 0.05 ~ "*")) %>% 
  mutate(y = df_all$clinical_ratio[1:14] + 2)

## AMR data ====
AMR <- read_tsv("DATA/final_dataset/AMR_info.tsv") %>% 
  select(-c(class, category, for_MDR_dicision)) %>%  
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  mutate_at(colnames(.)[-1], function(x){as.numeric(ifelse(x == "-", "0", "1"))}) %>% 
  left_join(x = df, y = .)

## calculate AMR gene number ====
AMR_incline <- AMR %>% 
  select(-c(sample_ID, clade, picked, source)) %>% 
  pivot_longer(cols = -cluster, names_to = "gene", values_to = "PA") %>% 
  group_by(cluster, gene) %>% 
  summarise(P_cluster_ratio = sum(PA)/n()) %>% 
  ungroup(gene) %>% 
  summarise(AMR_number = sum(P_cluster_ratio))

AMR_num <- df %>% 
  filter(picked == "Rep") %>% 
  left_join(AMR_incline, by = "cluster") %>% 
  mutate(clade = factor(x = clade, levels = 1:14))
 
## MDR data ====
MDR <- read_tsv("DATA/final_dataset/MDR_info.tsv") %>% 
  left_join(x = df, y = ., by = "sample_ID") %>% 
  select(sample_ID, clade, cluster, picked, resistance_num) %>% 
  mutate(resistance_num = replace_na(resistance_num, 0), 
         MDR = ifelse(resistance_num >= 3, "MDR", "notMDR")) 

## make data for x axis label ====
N_clade_all <- MDR %>% 
  group_by(clade, MDR) %>% 
  summarise(N = n()) %>% 
  pivot_wider(names_from = MDR, values_from = N) %>% 
  mutate(MDR = replace_na(MDR, 0), 
         sample_num = sum(MDR, notMDR),
         display_all = str_c('bgroup("(", frac(', MDR,', ', sample_num, '), ")")', sep = "")) 
  
N_clade_SNP <- MDR %>% 
  group_by(clade, cluster, MDR) %>% 
  summarise(N = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = MDR, values_from = N) %>% 
  mutate(MDRclade = ifelse(is.na(MDR), "notMDR", "MDR")) %>% 
  group_by(clade, MDRclade) %>% 
  summarise(N = n()) %>% 
  pivot_wider(names_from = MDRclade, values_from = N) %>% 
  mutate(MDR = replace_na(MDR, 0),
         sample_num = sum(MDR, notMDR),
         display_SNP = str_c("frac(", MDR, ", ", sample_num, ")", sep = "")) 

# graphical parameter ====
## ANI clade parameter ====
ANIbclade <- as.character(c(1:14))
ANIbcolor <- c("#dc143c", "#ff8c00", "#ff00ff", "#ffca80",	
               "#000000", "#ffb6c1", "#ffff00", "#005aff",	
               "#00bfff", "#00ffff", "#bfe4ff", "#03af7a", 
               "#00ff7f", "#a9a9a9") 

# draw graph ====
## clinical ratio bar plot ====
clinical_barplot <- df_all %>% 
  filter(clade != "all") %>% 
  ggplot(mapping = aes(x = clade, y = clinical_ratio, fill = clade)) +
  geom_hline(yintercept = tail(df_all$clinical_ratio, 1), linetype = "dashed") +
  geom_bar(stat = "identity")+ 
  geom_text(data = clinical_stat_result, mapping = aes(x = clade, y = y, label = show), 
            size = fontsize * 1.5, na.rm = TRUE) + 
  scale_fill_manual(values = ANIbcolor) + 
  scale_x_discrete(name = "", labels = NULL) + 
  scale_y_continuous(name = "clinical ratio (%)", limits = c(0,105)) + 
  guides(fill = "none") +  
  labs(tag = "A") +
  theme_bw() +
  theme(axis.title = element_text(size = fontsize_theme * 1.2), axis.text = element_text(size = fontsize_theme))

## AMR bee swarm plotの作成 ====
bp_AMR_num <- ggplot(data = AMR_num, mapping = aes(x = clade, y = AMR_number, fill = clade)) +
  geom_quasirandom(size = 1, shape = 21, colour = "#000000", stroke = .2) + 
  scale_fill_manual(values = ANIbcolor) + 
  scale_x_discrete(name = "", labels = NULL) + 
  scale_y_continuous(name = "AMR gene number", breaks = c(seq(0,max(AMR_num$AMR_number),4),26),
                     minor_breaks = NULL) + 
  guides(fill = "none") + 
  labs(tag = "B") +
  theme_bw() +
  theme(axis.title = element_text(size = fontsize_theme * 1.2), 
        axis.text = element_text(size = fontsize_theme))

## plot x axis label ====
X_label <- ggplot() +
  annotate(geom = "text", x = factor(1:14), y = -1, 
           label = 1:14, size = fontsize) + 
  annotate(geom = "text", x = factor(1:14), y = -3, 
           label = N_clade_SNP$display_SNP, size = fontsize, parse = TRUE) + 
  annotate(geom = "text", x = factor(1:14), y = -5, 
           label = N_clade_all$display_all, size = fontsize, parse = TRUE) +
  theme_void()

letter_plot <- clinical_barplot / bp_AMR_num / X_label
ggsave(filename = "RESULTS/clinical/F07ab_clade_clinical-AMR_plot.pdf", plot = letter_plot,
       device = "pdf", width = fig_width * 2, height = fig_height / 1.5, units = "mm", dpi = 1000)

# finalize with Illustrator
