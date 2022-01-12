# library ====
library(tidyverse)
library(data.table)

# data import ====
AMR <- read_tsv(file = "DATA/final_dataset/AMR_info.tsv")
sample_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper", clade = "ANI clade", cluster = "SNP10 cluster")

# data convert ====
## single gene distribution ====
### ESBL ====
ESBL_single <- AMR %>% 
  filter(category == "ESBL") %>% 
  select(-c(class, category, for_MDR_dicision)) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  pivot_longer(cols = -sample_ID, names_to = "gene", values_to = "variation") %>% 
  mutate(variation = ifelse(variation == "+", gene, variation)) %>% 
  pivot_wider(names_from = "gene", values_from = "variation")

### Carbapenemase ====
Car_single <- AMR %>% 
  filter(category == "Car") %>% 
  select(-c(class, category, for_MDR_dicision)) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  pivot_longer(cols = -sample_ID, names_to = "gene", values_to = "variation") %>% 
  mutate(variation = ifelse(variation == "+", gene, variation)) %>% 
  pivot_wider(names_from = "gene", values_from = "variation")

## mutiple gene distribution ====
### ESBL ====
ESBL_multi <- ESBL_single %>% 
  pivot_longer(cols = -sample_ID, names_to = "genes", values_to = "variation") %>% 
  mutate(genes = factor(x = genes, levels = c("TEM-1D", "SHV-OKP-LEN", "CTX-M-1", "CTX-M-2", "CTX-M-9"))) %>% 
  arrange(genes) %>% 
  filter(variation != "-") %>% 
  group_by(sample_ID) %>% 
  summarise(ESBL_pattern = str_c(genes, collapse = " + ")) %>% 
  mutate(multi_ESBL_pattern = ifelse(str_detect(ESBL_pattern, " "), ESBL_pattern, "-"),
         multi_ESBL = ifelse(multi_ESBL_pattern == "-", "-", "multi_ESBL"))

### Carbapenemase ====
Car_multi <- Car_single %>% 
  pivot_longer(cols = -sample_ID, names_to = "genes", values_to = "variation") %>% 
  mutate(genes = factor(x = genes, levels = c("KPC-1", "NDM-1", "SME-1", "GIM-1", "OXA-48", "IMP-1", "VIM-1"))) %>% 
  arrange(genes) %>% 
  filter(variation != "-") %>% 
  group_by(sample_ID) %>% 
  summarise(Car_pattern = str_c(genes, collapse = " + ")) %>% 
  mutate(multi_Car_pattern = ifelse(str_detect(Car_pattern, " "), Car_pattern, "-"),
         multi_Car = ifelse(multi_Car_pattern == "-", "-", "multi_Car"))

## combine single and multi data ====
ESBL_carba <- ESBL_single %>% 
  left_join(ESBL_multi) %>% 
  left_join(Car_single) %>% 
  left_join(Car_multi) %>% 
  mutate_all(~replace_na(., "-"))
  
## ESBL and Carbapenemase combination pattern ====
ESBL_carba_combination <- ESBL_carba %>% 
  select(sample_ID, ESBL_pattern, Car_pattern) %>% 
  mutate(ESBL_pattern = factor(ESBL_pattern, 
                               levels = c("TEM-1D", "CTX-M-1", "CTX-M-9", 
                                          "SHV-OKP-LEN", "TEM-1D + CTX-M-1", 
                                          "TEM-1D + CTX-M-2", "TEM-1D + SHV-OKP-LEN", 
                                          "TEM-1D + SHV-OKP-LEN + CTX-M-1", "-")),
         Car_pattern = factor(Car_pattern,
                              levels = c("KPC-1", "NDM-1", "SME-1", "GIM-1", 
                                         "OXA-48", "IMP-1", "VIM-1", 
                                         "NDM-1 + OXA-48", "-"))) %>% 
  arrange(ESBL_pattern, Car_pattern) %>% 
  pivot_longer(cols = -sample_ID, names_to = "pattern", values_to = "variation") %>% 
  group_by(sample_ID) %>% 
  summarise(combination_pattern = str_c(str_subset(string = variation, pattern = "^-$", negate = TRUE), collapse = " & ")) %>% 
  mutate(combination_pattern = ifelse(str_detect(string = combination_pattern, pattern = "&"), combination_pattern, "-"),
         combination = ifelse(str_detect(string = combination_pattern, pattern = "&"), "combination", "-"))
  
## data for analysis ====
ESBL_carba_count <- ESBL_carba %>% 
  left_join(ESBL_carba_combination) %>%
  mutate(ESBL_pattern = ifelse(ESBL_pattern == "-", "-", "ESBL"),
         Car_pattern = ifelse(Car_pattern == "-", "-", "Carbapenemase")) %>% 
  left_join(x = sample_info, y = .)

summary_category <- colnames(ESBL_carba_count)[-c(1:3)]

### strain ====
#### function ====
summarise_all <- function(x){
  result <- ESBL_carba_count %>% 
    mutate(all = "all") %>% 
    select(all, 
           category = all_of(x)) %>% 
    filter(category != "-") %>% 
    group_by(all, category) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = category, values_from = N)
}
summarise_clade <- function(x){
  result <- ESBL_carba_count %>% 
    select(clade,
           category = all_of(x)) %>% 
    filter(category != "-") %>% 
    group_by(clade, category) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = category, values_from = N)
}

#### count sample number ====
beta_all <- ESBL_carba_count %>% 
  mutate(all = "all") %>% 
  group_by(all) %>% 
  summarise(N = n()) %>% 
  left_join(summarise_all(x = summary_category[1])) %>% 
  left_join(summarise_all(x = summary_category[2])) %>% 
  left_join(summarise_all(x = summary_category[3])) %>% 
  left_join(summarise_all(x = summary_category[4])) %>% 
  left_join(summarise_all(x = summary_category[5])) %>% 
  left_join(summarise_all(x = summary_category[6])) %>% 
  left_join(summarise_all(x = summary_category[7])) %>% 
  left_join(summarise_all(x = summary_category[8])) %>% 
  left_join(summarise_all(x = summary_category[9])) %>% 
  left_join(summarise_all(x = summary_category[10])) %>% 
  left_join(summarise_all(x = summary_category[11])) %>% 
  left_join(summarise_all(x = summary_category[12])) %>% 
  left_join(summarise_all(x = summary_category[13])) %>% 
  left_join(summarise_all(x = summary_category[14])) %>% 
  left_join(summarise_all(x = summary_category[15])) %>% 
  left_join(summarise_all(x = summary_category[16])) %>% 
  left_join(summarise_all(x = summary_category[17])) %>% 
  left_join(summarise_all(x = summary_category[18])) %>% 
  left_join(summarise_all(x = summary_category[19])) %>% 
  left_join(summarise_all(x = summary_category[20])) %>% 
  transpose(keep.names = "mutation", make.names = "all")

beta_clade <- sample_info %>% 
  group_by(clade) %>% 
  summarise(N = n()) %>% 
  left_join(summarise_clade(x = summary_category[1])) %>% 
  left_join(summarise_clade(x = summary_category[2])) %>% 
  left_join(summarise_clade(x = summary_category[3])) %>% 
  left_join(summarise_clade(x = summary_category[4])) %>% 
  left_join(summarise_clade(x = summary_category[5])) %>% 
  left_join(summarise_clade(x = summary_category[6])) %>% 
  left_join(summarise_clade(x = summary_category[7])) %>% 
  left_join(summarise_clade(x = summary_category[8])) %>% 
  left_join(summarise_clade(x = summary_category[9])) %>% 
  left_join(summarise_clade(x = summary_category[10])) %>% 
  left_join(summarise_clade(x = summary_category[11])) %>% 
  left_join(summarise_clade(x = summary_category[12])) %>% 
  left_join(summarise_clade(x = summary_category[13])) %>% 
  left_join(summarise_clade(x = summary_category[14])) %>% 
  left_join(summarise_clade(x = summary_category[15])) %>% 
  left_join(summarise_clade(x = summary_category[16])) %>% 
  left_join(summarise_clade(x = summary_category[17])) %>% 
  left_join(summarise_clade(x = summary_category[18])) %>% 
  left_join(summarise_clade(x = summary_category[19])) %>% 
  left_join(summarise_clade(x = summary_category[20])) %>% 
  transpose(keep.names = "mutation", make.names = "clade")

#### output data ==== 
beta_summary <- left_join(beta_all, beta_clade) %>% 
  mutate(category = factor(c("all", 
                             rep("ESBL",6), 
                             rep("ESBL_multi", 5), 
                             rep("Car", 8), 
                             rep("Car_multi", 2), 
                             rep("combination", 11)),
                           levels = c("all", "ESBL", "ESBL_multi", "Car", "Car_multi", "combination")),
         order = c(1,
                   4,5,6,3,2,1,
                   3,4,2,5,1,
                   4,5,7,2,3,8,6,1,
                   2,1,
                   3,4,5,7,6,2,10,8,9,11,1)) %>% 
  arrange(category, order) %>% 
  select(category, mutation, all, as.character(1:14)) %>% 
  ungroup()

### SNP10 cluster ====
#### function ====
SNP_summarise_all <- function(x){
  result <- ESBL_carba_count %>% 
    mutate(all = "all") %>% 
    select(all, cluster,
           category = all_of(x)) %>% 
    filter(category != "-") %>% 
    group_by(all, cluster, category) %>% 
    summarise(N = n()) %>% 
    group_by(all, category) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = category, values_from = N)
}
SNP_summarise_clade <- function(x){
  result <- ESBL_carba_count %>% 
    select(clade, cluster,
           category = all_of(x)) %>% 
    filter(category != "-") %>% 
    group_by(clade, cluster, category) %>% 
    summarise(N = n()) %>% 
    group_by(clade, category) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = category, values_from = N)
}

#### count sample number ====
beta_all_SNP <- ESBL_carba_count %>% 
  mutate(all = "all") %>% 
  group_by(all, cluster) %>% 
  summarise(N = n()) %>% 
  group_by(all) %>% 
  summarise(N = n()) %>% 
  left_join(SNP_summarise_all(x = summary_category[1])) %>% 
  left_join(SNP_summarise_all(x = summary_category[2])) %>% 
  left_join(SNP_summarise_all(x = summary_category[3])) %>% 
  left_join(SNP_summarise_all(x = summary_category[4])) %>% 
  left_join(SNP_summarise_all(x = summary_category[5])) %>% 
  left_join(SNP_summarise_all(x = summary_category[6])) %>% 
  left_join(SNP_summarise_all(x = summary_category[7])) %>% 
  left_join(SNP_summarise_all(x = summary_category[8])) %>% 
  left_join(SNP_summarise_all(x = summary_category[9])) %>% 
  left_join(SNP_summarise_all(x = summary_category[10])) %>% 
  left_join(SNP_summarise_all(x = summary_category[11])) %>% 
  left_join(SNP_summarise_all(x = summary_category[12])) %>% 
  left_join(SNP_summarise_all(x = summary_category[13])) %>% 
  left_join(SNP_summarise_all(x = summary_category[14])) %>% 
  left_join(SNP_summarise_all(x = summary_category[15])) %>% 
  left_join(SNP_summarise_all(x = summary_category[16])) %>% 
  left_join(SNP_summarise_all(x = summary_category[17])) %>% 
  left_join(SNP_summarise_all(x = summary_category[18])) %>% 
  left_join(SNP_summarise_all(x = summary_category[19])) %>% 
  left_join(SNP_summarise_all(x = summary_category[20])) %>% 
  transpose(keep.names = "mutation", make.names = "all")

beta_clade_SNP <- sample_info %>% 
  group_by(clade, cluster) %>% 
  summarise(N = n()) %>% 
  group_by(clade) %>% 
  summarise(N = n()) %>% 
  left_join(SNP_summarise_clade(x = summary_category[1])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[2])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[3])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[4])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[5])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[6])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[7])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[8])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[9])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[10])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[11])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[12])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[13])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[14])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[15])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[16])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[17])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[18])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[19])) %>% 
  left_join(SNP_summarise_clade(x = summary_category[20])) %>% 
  transpose(keep.names = "mutation", make.names = "clade")

#### output data ==== 
beta_summary_SNP <- left_join(beta_all_SNP, beta_clade_SNP) %>% 
  mutate(category = factor(c("all", 
                             rep("ESBL",6), 
                             rep("ESBL_multi", 5), 
                             rep("Car", 8), 
                             rep("Car_multi", 2), 
                             rep("combination", 11)),
                           levels = c("all", "ESBL", "ESBL_multi", "Car", "Car_multi", "combination")),
         order = c(1,
                   4,5,6,3,2,1,
                   3,4,2,5,1,
                   4,5,7,2,3,8,6,1,
                   2,1,
                   3,4,5,7,6,2,10,8,9,11,1)) %>% 
  arrange(category, order) %>% 
  select(mutation, all, as.character(1:14)) %>% 
  ungroup()

### combine strain result and SNP10 cluster result ====
beta_summary_output <- beta_summary %>% 
  right_join(x = beta_summary_SNP, y = ., by = "mutation") %>% 
  mutate(category = as.character(category)) %>% 
  mutate_all(~replace_na(., 0)) %>% 
  mutate(
    all = str_c(all.x, " (", all.y, ")", sep = ""),
    `1` = str_c(`1.x`, " (", `1.y`, ")", sep = ""),
    `2` = str_c(`2.x`, " (", `2.y`, ")", sep = ""),
    `3` = str_c(`3.x`, " (", `3.y`, ")", sep = ""),
    `4` = str_c(`4.x`, " (", `4.y`, ")", sep = ""),
    `5` = str_c(`5.x`, " (", `5.y`, ")", sep = ""),
    `6` = str_c(`6.x`, " (", `6.y`, ")", sep = ""),
    `7` = str_c(`7.x`, " (", `7.y`, ")", sep = ""),
    `8` = str_c(`8.x`, " (", `8.y`, ")", sep = ""),
    `9` = str_c(`9.x`, " (", `9.y`, ")", sep = ""),
    `10` = str_c(`10.x`, " (", `10.y`, ")", sep = ""),
    `11` = str_c(`11.x`, " (", `11.y`, ")", sep = ""),
    `12` = str_c(`12.x`, " (", `12.y`, ")", sep = ""),
    `13` = str_c(`13.x`, " (", `13.y`, ")", sep = ""),
    `14` = str_c(`14.x`, " (", `14.y`, ")", sep = "")) %>% 
  select(category, mutation, all, as.character(1:14)) %>% 
  mutate_all(~ifelse(. == "0 (0)", "", .)) %>% 
  write_tsv(file = "RESULTS/AMR/AMR_summary_table_beta.tsv")
