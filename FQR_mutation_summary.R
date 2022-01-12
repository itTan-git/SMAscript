# library ====
library(tidyverse)
library(data.table)

# data import ====
AMR <- read_tsv(file = "DATA/final_dataset/AMR_info.tsv")
sample_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper", clade = "ANI clade", cluster = "SNP10 cluster")

# FQR summary table ====
## point mutation  ====
### data import ====
#### point mutation ====
FQR_PM <- AMR %>% 
  filter(class == "Gyrase point mutation") %>% 
  select(-c(class, category, for_MDR_dicision)) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  left_join(x = sample_info, y = .) %>% 
  mutate(gene_gyrA = case_when(gyrA_81 == "-" & gyrA_83 == "-" & gyrA_87 == "-" ~ "-",
                          TRUE ~ "gyrA"),
         gene_parC = ifelse(parC == "-", "-", "parC"),
         gene_parE = ifelse(parE == "-", "-", "parE"))

##### multi-mutation in gyrA ====
FQR_PM_gyrA <- FQR_PM %>% 
  select(sample_ID, starts_with("gyrA_")) %>% 
  pivot_longer(cols = -sample_ID, names_to = "position", values_to = "variation") %>% 
  group_by(sample_ID) %>% 
  summarise(multi_gyrA_pattern = str_c(str_subset(string = variation, pattern = "-", negate = TRUE), collapse = " + ")) %>% 
  mutate(multi_gyrA_pattern = ifelse(str_detect(multi_gyrA_pattern, " "), multi_gyrA_pattern, "-"),
         multi_gyrA = ifelse(multi_gyrA_pattern == "-", "-", "multi_gyrA"))

##### multi-mutation all gene ====
FQR_PM_multi_gene <- FQR_PM %>% 
  select(sample_ID, starts_with("gene_")) %>% 
  pivot_longer(cols = -sample_ID, names_to = "gene", values_to = "variation") %>% 
  group_by(sample_ID) %>% 
  summarise(multi_gene = str_c(str_subset(string = variation, pattern = "-", negate = TRUE), collapse = " + ")) %>% 
  mutate(multi_gene = ifelse(str_detect(multi_gene, " "), multi_gene, "-"))

#### quinolone resistance gene ====
quinolone <- AMR %>% 
  filter(class == "Quinolone") %>% 
  select(-c(class, category, for_MDR_dicision)) %>% 
  transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
  pivot_longer(cols = -sample_ID, names_to = "gene", values_to = "variation") %>% 
  mutate(variation = ifelse(variation == "-", 0, 1)) %>% 
  group_by(sample_ID) %>% 
  summarise(quinolone = sum(variation)) %>% 
  mutate(quinolone = ifelse(quinolone == 0, "-", "quinolone"))

##### combination with point mutation ====  
FQR_combination <- FQR_PM %>% 
  select(sample_ID, starts_with("gene_")) %>% 
  pivot_longer(cols = -sample_ID, names_to = "gene", values_to = "variation") %>% 
  group_by(sample_ID) %>% 
  summarise(multi_gene = str_c(str_subset(string = variation, pattern = "-", negate = TRUE), collapse = " + ")) %>% 
  left_join(quinolone) %>% 
  mutate(noPM = ifelse(multi_gene == "" & quinolone == "quinolone", "w/o FQR mutation", "-")) %>% 
  select(sample_ID, noPM)

#### data for analysis ====
FQR_PM_data <- FQR_PM %>% 
  left_join(FQR_PM_gyrA) %>% 
  left_join(FQR_PM_multi_gene) %>% 
  left_join(quinolone) %>% 
  left_join(FQR_combination)

summary_category <- colnames(FQR_PM_data)[-c(1:3)]

### strain ====
#### function ====
PM_summarise_all <- function(x){
  result <- FQR_PM_data %>% 
    mutate(all = "all") %>% 
    select(all, 
           gene = all_of(x)) %>% 
    filter(gene != "-") %>% 
    group_by(all, gene) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = gene, values_from = N)
}
PM_summarise_clade <- function(x){
  result <- FQR_PM_data %>% 
    select(clade,
           gene = all_of(x)) %>% 
    filter(gene != "-") %>% 
    group_by(clade, gene) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = gene, values_from = N)
}

#### count sample number ====
PM_all <- FQR_PM_data %>% 
  mutate(all = "all") %>% 
  group_by(all) %>% 
  summarise(N = n()) %>% 
  left_join(PM_summarise_all(x = summary_category[1])) %>% 
  left_join(PM_summarise_all(x = summary_category[2])) %>% 
  left_join(PM_summarise_all(x = summary_category[3])) %>% 
  left_join(PM_summarise_all(x = summary_category[4])) %>% 
  left_join(PM_summarise_all(x = summary_category[5])) %>% 
  left_join(PM_summarise_all(x = summary_category[6])) %>% 
  left_join(PM_summarise_all(x = summary_category[7])) %>% 
  left_join(PM_summarise_all(x = summary_category[8])) %>% 
  left_join(PM_summarise_all(x = summary_category[9])) %>% 
  left_join(PM_summarise_all(x = summary_category[10])) %>% 
  left_join(PM_summarise_all(x = summary_category[11])) %>% 
  left_join(PM_summarise_all(x = summary_category[12])) %>% 
  left_join(PM_summarise_all(x = summary_category[13])) %>% 
 transpose(keep.names = "mutation", make.names = "all")
  
PM_clade <- sample_info %>% 
  group_by(clade) %>% 
  summarise(N = n()) %>% 
  left_join(PM_summarise_clade(x = summary_category[1])) %>% 
  left_join(PM_summarise_clade(x = summary_category[2])) %>% 
  left_join(PM_summarise_clade(x = summary_category[3])) %>% 
  left_join(PM_summarise_clade(x = summary_category[4])) %>% 
  left_join(PM_summarise_clade(x = summary_category[5])) %>% 
  left_join(PM_summarise_clade(x = summary_category[6])) %>% 
  left_join(PM_summarise_clade(x = summary_category[7])) %>% 
  left_join(PM_summarise_clade(x = summary_category[8])) %>% 
  left_join(PM_summarise_clade(x = summary_category[9])) %>% 
  left_join(PM_summarise_clade(x = summary_category[10])) %>% 
  left_join(PM_summarise_clade(x = summary_category[11])) %>% 
  left_join(PM_summarise_clade(x = summary_category[12])) %>% 
  left_join(PM_summarise_clade(x = summary_category[13])) %>% 
  transpose(keep.names = "mutation", make.names = "clade")

#### output data ==== 
PM_summary <- left_join(PM_all, PM_clade) %>% 
  mutate(gene = case_when(
    str_detect(mutation, c("Gly81Asp|Ser83Ile|Ser83Arg|Asp87Tyr|Asp87Gly|Asp87Asn|gyrA")) ~ "gyrA",
    str_detect(mutation, c("Ser80Ile|Glu84Gly|Ser80Arg|Ala81Pro|His75Gln|Glu84Lys|Ala108Thr|parC")) ~ "parC",
    str_detect(mutation, c("Ser458Ala|Ile444Phe|Ser458Trp|parE")) ~ "parE",
    TRUE ~ "all"),
    position = as.numeric(str_extract(mutation, pattern = "[[:digit:]]+"))) %>% 
  mutate(gene = 
           ifelse(mutation == "gyrA + parC" | 
                    mutation == "gyrA + parE", "QRDR_multi",
                  ifelse(mutation == "Ser83Ile + Asp87Tyr" | 
                           mutation == "Ser83Ile + Asp87Asn" | 
                           mutation == "Ser83Ile + Asp87Gly" |
                           mutation == "multi_gyrA",  "gyrA_multi",
                         ifelse(mutation == "quinolone" | 
                                  mutation == "w/o FQR mutation", "quinolone", gene)))) %>% 
  mutate(position = replace_na(position, 0)) %>% 
  select(gene, position, mutation, all, as.character(1:14)) %>% 
  group_by(gene) %>% 
  arrange(position, -all, .by_group = TRUE) %>% 
  ungroup() %>% 
  select(-position)

### SNP10 cluster ====
#### function ====
PM_SNP_summarise_all <- function(x){
  result <- FQR_PM_data %>% 
    mutate(all = "all") %>% 
    select(all, cluster, gene = all_of(x)) %>% 
    filter(gene != "-") %>% 
    group_by(all, cluster, gene) %>% 
    summarise(N = n()) %>% 
    group_by(all, gene) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = gene, values_from = N)
}
PM_SNP_summarise_clade <- function(x){
  result <- FQR_PM_data %>% 
    select(clade, cluster, gene = all_of(x)) %>% 
    filter(gene != "-") %>% 
    group_by(clade, cluster, gene) %>% 
    summarise(N = n()) %>% 
    group_by(clade, gene) %>% 
    summarise(N = n()) %>% 
    pivot_wider(names_from = gene, values_from = N)
}

#### count sample number ====
PM_all_SNP <- FQR_PM_data %>% 
  mutate(all = "all") %>% 
  group_by(all, cluster) %>% 
  summarise(N = n()) %>% 
  group_by(all) %>% 
  summarise(N = n()) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[1])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[2])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[3])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[4])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[5])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[6])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[7])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[8])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[9])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[10])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[11])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[12])) %>% 
  left_join(PM_SNP_summarise_all(x = summary_category[13])) %>% 
  transpose(keep.names = "mutation", make.names = "all")
  
PM_clade_SNP <- sample_info %>% 
  group_by(clade, cluster) %>% 
  summarise(N = n()) %>% 
  group_by(clade) %>% 
  summarise(N = n()) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[1])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[2])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[3])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[4])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[5])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[6])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[7])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[8])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[9])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[10])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[11])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[12])) %>% 
  left_join(PM_SNP_summarise_clade(x = summary_category[13])) %>% 
  transpose(keep.names = "mutation", make.names = "clade")

#### output data ==== 
PM_summary_SNP <- left_join(PM_all_SNP, PM_clade_SNP) %>% 
  mutate(gene = case_when(
    str_detect(mutation, c("Gly81Asp|Ser83Ile|Ser83Arg|Asp87Tyr|Asp87Gly|Asp87Asn|gyrA")) ~ "gyrA",
    str_detect(mutation, c("Ser80Ile|Glu84Gly|Ser80Arg|Ala81Pro|His75Gln|Glu84Lys|Ala108Thr|parC")) ~ "parC",
    str_detect(mutation, c("Ser458Ala|Ile444Phe|Ser458Trp|parE")) ~ "parE",
    TRUE ~ "all"),
    position = as.numeric(str_extract(mutation, pattern = "[[:digit:]]+"))) %>% 
  mutate(gene = 
           ifelse(mutation == "gyrA + parC" | 
                    mutation == "gyrA + parE", "QRDR_multi",
                  ifelse(mutation == "Ser83Ile + Asp87Tyr" | 
                           mutation == "Ser83Ile + Asp87Asn" | 
                           mutation == "Ser83Ile + Asp87Gly" |
                           mutation == "multi_gyrA",  "gyrA_multi",
                         ifelse(mutation == "quinolone" | 
                                  mutation == "w/o FQR mutation", "quinolone", gene)))) %>% 
  mutate(position = replace_na(position, 0)) %>% 
  select(gene, position, mutation, all, as.character(1:14)) %>% 
  group_by(gene) %>% 
  arrange(position, -all, .by_group = TRUE) %>% 
  ungroup() %>% 
  select(-c(gene, position))

### combine strain result and SNP10 cluster result ====
PM_summary_output <- PM_summary %>% 
  left_join(x = PM_summary_SNP, y = ., by = "mutation") %>% 
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
  select(gene, mutation, all, as.character(1:14)) %>% 
  mutate_all(~ifelse(. == "0 (0)", "", .)) %>% 
  write_tsv(file = "RESULTS/AMR/AMR_summary_table_FQR.tsv")
  