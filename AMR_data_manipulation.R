# script option ====
## library ====
library(tidyverse)
library(data.table)

## function ====
get_cat_data <- function(data, MDR){
  resist_col_name <- MDR
  cat_data <- data %>% 
    filter(for_MDR_dicision == MDR) %>% 
    select(-c(class, category, for_MDR_dicision)) %>% 
    transpose(keep.names = "sample_ID", make.names = "gene name") %>% 
    pivot_longer(cols = -sample_ID, names_to = "gene", values_to = "PA") %>% 
    filter(PA != "-") %>% 
    group_by(sample_ID) %>% 
    summarise(!!resist_col_name := "resist")
}

# data import ====
AMR <- read_tsv(file = "DATA/final_dataset/AMR_info.tsv")
  
# MDR decision ====
## make each categorical table ====
# category = "Agly", "Qnr", "Phe", "ST-resi", "Tet", "carBla", "ESBL"
Agly <- get_cat_data(data = AMR, MDR = "Agly")
Qnr <- get_cat_data(data = AMR, MDR = "Qnr")
Phe <- get_cat_data(data = AMR, MDR = "Phe")
STresi <- get_cat_data(data = AMR, MDR = "ST-resi")
Tet <- get_cat_data(data = AMR, MDR = "Tet")
carBla <- get_cat_data(data = AMR, MDR = "carBla")
ESBL <- get_cat_data(data = AMR, MDR = "ESBL")

## combine all data ====
resist_distribution <- Agly %>% 
  full_join(Qnr) %>%  full_join(Phe) %>% 
  full_join(STresi) %>%  full_join(Tet) %>% 
  full_join(carBla) %>%  full_join(ESBL) %>% 
  arrange(sample_ID)

## count resist category number ====
resist_num <- resist_distribution %>% 
  pivot_longer(cols = -sample_ID, names_to = "category", values_to = "resistance") %>% 
  filter(!is.na(resistance)) %>% 
  group_by(sample_ID) %>% 
  summarise(resistance_num = n()) %>% 
  ungroup()

## all sample resistance data ====
MDR_distribution <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  select(sample_ID = "Strain ID in this paper") %>% 
  left_join(resist_distribution) %>% 
  left_join(resist_num) %>% 
  write_tsv(file = "DATA/final_dataset/MDR_info.tsv")
