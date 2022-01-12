# clade core gene list

# library ====
library(tidyverse)
#library(RVAideMemoire)

# output ====
output_path <- "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/core_gene/"

# data manipulate ====
## for scoary ====
sample_info <- read_tsv(file = "DATA/final_dataset/sample_info.tsv") %>% 
  filter(`Representative strain in each SNP10 cluster` == "Rep") %>% 
  select(strainID = `Strain ID in this paper`,
         ANI = `ANI clade`) %>% 
  mutate(clade = case_when(ANI == 1 | ANI == 2 ~ 1,
                           TRUE ~ 0)) %>% 
  select(Name = strainID, clade, ANI) %>% 
  write_csv(x = select(sample_info, Name, clade), file = str_c(output_path, "strain_clade_info.csv"), na = "")

pangenome_info <- read_csv(file = "DATA/pangenome_analysis/gene_presence_absence.csv") %>% 
  select(c("Gene","Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", 
            "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", 
            "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"))

pangenome_preab <- read_csv(file = "DATA/pangenome_analysis/gene_presence_absence.csv") %>% 
  select(-c("Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", 
            "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", 
            "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc")) %>% 
  pivot_longer(cols = -Gene, names_to = "strain", values_to = "locusID") %>% 
  mutate(strain = str_split(string = strain, pattern = "_", simplify = TRUE)[,1]) %>% 
  filter(strain %in% sample_info$Name) %>% 
  mutate(pre_ab = case_when(is.na(locusID) ~ 0, TRUE ~ 1)) %>% 
  group_by(Gene) %>% 
  filter(sum(pre_ab) > 0) %>% 
  select(-pre_ab) %>% 
  pivot_wider(names_from = "strain", values_from = "locusID") %>% 
  right_join(x = pangenome_info, y = ., by = "Gene") %>% 
  write_csv(file = str_c(output_path, "strain_rename_gene_presence_absence.csv"), na = "")
  select(-c("Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", 
            "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", 
            "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc")) %>% 
  pivot_longer(cols = -Gene, names_to = "strain", values_to = "pre_ab") %>% 
  left_join(x = select(sample_info, Name, ANI), y = ., by = c("Name" = "strain")) %>% 
  mutate(pre_ab = case_when(is.na(pre_ab) ~ 0, TRUE ~ 1))

# gene filtering ====
## get accessory gene ====
accessory_gene <- pangenome_preab %>% 
  group_by(Gene) %>% 
  nest() %>% 
  mutate(presence_ratio = map_dbl(data, function(x){
    sum(x$pre_ab) / nrow(x)
  })) %>% 
  filter(presence_ratio < .99) %>% 
  select(-presence_ratio) %>% 
  unnest(cols = c(data))

# gene presence absence in SNP10 strain ====
all_gene_preab <- accessory_gene %>% 
  group_by(Gene) %>% 
  summarise(ANI = 15, # 15 = all
            strain_num = n(),
            n_presence = sum(pre_ab),
            n_absence = strain_num - n_presence,
            ratio_presence = n_presence / strain_num * 100) %>% 
  ungroup() %>% 
  select(Gene, clade = ANI, presence = n_presence, absence = n_absence, ratio_presence)
  
# gene presence absence each clade data ====
clade_gene_preab <- accessory_gene %>% 
  group_by(Gene, ANI) %>% 
  summarise(strain_num = n(),
            n_presence = sum(pre_ab),
            n_absence = strain_num - n_presence,
            ratio_presence = n_presence / strain_num * 100) %>% 
  ungroup() %>% 
  select(Gene, clade = ANI, presence = n_presence, absence = n_absence, ratio_presence)

# gene presence absence in clades except for clade 1 and 2 ====
other_gene_preab <- accessory_gene %>% 
  filter(ANI != 1 & ANI != 2) %>% 
  group_by(Gene) %>% 
  summarise(ANI = 16, # 16 = clade 3-14
            strain_num = n(),
            n_presence = sum(pre_ab),
            n_absence = strain_num - n_presence,
            ratio_presence = n_presence / strain_num * 100) %>% 
  ungroup() %>% 
  select(Gene, clade = ANI, presence = n_presence, absence = n_absence, ratio_presence)

# get clade concentrated gene list ====
## clade 1 ====
conc_clade1_list <- clade_gene_preab %>% 
  filter(clade == 1, ratio_presence > 90) %>% 
  left_join(all_gene_preab, by = "Gene", suffix = c("_clade", "_all")) %>% 
  filter(ratio_presence_clade > ratio_presence_all) %>% 
  pull(Gene)

## clade 2 ====
conc_clade2_list <- clade_gene_preab %>% 
  filter(clade == 2, ratio_presence > 90) %>% 
  left_join(all_gene_preab, by = "Gene", suffix = c("_clade", "_all")) %>% 
  filter(ratio_presence_clade > ratio_presence_all) %>% 
  pull(Gene)

## compare gene list ====
### all gene ====
conc_all_clade1and2_list <- union(conc_clade1_list, conc_clade2_list)

### share ====
conc_clade1and2_list <- intersect(conc_clade1_list, conc_clade2_list)

### clade 1 only ====
conc_clade1only_list <- setdiff(conc_clade1_list, conc_clade2_list)

### clade 2 only ==== 
conc_clade2only_list <- setdiff(conc_clade2_list, conc_clade1_list)

## import scoary result ====
scoary_result <- read_csv(file = str_c(output_path, "clade_29_12_2021_0108.results.csv"), col_names = TRUE) %>% 
  select(Gene, `Non-unique Gene name`, Annotation, Naive_p, Bonferroni_p, Benjamini_H_p)

## import SM39 and Db11 annotation ID data ====
SM39_locus_ID <- read_tsv(file = "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/core_gene/add_annotation/DL536_SM39_convert.tsv")
Db11_locus_ID <- read_tsv(file = "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/core_gene/add_annotation/DL020_Db11_convert.tsv")

### pangenome data to annotation data ====
SM39_Db11_annotation <- read_csv(file = "DATA/pangenome_analysis/gene_presence_absence.csv") %>% 
  select(-c("Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", 
            "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", 
            "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc")) %>% 
  pivot_longer(cols = -Gene, names_to = "strain", values_to = "locusID") %>% 
  mutate(strain = str_split(string = strain, pattern = "_", simplify = TRUE)[,1]) %>% 
  filter(strain == "DL536" | strain == "DL020") %>% 
  filter(Gene %in% conc_clade1and2_list) %>% 
  mutate(locusID = str_split(string = locusID, pattern = "___", simplify = TRUE)[,1]) %>% 
  pivot_wider(names_from = "strain", values_from = "locusID") %>% 
  left_join(Db11_locus_ID, by = c("DL020" = "roary_locus_tag")) %>% 
  left_join(SM39_locus_ID, by = c("DL536" = "roary_locus_tag"), suffix = c("_Db11", "_SM39")) %>% 
  mutate(`Annotation in the SM39/Db11 genomes` = case_when(!is.na(genbank_locus_tag_SM39) ~ genbank_product_SM39,
                                                           TRUE ~ genbank_product_Db11))

### no product data ====
SM39_Db11_annotation_noDATA <- SM39_Db11_annotation %>% 
  filter(is.na(genbank_locus_tag_SM39)) %>% 
  arrange(Gene)

output_clade_preab <- clade_gene_preab %>% 
  filter(Gene %in% conc_clade1and2_list) %>% 
  filter(clade == 1 | clade == 2) %>% 
  group_by(Gene) %>% 
  mutate(`clade 1 and clade 2` = sum(presence) / sum(presence, absence) * 100) %>% 
  ungroup() %>% 
  select(-c(presence, absence)) %>% 
  pivot_wider(names_from = "clade", values_from = "ratio_presence", names_prefix = "clade ") %>% 
  left_join(other_gene_preab, by = "Gene") %>% 
  left_join(all_gene_preab, by = "Gene", suffix = c("_other", "_all")) %>% 
  select(Gene, `clade 1 and clade 2`, `clade 1`, `clade 2`, 
         `others (not-clade 1/2)` = ratio_presence_other) %>% 
  left_join(scoary_result, by = "Gene") %>% 
  left_join(SM39_Db11_annotation, by = "Gene") %>% 
  mutate(Db11 = str_replace_na(genbank_locus_tag_Db11))
  select(Gene, `Non-unique Gene name`, Annotation, `p value*` = Bonferroni_p,
         `clade 1`, `clade 2`, `clade 1 and clade 2`, `others (not-clade 1/2)`,
         Db11_order = genbank_locus_tag_Db11, SM39_order = genbank_locus_tag_SM39,
         `Annotation in the SM39/Db11 genomes`) %>% 
  mutate(`Annotation in the SM39/Db11 genomes` = URLdecode(`Annotation in the SM39/Db11 genomes`)) %>% 
  filter(`p value*` < 0.01) %>% 
  arrange(`p value*`) %>% 
  write_tsv(file = str_c(output_path, "scoary_1and2_vs_others.tsv"), na = "", )
  