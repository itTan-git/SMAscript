# convert annotation file used in Roary to genbank

# library ====
library(tidyverse)
library(bedtoolsr)

# file path ====
file_path <- "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/core_gene/add_annotation/"

# SM39 ====
## data import ====
SM39_genbank <- read_tsv(file = str_c(file_path, "DL536_SM39_Genbank.gff"), comment = "#",
                       col_names = c("seqname", "source", "feature", "start", 
                                     "end", "score", "strand", "frame", "attribute")) %>% 
  mutate(attribute = str_replace_all(string = attribute, pattern = "'", replacement = "%27"),
         seqname = case_when(seqname == "AP013063.1" ~ "sequence1", 
                             seqname == "AP013064.1" ~ "sequence2",
                             TRUE ~ "sequence3")) %>% 
  filter(str_detect(string = attribute, pattern = "locus_tag=")) %>% 
  filter(feature %in% c("CDS", "tRNA", "rRNA"))

SM39_roary <- read_tsv(file = str_c(file_path, "DL536_GCA_000828775.fasta_dfast.gff"), n_max = 4997, comment = "#",
                       col_names = c("seqname", "source", "feature", "start", 
                                     "end", "score", "strand", "frame", "attribute")) %>% 
  mutate(attribute = str_replace_all(string = attribute, pattern = "'", replacement = "%27")) %>% 
  filter(str_detect(string = attribute, pattern = "locus_tag="))
  

## intersect 2 files ====
SM39_intersect <- bt.intersect(a = SM39_roary, b = SM39_genbank, wa = TRUE, wb = TRUE) %>% 
  `colnames<-`(c(str_c("roary", colnames(SM39_roary), sep = "_"), str_c("genbank", colnames(SM39_genbank), sep = "_"))) %>% 
  mutate(roary_locus_tag = str_extract(string = roary_attribute, pattern = "(?<=locus_tag=).+?(?=;)"),
         genbank_locus_tag = str_extract(string = genbank_attribute, pattern = "(?<=locus_tag=).+"),
         genbank_locus_tag = str_remove(string = genbank_locus_tag, pattern = ";.+")) %>% 
  group_by(roary_locus_tag) %>% 
  mutate(diff_start = abs(roary_start - genbank_start),
         diff_end = abs(roary_end - genbank_end)) %>% 
  nest() %>% 
  mutate(genbank_locus_tag = map_chr(data, function(x){
    min_diff_start <- x %>% 
      filter(diff_start == min(diff_start))
    min_diff_end <- x %>% 
      filter(diff_end == min(diff_end))
    region_start <- min(min_diff_start$genbank_start, min_diff_start$genbank_end,
                        min_diff_end$genbank_start, min_diff_end$genbank_end)
    region_end <- max(min_diff_start$genbank_start, min_diff_start$genbank_end,
                        min_diff_end$genbank_start, min_diff_end$genbank_end)
    locus_tag <- x %>%  
      filter(genbank_start >= region_start, genbank_end <= region_end) %>% 
      arrange(genbank_locus_tag) %>% 
      pull(genbank_locus_tag)
    return(str_c(locus_tag, collapse = ";"))
    })) %>% 
  select(roary_locus_tag, genbank_locus_tag)

## add annotation from genbank file ====
SM39_genbank_annotation <- SM39_genbank %>% 
  mutate(genbank_locus_tag = str_extract(string = attribute, pattern = "(?<=locus_tag=).+"),
         genbank_locus_tag = str_remove(string = genbank_locus_tag, pattern = ";.+"),
         genbank_product = case_when(
           feature == "CDS" ~ str_extract(string = attribute, pattern = "(?<=product=).+?(?=;)"),
           TRUE ~ str_extract(string = attribute, pattern = "(?<=product=).+")),
         genbank_product = case_when(is.na(genbank_product) ~ str_extract(string = attribute, pattern = "(?<=Note=).+?(?=;)"),
                                     TRUE ~ genbank_product))

## output converted annotation ID ====
SM39_intersect_annotation <- SM39_intersect %>% 
  left_join(SM39_genbank_annotation, by = "genbank_locus_tag") %>% 
  select(roary_locus_tag, genbank_locus_tag, genbank_product) %>% 
  write_tsv(file = str_c(file_path, "DL536_SM39_convert.tsv"), na = "")

# Db11 ====
## data import ====
Db11_genbank <- read_tsv(file = str_c(file_path, "DL020_Db11_Genbank.gff"), comment = "#",
                       col_names = c("seqname", "source", "feature", "start", 
                                     "end", "score", "strand", "frame", "attribute")) %>% 
  mutate(attribute = str_replace_all(string = attribute, pattern = "'", replacement = "%27"),
         seqname = "sequence1") %>% 
  filter(str_detect(string = attribute, pattern = "locus_tag=")) %>% 
  filter(feature %in% c("CDS", "tRNA", "rRNA"))

Db11_roary <- read_tsv(file = str_c(file_path, "DL020_GCA_000513215.fasta_dfast.gff"), n_max = 4787, comment = "#",
                       col_names = c("seqname", "source", "feature", "start", 
                                     "end", "score", "strand", "frame", "attribute")) %>% 
  mutate(attribute = str_replace_all(string = attribute, pattern = "'", replacement = "%27")) %>% 
  filter(str_detect(string = attribute, pattern = "locus_tag="))

## intersect 2 files ====
Db11_intersect <- bt.intersect(a = Db11_roary, b = Db11_genbank, wa = TRUE, wb = TRUE) %>% 
  `colnames<-`(c(str_c("roary", colnames(Db11_roary), sep = "_"), str_c("genbank", colnames(Db11_genbank), sep = "_"))) %>% 
  mutate(roary_locus_tag = str_extract(string = roary_attribute, pattern = "(?<=locus_tag=).+?(?=;)"),
         genbank_locus_tag = str_extract(string = genbank_attribute, pattern = "(?<=locus_tag=).+"),
         genbank_locus_tag = str_remove(string = genbank_locus_tag, pattern = ";.+")) %>% 
  group_by(roary_locus_tag) %>%
  mutate(diff_start = abs(roary_start - genbank_start),
         diff_end = abs(roary_end - genbank_end)) %>% 
  nest() %>% 
  mutate(genbank_locus_tag = map_chr(data, function(x){
    min_diff_start <- x %>% 
      filter(diff_start == min(diff_start))
    min_diff_end <- x %>% 
      filter(diff_end == min(diff_end))
    region_start <- min(min_diff_start$genbank_start, min_diff_start$genbank_end,
                        min_diff_end$genbank_start, min_diff_end$genbank_end)
    region_end <- max(min_diff_start$genbank_start, min_diff_start$genbank_end,
                        min_diff_end$genbank_start, min_diff_end$genbank_end)
    locus_tag <- x %>%  
      filter(genbank_start >= region_start, genbank_end <= region_end) %>% 
      arrange(genbank_locus_tag) %>% 
      pull(genbank_locus_tag)
    return(str_c(locus_tag, collapse = ";"))
    })) %>% 
  select(roary_locus_tag, genbank_locus_tag)

## add annotation from genbank file ====
Db11_genbank_annotation <- Db11_genbank %>% 
  mutate(genbank_locus_tag = str_extract(string = attribute, pattern = "(?<=locus_tag=).+"),
         genbank_locus_tag = str_remove(string = genbank_locus_tag, pattern = ";.+"), 
         genbank_product = str_extract(string = attribute, pattern = "(?<=product=).+?(?=;)"))

## output converted annotation ID ====
Db11_intersect_annotation <- Db11_intersect %>% 
  left_join(Db11_genbank_annotation, by = "genbank_locus_tag") %>% 
  select(roary_locus_tag, genbank_locus_tag, genbank_product) %>% 
  write_tsv(file = str_c(file_path, "DL020_Db11_convert.tsv"), na = "")
