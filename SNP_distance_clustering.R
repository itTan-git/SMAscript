# SNP distance clustering
# SNP dist = 5, 10, 20, 30

## library ====
library(tidyverse)
library(ggtree)

## import SNP distance data ====
snpdist <- read.table(file = "DATA/snp_distance_clustering/roary_CoreGeneAllignment_snpdist.tsv", 
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
sample_list <- rownames(snpdist) %>% 
  as_tibble() %>%  arrange() %>% `names<-`("sample")

## clustering based on SNP distance and get dendrogram ====
clust <- hclust(as.dist(snpdist), method = "ward.D2")
tree <- ggtree(tr = clust)

## get max SNP distance within node ==== 
### make node list ====
node_list <- tree$data %>% 
  filter(!isTip) %>%  
  select(node) 
node_list$max_snpdist <- NA 

### check max SNP distance ====
for (node in 1:nrow(node_list)){
  temp_clade <- tree %>% 
    groupClade(.node = node_list$node[node]) 
  temp_clade_tip <- temp_clade$data %>% 
    filter(group == 1 & isTip == T) %>% 
    pull(label) 
  temp_clade_ANI <- snpdist[is.element(rownames(snpdist), temp_clade_tip), 
                            is.element(rownames(snpdist), temp_clade_tip), 
                            drop = F] 
  node_list$max_snpdist[node] <- max(temp_clade_ANI) 
}

## get node list which SNP distance is lower than threshold ====
node_list_threshold <- node_list %>% 
  mutate(
    SNP5 = case_when(max_snpdist <= 5 ~ node),
    SNP10 = case_when(max_snpdist <= 10 ~ node),
    SNP20 = case_when(max_snpdist <= 20 ~ node),
    SNP30 = case_when(max_snpdist <= 30 ~ node)
    )

result <- sample_list %>% 
  mutate(strainID = str_split(string = sample, pattern = "_", simplify = TRUE)[, 1])
  
for (snp in 1:4){
  temp_node_list <- node_list_threshold %>% 
    select(1, snp+2) %>% 
    `colnames<-`(c("node", "SNP")) %>% 
    filter(!is.na(SNP))
    
  ## add node ID to strain ====
  snpdist_node <- data.frame(label = NULL, node = NULL) 
  for (node in 1:nrow(temp_node_list)){
    temp_clade <- tree %>% 
      groupClade(.node = temp_node_list$node[node]) 
    temp_clade_tip <- temp_clade$data %>% 
      filter(group == 1 & isTip == T) %>% 
      select(label) %>%  
      mutate(node = temp_node_list$node[node]) 
    snpdist_node <- rbind(snpdist_node, temp_clade_tip) 
  }
  
  ## add minimum node ID to strain ====
  snpdist_cluster <- snpdist_node %>% 
    group_by(label) %>% 
    summarise(cluster = min(node)) %>% 
    mutate(cluster = as.character(factor(cluster, labels = 1:length(unique(cluster))))) 
  
  sample_cluster <- sample_list %>% 
    left_join(snpdist_cluster, by = c("sample" = "label")) %>% 
    arrange(sample)
  
  temp_sample_cluster <- sample_cluster %>% 
    mutate(clusterID = case_when(!is.na(cluster) ~ str_c("c", str_pad(string = cluster, width = 3, side = "left", pad = "0")),
                                  TRUE ~ str_split(string = sample, pattern = "_", simplify = TRUE)[, 1])) %>% 
    select(-cluster)
  
  result <- result %>% 
    left_join(temp_sample_cluster, by = "sample")
}

## output ====
result_out <- result %>% 
  select(-sample) %>% 
  `colnames<-`(c("strainID", "clusterID_SNP5", "clusterID_SNP10", "clusterID_SNP20", "clusterID_SNP30")) %>% 
  write_tsv(file = "LETTER/submission/Microbiology_Genomics/response_to_reviewers/additional_analysis/SNPclustering/SNP_dist_clustering.tsv", na = "")



