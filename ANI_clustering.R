# clustering based all-to-all ANI

# library
library(tidyverse)
library(ggtree)

# ANI clustering dendrogram data ====
## import matrix data ====
ono_allvsall <- read_tsv(file = "DATA/ANIm_clustering/ono_tree_matrix.tsv") 
ono_allvsall_df <- ono_allvsall %>% 
  select(-allvsall) %>%  data.frame() %>% 
  `rownames<-`(ono_allvsall$allvsall)
sampleID <- str_split(ono_allvsall$allvsall, "_", simplify = T)
rownames(ono_allvsall_df) <- sampleID[,1] ; colnames(ono_allvsall_df) <- sampleID[,1] 
ono_tree_matrix <- as.matrix(ono_allvsall_df) 
ono_tree_dist <- dist(ono_tree_matrix, method = "euclidean") 

### draw tree complete method(default) ====
ono_tree_hclust <- ape::as.phylo(hclust(ono_tree_dist, method = "complete"))
write.tree(phy = ono_tree_hclust, file = "RESULTS/ANIm_clustering/ANI_clustering_tree.nwk")

# minimum ANI within node ====
## draw dendrogram ====
tree <- ono_tree_hclust
gg_tree <- ggtree(tr = tree, ladderize = T, layout = "rectangular") +
  geom_tiplab(align = F, cex = .5, angle = 270) + 
  layout_dendrogram() 

## calculate minimum ANI in each node ====
### make node list ====
node_list <- gg_tree$data %>% 
  filter(!isTip) %>%  
  select(node) 

minANI <- NULL
### calculation ====
for (N in 1:length(node_list$node)){
  temp_clade <- gg_tree %>% 
    groupClade(.node = node_list$node[N]) 
  temp_clade_tip <- temp_clade$data %>% 
    filter(group == 1 & isTip == T) %>% 
    pull(label) 
  temp_clade_ANI <- ono_allvsall_df[is.element(rownames(ono_allvsall_df), temp_clade_tip), 
                                    is.element(rownames(ono_allvsall_df), temp_clade_tip), 
                                    drop = F] 
  minANI <- c(minANI, min(temp_clade_ANI))
}
node_list <- mutate(node_list, minimum_ANI = minANI) %>% 
  write_tsv(file = "RESULTS/ANIm_clustering/minimumANI_eachnode.tsv")

# add clade info ====
## filtering node for display ====
threshold_node <- node_list %>% 
  filter(minimum_ANI < 97.5) %>% 
  pull(node)
clade_node <- NULL
### get node list which have "clade_node" in parent ====
for (i in 1:length(threshold_node)){
  temp_node <- gg_tree$data %>% 
    filter(parent == threshold_node[i] & node != threshold_node[i]) %>% 
    pull(node)
  clade_node <- c(clade_node, temp_node)
}
clade_node <- sort(unique(clade_node))

## make group with clade_node ====
clade_group_tree <- gg_tree %>% 
  groupClade(clade_node) %>% 
  ggtree::rotate(776) %>% ggtree::rotate(781) %>% ggtree::rotate(788) %>% ggtree::rotate(785) 
clade_info <- clade_group_tree$data %>% 
  filter(isTip) %>% 
  arrange(y) %>% 
  mutate(group = as.numeric(group)) %>% 
  mutate(group = factor(group, levels = unique(group), labels = c(1:length(unique(group))))) %>% 
  select(sample_ID = label, clade = group) %>% 
  arrange(sampleID) %>%  
  write_tsv(file = "RESULTS/ANIm_clustering/ANIclade_info.tsv")
