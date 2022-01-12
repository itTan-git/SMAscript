# clustering based on pangenome presence/absence data

# library ====
library(tidyverse)
library(gplots)
library(ape)

# file path ====
file_path <- "/Users/taniguchi/git/Ono_data/DATA/pangenome_analysis/"

# data import ====
data_ori <- read_tsv(file = str_c(file_path, "sample_data_gene_presence_absence.tsv", sep = ""), col_names = T, na = "")
data_ori[is.na(data_ori)] <- 0
cud_color <- c( "#dc143c", "#ffa07a", "#bd7093", "#ffb6c1",
                "#ff00ff", "#9400d3", "#ffff00", "#00ff7f",
                "#00ffff", "#00bfff", "#ff8c00", "#a9a9a9")

data <- data_ori[,-c(2:4,7:11)] 
data[,4:ncol(data)][data[,4:ncol(data)]!=0] <- 1 
data[,4:ncol(data)] <- mutate_if(data[,4:ncol(data)], is.character, as.numeric) 
data_use <- data

## manipulate data for clustering ====
data_use_for_hc <- as.data.frame(data_use[,-c(1:3)]) 
rownames(data_use_for_hc) <- data_use$ID
hc_data_use <- hclust(d = dist(x = data_use_for_hc, method = "euclidean"), method = "ward.D2")

phy_data_use <- as.phylo(hc_data_use)
write.tree(phy = phy_data_use, file = str_c(file_path, "sample_pangenome_clustering.nwk", sep = ""))
