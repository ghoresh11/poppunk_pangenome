library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)


setwd("/Users/gh11/poppunk_pangenome/4_pairwise_roary/classify_genes/")

all_cogs = read.table("cog_per_gene.csv", sep = ",", header = T, comment.char = "", quote = "", stringsAsFactors = F)
cols = read.table(file = "cog_descs.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)

## add COGs that are mixed...
missing_cogs = unique(all_cogs$COGs)
missing_cogs = missing_cogs[which(!check_missing_cogs %in% cols$COG)]
cols = rbind(cols,
             data.frame(COG = missing_cogs,
                        Desc = rep("Mixed", length(missing_cogs)),
                        Title = rep("Mixed", length(missing_cogs)),
                        Col = rep("black", length(missing_cogs)), stringsAsFactors = F))
all_cogs$COGs = factor(all_cogs$COGs, cols$COG)
o = c("Real_core", "Ubiquitous", "Missing_in_one",
      "40_45_Cluster_specific", "25_39_Cluster_specific", "10_24_Cluster_specific",
      "2_9_Cluster_specific",  "Cluster_specific" ,
      "40_46_varied","25_39_varied","10_24_varied" , "2_9_varied","Intermediate_and_rare",
      "Multicluster_intermediate","Cluster_specific_intermediate",
      "Multicluster_rare","Cluster_specific_rare" )
all_cogs$Cat = factor(all_cogs$Cat, o)

## add facets based on general distribution
all_cogs$Group = factor(all_cogs$Group, c("Core", "Cluster specific core", "Varied", "Intermediate","Rare"))

## COG category accross the gene classes
ggplot(all_cogs, aes(x = Cat, fill = COGs)) + geom_bar(color = "black", position = "fill") +
  scale_fill_manual(values = cols$Col) + theme_classic(base_size = 14)+ 
  facet_grid(. ~ Group, space = "free_x", switch = "x", scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


