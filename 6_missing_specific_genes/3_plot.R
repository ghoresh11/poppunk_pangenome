library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)


setwd("/Users/gh11/poppunk_pangenome/6_missing_specific_genes/")


#### MAIN ####

## 1. Using analysis on merged genes to remove wrong ones
specific = read.csv("specific_genes_detailed.csv", sep ="\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
missing = read.csv("missing_genes_detailed.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")

details = read.table("specific_gene_types.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F)
specific_details = details[match(specific$gene, details$Specific_gene),]
specific = cbind(specific, specific_details)
specific$Type[which(specific$Type == "truly cluster specific" & specific$classification == "unknown")] = "unknown function"

phylogroup_order =  c("B2", "F", "D","E","A","C","B1","U") ## define
type_order = c("wrong","truncated","longer","same functional group","truly cluster specific", "unknown function")
cols = c("#332288","#CC6677","#DDCC77","#88CCEE","#44AA99", "#d3d3d3")

specific$cluster = as.character(specific$cluster)
specific$Type = factor(specific$Type, type_order)
specific$phylogroup = factor(specific$phylogroup,phylogroup_order)
ggplot(specific, aes(x = cluster, fill = Type)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols)


## do the same for missing
missing_details = details[match(missing$gene, details$Missing_gene),]
missing = cbind(missing, missing_details)
missing$Type[which(is.na(missing$Type))] = "truly missing in cluster"
missing$Type[which(missing$Type == "longer")] = "truncated*"
missing$Type[which(missing$Type == "truncated")] = "longer*"
missing$Type[which(missing$Type == "truly missing in cluster" & missing$classification == "unknown")] = "unknown function"

type_order = c("wrong","truncated*","longer*","same functional group","truly missing in cluster", "unknown function")
missing$cluster = as.character(missing$cluster)
missing$Type = factor(missing$Type, type_order)
missing$phylogroup = factor(missing$phylogroup,phylogroup_order)
for (cluster in unique(specific$cluster)) { ## add clusters with no missing genes
  if (!cluster %in% missing$cluster){
    add = missing[1,]
    add$cluster[1] = cluster
    add$num_genes[1] = 0
    add$phylogroup[1] = specific$phylogroup[specific$cluster == cluster][1]
    for (i in 4:length(add)) {
      add[1,i] = NA
    }
    missing = rbind(missing, add)
  }
}
ggplot(missing, aes(x = cluster, fill = Type, color = Type)) + geom_bar() + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols,
                    labels = c("wrong","truncated","longer","same functional group","truly missing in cluster","unknown function","")) + 
  scale_color_manual(values = c(rep("black", 6),NA), guide = F)


functional_groups = read.table("functional_clusters.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F, quote = "")

## specific:
specific_functional_groups = functional_groups[which(functional_groups$Type %in% c("both","specific") & functional_groups$Num_clusters > 1),]

 

