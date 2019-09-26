library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(ape)
library(ggtree)
library(dplyr)

setwd("/Users/gh11/poppunk_pangenome/6_missing_specific_genes/")


#### MAIN ####

phylogroup_order =  c("B2", "F", "D","E","A","C","B1","U") ## define
type_order = c("wrong", "truncated","longer","same functional group","true")
cols = c("black", "#CC6677","#DDCC77","#88CCEE","#44AA99")

## 1. Using analysis on merged genes to remove wrong ones
specific = read.csv("specific_genes_detailed.csv", sep ="\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
specific_details = read.table("specific_gene_types.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F, quote = "")
#specific_details$Type[which(specific_details$Type == "wrong")] = "truncated"
specific_details = specific_details[match(specific$gene, specific_details$Gene),]
specific = cbind(specific, specific_details)


specific$cluster = as.character(specific$cluster)
specific$Type = factor(specific$Type, type_order)
specific$phylogroup = factor(specific$phylogroup,phylogroup_order)
ps1 = ggplot(specific, aes(x = cluster, fill = Type)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols)
ps = ggplot(specific, aes(x = cluster)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

## plot it proportionally? -> doesn't tell me much
specific_prop = data.frame(prop.table(table(specific$cluster, specific$Type),1))
specific_prop = cbind(specific_prop, phylogroup = specific$phylogroup[match(specific_prop$Var1, specific$cluster)])
ggplot(specific_prop, aes(x = Var1, fill = Var2, y = Freq)) + geom_bar(stat = "identity", color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols)
counts = data.frame(table(specific$cluster))


missing = read.csv("missing_genes_detailed.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
missing_details = read.table("missing_gene_types.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F, quote = "")
#missing_details$Type[which(missing_details$Type == "wrong")] = "longer"
missing_details = missing_details[match(missing$gene, missing_details$Gene),]
missing = cbind(missing, missing_details)

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
pm1 = ggplot(missing, aes(x = cluster, fill = Type, color = Type)) + geom_bar() + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols,
                    labels = c("truncated","longer","same functional group","truly missing in cluster","")) + 
  scale_color_manual(values = c(rep("black", 6),NA), guide = F) 
pm = ggplot(missing, aes(x = cluster)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



grid.arrange(ps, pm, ncol = 1)+ theme(legend.position = "None")

legend = as_ggplot(get_legend(ps1))
pm1 = pm1 + ggtitle("Missing genes") + theme(legend.position = "None")
ps1 = ps1 + ggtitle("Specific genes")+ theme(legend.position = "None")

grid.arrange(ps1, pm1, ncol = 1)



## should have the linear regression between missing and specific
grid.arrange(pm, ps, NA, legend, layout_matrix = rbind(c(1,1,3),
                                                       c(1,1,3),
                                                     c(2,2,4)))

## look at connection between number of fused genes/truncated genes and cluster size
## more genomes = more fused genes, which means it's a probability thing rather
## than an actual thing, alot of the fused genes are called in correctly
size_num_genes = data.frame(table(specific$Type, specific$cluster), stringsAsFactors = F)
colnames(size_num_genes) = c("type","cluster","count")
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",",
                           header = T, stringsAsFactors = F)
size_num_genes = cbind(size_num_genes, size = cluster_sizes$Size[match(size_num_genes$cluster, cluster_sizes$Cluster)])
for (kind in c("longer","truncated")) {
  curr = size_num_genes[which(size_num_genes$type == kind),]
  print(ggplot(curr, aes(x = size, y = count)) + geom_point() + ggtitle(kind) + 
    geom_smooth(method='lm'))
}

## check if it's true that some longer genes are not representative
gene_lengths = read.table("gene_lengths.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
longer_genes = specific$gene[which(specific$Type == "longer")]
longer_genes = gene_lengths[which(gene_lengths$Group %in% longer_genes),]
longer_genes = cbind(longer_genes, is_wrong = rep(T, dim(longer_genes)[1]))
for (gene in unique(longer_genes$Group)){
  curr = longer_genes[which(longer_genes$Group == gene),]
  mean_length_in_cluster = curr$Mean[which(curr$Group == curr$Gene)]
  sd_length_in_cluster = curr$Sd[which(curr$Group == curr$Gene)]
  other_lengths = curr$Mean[which(curr$Group != curr$Gene)]
  if (length(which(other_lengths < (mean_length_in_cluster - sd_length_in_cluster))) > 0) {
    longer_genes$is_wrong[which(longer_genes$Group == gene)] = F
  }
}
## correct the specific and missing data frames and replot them with the new values
wrong_genes = unique(longer_genes$Group[which(longer_genes$is_wrong)])
new_wrong = as.character(unlist(specific$Type))
new_wrong[which(specific$gene %in% wrong_genes)] = "new wrong"
specific = cbind(specific, new_wrong = new_wrong)
specific$new_wrong = factor(specific$new_wrong, c("new wrong",type_order))
ggplot(specific, aes(x = cluster, fill = new_wrong)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("#d3d3d3",cols), guide = F)
size_num_genes = data.frame(table(specific$new_wrong, specific$cluster), stringsAsFactors = F)
colnames(size_num_genes) = c("type","cluster","count")
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",",
                           header = T, stringsAsFactors = F)
size_num_genes = cbind(size_num_genes, size = cluster_sizes$Size[match(size_num_genes$cluster, cluster_sizes$Cluster)])
for (kind in c("longer","truncated")) {
  curr = size_num_genes[which(size_num_genes$type == kind),]
  print(ggplot(curr, aes(x = size, y = count)) + geom_point() + ggtitle(kind))
}



## tell a story on cluster 14, a non O157 EHEC that has genes from all categories
forteen = specific[which(specific$cluster == 14),]
trunc = forteen[which(forteen$new_wrong == "truncated"),]
longer = forteen[which(forteen$new_wrong == "longer"),]
same = forteen[which(forteen$new_wrong == "same functional group"),]
true = forteen[which(forteen$new_wrong %in% c("true", "unknown function")),]
wrong = forteen[which(forteen$new_wrong %in% c("wrong")),]


true_specific = specific[which(specific$new_wrong %in% c("same functional group", "truncated")),]
ggplot(true_specific, aes(x = cluster)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
