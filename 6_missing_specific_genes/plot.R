library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)


setwd("/Users/gh11/poppunk_pangenome/4_pairwise_roary/missing_specific_genes/")


#### MAIN ####


specific = read.csv("specific_genes_detailed.csv", sep ="\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
missing = read.csv("depleted_genes_detailed.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
specific_is_wrong_split = rep("No", dim(specific)[1])
missing_is_wrong_split = rep("No", dim(missing)[1])

cats = read.table("Cats.csv",sep = ",", header = F, stringsAsFactors = F, comment.char = "", quote = "")
specific = cbind(specific, cat2 = cats$V2[match(specific$classification, cats$V1)])

for (cluster in unique(c(missing$cluster,specific$cluster))) {
  example_specific = specific[which(specific$cluster == cluster),]
  example_missing = missing[which(missing$cluster == cluster),]
  for (gene in example_specific$gene){
    index = which(grepl(x = example_missing$gene, pattern = gene, fixed = T))
    if (length(index) == 1) {
      specific_is_wrong_split[which(specific$gene == gene)] = "Yes"
      missing_is_wrong_split[which(missing$gene %like% gene)] = "Yes"
    } 
    
  }
  for (gene in example_missing$gene){
    index = which(grepl(x = example_specific$gene, pattern = gene, fixed = T))
    if (length(index) == 1) {
      missing_is_wrong_split[which(missing$gene == gene)] = "Yes"
      specific_is_wrong_split[which(specific$gene %like% gene)] = "Yes"
    }
  }
}
specific = cbind(specific, is_wrong = specific_is_wrong_split)
missing = cbind(missing, is_wrong = missing_is_wrong_split)

phylogroup_order =  c("B2", "F", "D","E","A","C","B1","U")

specific$cluster = as.character(specific$cluster)
specific$phylogroup = factor(specific$phylogroup,phylogroup_order)
specific = specific[-which(specific$is_wrong == "Yes"),]
ggplot(specific, aes(x = cluster)) + geom_bar() + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


### plotting heatmaps of genes along the tree
tree = read.tree("/Users/gh11/poppunk_pangenome//6_AMR_vir_plasmid/smaller_tree/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "NC_011740")
tree = drop.tip(tree, tip =  "NC_011740")
is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]
o = tree$tip.label[ordered_tips]
o = paste("X",o, sep = "")

## Check how many genes are in the functional groups I defined and which members they have
functional_groups = read.table("specific_functional_clusters_detailed.csv", header = T, comment.char = "",
                               stringsAsFactors = F, sep = ",")

## total number of enriched genes that are also in a functional group
o = o[-which(!o %in% colnames(functional_groups))]
sum(rowSums(functional_groups[,5:50])) / dim(specific)[1]
functional_groups.m = melt(functional_groups, id.vars = c("CC","Common_terms", "Title","Cat"))
functional_groups.m$variable = factor(functional_groups.m$variable, o)
#functional_groups.m$Common_terms = factor(functional_groups.m$Common_terms, functional_groups.m$Common_terms)

ggplot(functional_groups.m, aes(x = Common_terms, y = variable, fill = value)) + geom_tile(color = "black") + 
  scale_fill_gradient(low = "white", high = "black",limits = c(0,5)) +
  facet_grid(~Cat, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(labels = rep("",dim(functional_groups.m)[1])) 


missing$cluster = as.character(missing$cluster)
missing$phylogroup = factor(missing$phylogroup, phylogroup_order)
missing = missing[-which(missing$is_wrong == "Yes"),]
ggplot(missing, aes(x = cluster, fill = is_wrong)) + geom_bar() + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0))+scale_fill_brewer(palette = "Set2") + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



gene_clean = data.frame(table(substring(text = specific$gene, first = 1, last = 4)))


## testing this with clusters 1-9 to see what I can understand from my work this morning:
corrected = read.table(file = "specific_genes_detailed_corrected.csv", sep = ',', header = T, stringsAsFactors = F, quote = "")
corrected_1 = corrected[which(corrected$cluster == 6),]
table(corrected_1$classification)

