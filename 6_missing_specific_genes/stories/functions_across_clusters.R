library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(reshape2)
library(ape)


setwd("/Users/gh11/poppunk_pangenome/6_missing_specific_genes/stories/")

graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "", stringsAsFactors = F)
graphics = graphics[-which(graphics$Cluster %in% c(21, 43, 49)),]
cluster_order = graphics$Cluster

name = "outer_membrane_usher_protein"
cluster_specific = read.table(paste(name, "_cluster_specific.csv", sep = ""), sep = ",",
                              comment.char = "", stringsAsFactors = F, header = T)
cluster_specific = cbind(cluster_specific, phylo = graphics$Phylogroup[match(cluster_specific$Cluster, graphics$Cluster)])
cluster_specific$Cluster = factor(cluster_specific$Cluster, cluster_order)
ggplot(cluster_specific, aes(x = Cluster, y = Count, fill = Type)) + geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = rev(brewer.pal(n = 4, "Blues")[-1]), guide = F) + theme_classic(base_size = 14)  + ggtitle("Cluster specific") +
  scale_y_continuous(breaks = 1:max(cluster_specific$Count)) + coord_flip()

  facet_grid(Type~phylo, scales = "free",  space = "free",switch = "x")


multicluster =  read.table(paste(name, "_multicluster.csv", sep = ""), sep = ",",
                           comment.char = "", stringsAsFactors = F, header = T)
multicluster = cbind(multicluster, phylo = graphics$Phylogroup[match(multicluster$Cluster, graphics$Cluster)])
multicluster$Cluster = factor(multicluster$Cluster, cluster_order)
multicluster_matrix = data.frame(dcast(data = multicluster,formula = Cluster~Gene,fun.aggregate = sum,value.var = "Freq"), stringsAsFactors = F)
rownames(multicluster_matrix) = multicluster_matrix[,1]
multicluster_matrix = multicluster_matrix[,-1]
res = hclust(dist(t(multicluster_matrix)))
order1 = res$labels[res$order]
order = sapply(FUN = gsub, X= order1, pattern = ".", replacement = "*", fixed = T)

multicluster$Gene = factor(multicluster$Gene, order)
B = ggplot(multicluster, aes(x = Cluster, y = Gene, fill = Freq)) + geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "black", guide = F) + 
  theme_classic(base_size = 14) + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x") + ggtitle("Multicluster") 



## plot on a tree
tree = read.tree("../../7_AMR_vir_plasmid/smaller_tree/raxml_tree_mod.nwk")
plot(tree)
tree = root(tree,outgroup = "NC_011740")
tree = drop.tip(tree, tip =  "NC_011740")
for (i in c(21,43,49)){
  tree = drop.tip(tree, tip =  as.character(i))
}

res = multicluster_matrix
res = res[match(tree$tip.label, rownames(res)),]

ggtree(tree) + geom_tiplab(align = T)
## change the order of the columns
#res = res[,order(colSums(res), decreasing = T)]
#res = res[,-which(sapply(FUN = max, X = res)<0.1)] ## to remove very rare genes
p = ggtree(tree)  +
  theme(legend.position="right") + geom_tiplab(align = T)
p = gheatmap(p = p, data = res, offset = 0.1, width=5, color = "black", 
             high = "#023858", low = "#fff7fb", colnames_angle = 90, colnames = T, colnames_position = "top",
             font.size = 3, colnames_level = order1, colnames_offset_y = 5)
p

## correlate freq of particular immunity gene with pathotype
pathotypes = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_sts.csv",
                        sep = ",", header = T, stringsAsFactors = F)
multicluster = cbind(multicluster, path = pathotypes$Pathotype[match(multicluster$Cluster, pathotypes$cluster)])
multicluster = cbind(multicluster, mdr = pathotypes$MDR[match(multicluster$Cluster, pathotypes$cluster)])

curr = multicluster[multicluster$Gene == "group_2383*",]
ggplot(curr, aes(x = path, y = Freq)) + geom_boxplot()
ggplot(curr, aes(x = mdr, y = Freq)) + geom_boxplot()


