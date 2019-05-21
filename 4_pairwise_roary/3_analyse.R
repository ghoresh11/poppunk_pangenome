library(ggplot2)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/6_pairwise_roary/")

## PLOTTING VARIABLES ###
cogs = read.table("../3_eggnog/cog_descs.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")


variable_order = c("rare","inter","core")

### Read in the outputs of the python script ###



presence_absence = read.table("030519_complete_presence_absence.csv", header = F, 
                              stringsAsFactors = F, comment.char = "", quote = "", sep = ",")
## total number of gene clusters: 
## basic plots of gene frequencies




freqs = read.table("030519_freqs.csv", header = T, row.names = NULL,
                   stringsAsFactors = F, comment.char = "", quote = "", sep =",")

#melted_freqs = read.table("melted_gene_freqs.csv", header = T, stringsAsFactors = F, comment.char = "", quote = "", sep = ",")


## TO draw a PCA plot of the genes to see if they can be separated according to specific rules
for_pca = data.frame(t(freqs[,-c(1,2)]))
write.table(x = for_pca, file = "/Users/gh11/Downloads/tsne_python/freqs_for_tsne.txt",
            col.names = F, row.names = F, quote = F)
# remove = c()
# for (i in 1:dim(for_pca)[2]) {
#   if (length(unique(for_pca[,i])) == 1) {
#     remove = c(remove, i)
#   }
# }
# for_pca = for_pca[,-remove]

cluster_properties = read.table("cluster_props.csv", sep = ",", header = T,
                                stringsAsFactors = F)
cluster_properties = cluster_properties[-which(cluster_properties$Cluster == 21),]

### PCA plot of the clusters -> what are the relationships between the clusters based on the frequencies of all genes
freqs.pca = prcomp(freqs[,-c(1,2)] , center = T, scale. = T)
#freqs.tsne = tsne(freqs[,-c(1,2)])
freqs.pca = data.frame(freqs.pca$rotation)
freqs.pca = cbind(freqs.pca, cluster_properties)
cluster_colors = c(
  brewer.pal(n = 8, "Blues")[-1],
  brewer.pal(n = 8, "Reds")[-1],
  brewer.pal(n = 8, "Greens")[-1],
  brewer.pal(n = 8, "Greys")[-1],
  brewer.pal(n = 8, "Purples")[-1],
  brewer.pal(n = 8, "PuOr")[1:4]
)
freqs.pca$Cluster = factor(freqs.pca$Cluster, 1:39)
ggplot(freqs.pca, aes(x = PC2, y = PC1, color = Cluster)) + geom_point(alpha = 0.8, size = 5) +
  scale_color_manual(values =cluster_colors, guide = F) + theme_bw(base_size = 16) +
  geom_text(aes(label=Cluster),hjust=-0.2, vjust=-0.2)


## View tsne results
tsne_results = read.table("/Users")


## retain only genes which have a mean of less than 95%
freqs_for_acc_pca = freqs[,-c(1,2)]
remove = which(unlist(as.numeric(lapply(data.frame(t(freqs_for_acc_pca)), FUN = mean))) < 0.15)
freqs_for_acc_pca = freqs_for_acc_pca[-remove,]
freqs_for_acc_pca = prcomp(freqs_for_acc_pca, center = T, scale. = T)
summary(freqs_for_acc_pca)
freqs_for_acc_pca = data.frame(freqs_for_acc_pca$rotation)
freqs_for_acc_pca = cbind(freqs_for_acc_pca, cluster_properties)
freqs_for_acc_pca$Cluster = factor(freqs_for_acc_pca$Cluster, 1:39)
ggplot(freqs_for_acc_pca, aes(x = PC1, y = PC2, color = Cluster)) + geom_point(alpha = 0.8, size = 5) +
  scale_color_manual(values =cluster_colors, guide = F) + theme_bw(base_size = 16) +
  geom_text(aes(label=Cluster),hjust=-0.2, vjust=-0.2)





between_cluster_dists = read.table("../2_dists_roary_analysis/between_cluster_dist.csv", sep = ",", 
                                   header = T, stringsAsFactors = F)
between_cluster_dists = cbind(between_cluster_dists,
                              PC1 = rep(0, dim(between_cluster_dists)[1]),
                              PC2 = rep(0, dim(between_cluster_dists)[1]))
between_cluster_dists = between_cluster_dists[-which(between_cluster_dists$cluster1 == 21 | between_cluster_dists$cluster2 == 21),]

## compare PCA to poppunk distances
for (i in 1:38) {
  for (j in (i+1):39) {
    if (i == 21 || j == 21) {next}
    index = which(between_cluster_dists$cluster1 == i &
                    between_cluster_dists$cluster2 == j)
    between_cluster_dists$PC1[index] = abs(freqs.pca$PC1[i] - freqs.pca$PC1[j])
    between_cluster_dists$PC2[index] = abs(freqs.pca$PC2[i] - freqs.pca$PC2[j])
    
  }
}

ggplot(between_cluster_dists, aes(x = core_median, y = PC2)) + geom_point(size = 2, alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Core poppunk distance") + ylab("PC2 pairwise distance")

ggplot(between_cluster_dists, aes(x = accessory_median, y = PC2)) + geom_point(size = 2, alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Accessory poppunk distance") + ylab("PC2 pairwise distance")

## Trying to answer these questions:
#1. How many rare/core genes are shared?
#2. How many core genes are shared?
gene_classes = data.frame(gene = freqs$Gene,
                          cog = freqs$COG,
                          core = rep(0, dim(freqs)[1]),
                          inter = rep(0, dim(freqs)[1]),
                          rare = rep(0, dim(freqs)[1]),
                          absent = rep(0, dim(freqs)[1]), stringsAsFactors = F)
for (i in 1:dim(freqs)[1]){
  curr_vec = freqs[i,-c(1,2)]
  gene_classes$core[i] = length(which(curr_vec >= 0.95))
  gene_classes$inter[i] = length(which(curr_vec < 0.95 & curr_vec >= 0.15))
  gene_classes$rare[i] =  length(which(curr_vec < 0.15 & curr_vec > 0))
}


write.table(x = gene_classes, file  = "030519_gene_classes.csv", sep = ",", col.names = T, row.names = F, quote = F)
gene_classes = read.table("030519_gene_classes.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "")

### check if the lack of a gene is uniformly distributed or if it's always the same cluster
## Good to have these plots, but not sure how to interpret them...
for (j in 1:37){
  thirty_eight = gene_classes$gene[which(gene_classes$core == j)]
  thirty_eight = freqs[which(freqs$Gene %in% thirty_eight),]
  vec = c()
  for (i in 1:dim(thirty_eight)[1]) {
    vec = c(vec,which(thirty_eight[i,-c(1,2)] > 0.15))
  }
  df = data.frame(freq = table(vec))
  df$freq.vec = factor(df$freq.vec, 1:39)
  p = ggplot(df, aes(x = freq.vec, y = freq.Freq)) + geom_bar(stat = "identity") + 
    theme_classic(base_size = 16) + 
    xlab("Cluster") + ylab("Number of genes not in core") + ggtitle(paste(j, "/39",sep= ""))
  print(p)
}


## This plots how many genes are core in 39/39, 38/39 etc -> but doesn't show what their other categories are
for (variable in variable_order) {
  index = which(colnames(gene_classes) == variable)
  test = gene_classes[which(gene_classes[,index]>0),] # gene was core at least once
  
  test_table = data.frame(table(test$cog, test[,index]))
  for (cog in cogs$COG){
    if (!cog %in% test$Var1) {
      test_table = rbind(test_table, data.frame(Var1 = cog, Var2 = 1, Freq = 0))
    }
  }
  test_table$Var1 = factor(test_table$Var1, cogs$COG)
  p1 = ggplot(test_table, aes(x = Var2, y = Freq, fill = Var1)) + geom_bar(stat = "identity", color = "black", size = 0.1) +
    ggtitle(variable) + scale_fill_manual(values = cogs$Col, guide = F) + theme_classic(base_size = 16) +
    ylab("Genes") + xlab(paste("Number of groups in which gene is present as", variable))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  check = 3:5
  check = check[-which(check == index)]
  test$absent[which(test[,check[1]] > 0 & test[,check[2]] > 0 )] = "mix"
  test$absent[which(test[,check[1]] > 0 & test[,check[2]] == 0 )] = colnames(test)[check[1]]
  test$absent[which(test[,check[1]] == 0 & test[,check[2]] > 0 )] = colnames(test)[check[2]]
  
  t = data.frame(table(test[,c(index,6)]))
  colnames(t) = c("Cluster", "variable", "value")
  t = rbind(t, data.frame(
    Cluster = rep(1,3),
    variable = variable_order,
    value = rep(0,3), stringsAsFactors = F
  ))
  
  t$Cluster = factor(t$Cluster, 1:39)
  t$variable = factor(t$variable, c(0, "mix", variable_order))
  
  p2 = ggplot(t, aes(x = Cluster, fill = variable, y = value)) + geom_bar(stat = "identity", size = 0.1, color = "black") +
    theme_classic(base_size = 16) + scale_fill_manual(values = c("white","#d3d3d3", brewer.pal(n=4, "Blues")[-1]),
                                                      labels = c("absent","either","rare","intermediate","core"), 
                                                      guide = F) +
    xlab("Genes") + xlab("Number of groups in which gene is present as something else") + ylab("Genes")
  grid.arrange(p1, p2, nrow = 2, ncol = 1) 
}



## Analysing entirely rare variants
rare = gene_classes[which(gene_classes$rare>0 & gene_classes$core == 0 & gene_classes$inter == 0),]
length(which(rare$rare == 1))
length(which(rare$cog %in% c("S","?")))
rare_table = data.frame(table(rare$cog, rare$rare))
rare_table$Var1 = factor(rare_table$Var1, cogs$COG)
ggplot(rare_table, aes(x = Var2, y = Freq, fill = Var1)) +geom_bar(stat = "identity", color = "black", size = 0.1) +
  scale_fill_manual(values = cogs$Col) + theme_classic(base_size = 16) +
  ylab("Genes") + xlab("Number of groups in which gene is present") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position="bottom")





## so 1,344 are core only to one group, but rare or inter in a different group
length(which(gene_classes$core == 38))
## it would be interesting to randomly remove one group from the population and see how much that changes the number of
## core genes, it should be robust to that which to me means this is the "TRUE" core of E. coli

## how many core genes are rare/intermediate in other clusters?
t = melt(core[,3:5], id.vars = "core")
t$core = factor(t$core, 1:38)
t$variable = factor(t$variable, c("rare","inter"))
ggplot(t, aes(x = core, fill = variable, y = value)) + geom_bar(stat = "identity") +
  theme_classic(base_size = 16) + scale_fill_manual(values = brewer.pal(n = 4, "Greys")[2:3]) + 
  ggtitle("core")







### find genes that are enriched or depleted in any of the groups
mean = unlist(lapply(for_pca, median))
hist(mean, breaks = 50)
variance =  unlist(lapply(for_pca, var))
plot(density(variance))
freqs = cbind(freqs, mean)
out_dir = "/Users/gh11/Submissions/ecoli_mobilome/figures/complete_presence_absence/"
count_genes = data.frame(Cluster =  rep("?",39*54),
                         COG =  rep(unique(freqs$COG), 39),
                         enriched = rep(0,39*54 ),
                         depleted = rep(0,39*54 ), stringsAsFactors = F)


for (i in 1:39){
  count_genes$Cluster[(1 + 54*(i-1)):(54*i)] = i
  enriched = which(freqs[,i+2] > freqs$mean + 0.2)
  depleted = which(freqs[,i+2] < freqs$mean - 0.2)
  for (e in enriched){
    cog = freqs$COG[e]
    count_genes$enriched[which(count_genes$Cluster == i & count_genes$COG == cog)] = 
      count_genes$enriched[which(count_genes$Cluster == i & count_genes$COG == cog)] + 1
  }
  for (d in depleted){
    cog = freqs$COG[d]
    count_genes$depleted[which(count_genes$Cluster == i & count_genes$COG == cog)] = 
      count_genes$depleted[which(count_genes$Cluster == i & count_genes$COG == cog)] + 1
  }
  
  # col = rep("no", dim(freqs)[1])
  # col[which(freqs[,i+2] > freqs$mean + 0.2 | freqs[,i+2] < freqs$mean - 0.2)] = "yes"
  # curr = data.frame(COG = freqs$COG, freq = freqs[,i+2], col = col)
  # p = ggplot(curr, aes(x = mean, y = freq, fill = col)) + 
  #   geom_point(size = 3, alpha = 0.7, pch = 21, color = "black")  +
  #   scale_fill_manual(values = c("#d3d3d3", "purple")) + theme_bw(base_size = 16) +
  #   xlab("Mean frequency") + ylab("Frequency in Cluster") + ggtitle(i)
  # fileout = paste(out_dir, i, ".pdf", sep = "")
  # ggsave(p, filename = fileout, height = 5, width = 6)
  
}
count_genes$COG = factor(count_genes$COG, cogs$COG)
count_genes$Cluster = factor(count_genes$Cluster, 1:39)
ggplot(count_genes, aes(x = Cluster, y = enriched, fill = COG)) + 
  geom_bar(stat = "identity") + theme_classic(base_size = 16) +
  scale_fill_manual(values = cogs$Col, guide = F) + ggtitle("Enriched genes")

ggplot(count_genes, aes(x = Cluster, y = depleted, fill = COG)) + 
  geom_bar(stat = "identity") + theme_classic(base_size = 16) +
  scale_fill_manual(values = cogs$Col, guide = F) + ggtitle("Depleted genes")
