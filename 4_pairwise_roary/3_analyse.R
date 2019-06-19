library(ggplot2)
library(RColorBrewer)
library(vegan)
library(reshape)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/4_pairwise_roary/")

variable_order = c("rare","inter", "core")


## graphics
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "", stringsAsFactors = F)

# 
# presence_absence = read.table("170619/complete_presence_absence.csv", header = F, 
#                               stringsAsFactors = F, comment.char = "", quote = "", sep = ",")


freqs = read.table("180619/freqs.csv", header = T,row.names = 1,
                   stringsAsFactors = F, comment.char = "", quote = "", sep =",")
o = as.numeric(unlist(lapply(colnames(freqs), gsub, pattern = "X", replacement = "")))
graphics = graphics[match(o, graphics$Cluster),]

#melted_freqs = read.table("melted_gene_freqs.csv", header = T, stringsAsFactors = F, comment.char = "", quote = "", sep = ",")


## TO draw a PCA plot of the genes to see if they can be separated according to specific rules
for_pca = t(freqs)
# write.table(x = for_pca, file = "/Users/gh11/Downloads/tsne_python/freqs_for_tsne.txt",
#             col.names = F, row.names = F, quote = F)
remove = c()
for (i in 1:dim(for_pca)[2]) {
  if (length(unique(for_pca[,i])) == 1) {
    remove = c(remove, i)
  }
}
for_pca = for_pca[,-remove]

### PCA plot of the clusters -> what are the relationships between the clusters based on the frequencies of all genes
freqs.pca = prcomp(freqs , center = T, scale. = T)
summary(freqs.pca)
freqs.pca = data.frame(freqs.pca$rotation)
freqs.pca = cbind(freqs.pca, Cluster = as.character(o))
freqs.pca$Cluster = factor(freqs.pca$Cluster, o)


## PCA plot
ggplot(freqs.pca, aes(x = PC2, y = PC1, color = Cluster, shape = Cluster)) + geom_point(size = 3.5, stroke = 1, alpha = 0.7) +
  scale_color_manual(values = graphics$Color, guide = F) + theme_bw(base_size = 16) +
  geom_text(aes(label=Cluster),hjust=-0.3, vjust=-0.3) +
  scale_shape_manual(values =  graphics$Shape, guide = F)


## retain only genes which have a mean of less than 95%
## check if there is a different division based on rare vs core genes
## change here for genes in a specific "type" -> rare/core/inter
freqs_for_acc_pca = for_pca
means = apply(freqs_for_acc_pca, 2, FUN = mean)
means = as.numeric(unlist(means))
names(means) = colnames(freqs_for_acc_pca)
keep = which(means < 0.95 & means > 0.15)
#keep = which(means >= 0.95)
#keep = which(means  <= 0.15)
freqs_for_acc_pca = t(freqs_for_acc_pca[,keep])
freqs_for_acc_pca = prcomp(freqs_for_acc_pca, center = T, scale. = T)
summary(freqs_for_acc_pca)
freqs_for_acc_pca = data.frame(freqs_for_acc_pca$rotation)
freqs_for_acc_pca = cbind(freqs_for_acc_pca,Cluster =  o)
freqs_for_acc_pca$Cluster = factor(freqs_for_acc_pca$Cluster, o)
ggplot(freqs_for_acc_pca, aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) + geom_point(size = 3.5, stroke = 1, alpha = 0.7)  +
  scale_color_manual(values =graphics$Color, guide = F) + theme_bw(base_size = 16) +
  geom_text(aes(label=Cluster),hjust=-0.2, vjust=-0.2) +
  scale_shape_manual(values =graphics$Shape, guide = F)


## Trying to answer these questions:
#1. How many rare/core genes are shared?
#2. How many core genes are shared?
gene_classes = data.frame(gene = rownames(freqs),
                          core = rep(0, dim(freqs)[1]),
                          inter = rep(0, dim(freqs)[1]),
                          rare = rep(0, dim(freqs)[1]), stringsAsFactors = F)
for (i in 1:dim(freqs)[1]){
  curr_vec = freqs[i,]
  gene_classes$core[i] = length(which(curr_vec >= 0.95))
  gene_classes$inter[i] = length(which(curr_vec < 0.95 & curr_vec >= 0.15))
  gene_classes$rare[i] =  length(which(curr_vec < 0.15 & curr_vec > 0))
}


write.table(x = gene_classes, file  = "180619/gene_classes.csv", sep = ",", col.names = T, row.names = F, quote = F)
gene_classes = read.table("180619/gene_classes.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "", row.names = 1)

gene_classes = cbind(gene_classes, total_presence = rowSums(gene_classes))
total_presence = data.frame(table(gene_classes$total_presence))



## examples of genes
curr = data.frame(name = o,
                  freq = unlist(freqs[which(rownames(freqs) == "group_6968"),]), stringsAsFactors = F)
curr$name = factor(curr$name, o)
ggplot(curr, aes(x = name, y = freq)) + geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Frequency") + xlab("Cluster")

## this is what I would call the core -> these genes are found across all the cluster
## but there are a few clusters that are missing some key genes

df = data.frame(variable = character(0), type = character(0),value = numeric(0),stringsAsFactors = F)
for (size in 1:47) {
  omni = gene_classes[which(gene_classes$total_presence == size),]
  omni = data.frame(table(omni$core, omni$inter, omni$rare), stringsAsFactors = F)
  colnames(omni) = c("core","inter","rare","freq")
  #omni = omni[-which(omni$freq == 0),]
  for (i in 1:3){
    omni[,i] = as.numeric(as.character(omni[,i]))
    omni[,i] = (omni[,i]/size)*omni$freq
  }
  omni_percent = colSums(omni)[-4] / length(which(gene_classes$total_presence == size))
  omni_percent = data.frame(variable = rep(paste(size,"/47",sep=""),3),type = names(omni_percent), value = omni_percent)
  df = rbind(df, omni_percent)
}
p2 = ggplot(df, aes(x = variable, y = value, fill = type)) + geom_bar(stat = "identity", color = "black") +
  theme_classic(base_size = 14) + scale_fill_manual(values = rev(brewer.pal(n=4,"Blues"))) + theme(legend.position = "bottom")+
  xlab("Number of clusters in which gene is present") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Fraction of genes")

total_presence$Var1 = unique(df$variable)
p1 = ggplot(total_presence, aes(x = total_presence$Var1, y = total_presence$Freq)) + geom_bar(stat = "identity") +
  xlab("Number of clusters in which gene is present") + ylab("Number of genes") +
  theme_bw(base_size = 16)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


p1
grid.arrange(p1,p2)
### check if the lack of a gene is uniformly distributed or if it's always the same cluster
# ## Good to have these plots, but not sure how to interpret them...
# for (j in 1:37){
#   thirty_eight = gene_classes$gene[which(gene_classes$core == j)]
#   thirty_eight = freqs[which(freqs$Gene %in% thirty_eight),]
#   vec = c()
#   for (i in 1:dim(thirty_eight)[1]) {
#     vec = c(vec,which(thirty_eight[i,-c(1,2)] > 0.15))
#   }
#   df = data.frame(freq = table(vec))
#   df$freq.vec = factor(df$freq.vec, 1:39)
#   p = ggplot(df, aes(x = freq.vec, y = freq.Freq)) + geom_bar(stat = "identity") + 
#     theme_classic(base_size = 16) + 
#     xlab("Cluster") + ylab("Number of genes not in core") + ggtitle(paste(j, "/39",sep= ""))
#   print(p)
# }


ggplot(gene_classes, aes(x = core, y = rare)) + geom_point(position = "jitter", alpha = 0.4) +
  xlab("Number of clusters where gene is core") + ylab("Number of clusters where gene is rare")


## This plots how many genes are core in 39/39, 38/39 etc -> but doesn't show what their other categories are
for (variable in variable_order) {
  index = which(colnames(gene_classes) == variable)
  test = gene_classes[which(gene_classes[,index]>0),] # gene was core at least once
  
  test_table = data.frame(table(test$cog, test[,index]))
  for (cog in cogs$COG){
    if (!cog %in% test$Var1) {
      test_table = rbind(test_table, data.frame( Var2 = 1, Freq = 0))
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
means = apply(freqs, 1, FUN = median)
freqs = cbind(freqs, means)
# hist(means, breaks = 50)
# variance =  apply(freqs, 1, FUN = var)
# plot(density(variance))

counts = data.frame(cluster = numeric(0),
                    num_depeleted = numeric(0),
                    num_enriched = numeric(0), stringsAsFactors = F)

for (i in 1:(dim(freqs)[2]-1)){
  # col = rep("same", dim(freqs)[1])
  # col[ which(abs(freqs$means - freqs[,i]) > 0.2)] = "diff"
  # df = data.frame(median = freqs$means,
  #                 cluster = freqs[,i], col = col, stringsAsFactors = F)
  # p = ggplot(df, aes(x = median, y = cluster, color = col)) + geom_point(size = 2, alpha = 0.8) +
  #   theme_bw(base_size = 12) + xlab("Mean frequency") + ylab("Frequency in cluster") +
  #   scale_color_manual(values = c("red","black"), guide = F) + ggtitle(i)
  # fileout = paste("/Users/gh11/Submissions/my_thesis/Chapter3/figures/6_complete_presence_absence/freqs_per_cluster/",i, ".pdf",
  #                 sep = "")
  # ggsave(plot = p, filename = fileout, width = 6, height = 6)
  
  counts = rbind(counts, data.frame(cluster = o[i], 
                                    num_depleted = length(which(freqs$means - freqs[,i] > 0.8)), 
                                    num_enriched = length(which(freqs[,i] - freqs$means > 0.8))))
}
counts$cluster = factor(counts$cluster,o)
p1 = ggplot(counts, aes(x = cluster, y = num_enriched)) + 
  geom_bar(stat = "identity") + theme_classic(base_size = 16) + ggtitle("Enriched genes") +
  xlab("Genes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

p2 = ggplot(counts, aes(x = cluster, y = num_depleted)) + 
  geom_bar(stat = "identity") + theme_classic(base_size = 16) + ggtitle("Depleted genes") +
  ylab("Genes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

grid.arrange(p1, p2)

length(which(freqs$X12 - freqs$means > 0.9))


