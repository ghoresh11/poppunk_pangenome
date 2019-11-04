library(ggplot2)
library(RColorBrewer)
library(vegan)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(reshape2)
library(ape)

setwd("/Users/gh11/poppunk_pangenome/6_missing_specific_genes//")

gene_classes = read.table("../4_pairwise_roary/231019_corrected/gene_classes.csv",
                          sep = ",", comment.char = "", stringsAsFactors = F, header = T,
                          quote = "")
gene_assignments = read.table("../5_classify_genes/gene_classification.csv", sep = ",", header = F,
                              comment.char = "", quote = "", stringsAsFactors = F)
freqs = read.table("../4_pairwise_roary/231019_corrected//freqs.csv", header = T,row.names = 1,
                   stringsAsFactors = F, comment.char = "", quote = "", sep =",")
freqs = freqs[,-which(colnames(freqs) == "X50")]
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "", stringsAsFactors = F)

## genes which are missing only in one cluster
missing_in_one = gene_assignments$V1[which(gene_assignments$V2 == "Missing in one")]
## for the ubiq genes, which clusters are they missing from
missing_in_one = freqs[which(rownames(freqs) %in% missing_in_one),]


summary_depleted = data.frame(cluster =  colnames(missing_in_one),
                              num_genes = rep(0, dim(missing_in_one)[2]),
                              missing_genes = rep("", dim(missing_in_one)[2]), stringsAsFactors = F)
for (i in 1:dim(missing_in_one)[1]){
  not_core = which(missing_in_one[i,]<0.15) ## which clusters are these genes in very low frequency
  for (j in not_core){
    summary_depleted$missing_genes[j] = paste(summary_depleted$missing_genes[j], rownames(missing_in_one)[i], sep = "/")
    summary_depleted$num_genes[j] = summary_depleted$num_genes[j] + 1
  }
}

summary_depleted$cluster =  sapply(summary_depleted$cluster, FUN = gsub, pattern = "X", replacement = "", fixed = T)
summary_depleted = cbind(summary_depleted, phylogroup = graphics$Phylogroup[match(summary_depleted$cluster, graphics$Cluster)])
summary_depleted$cluster = factor(summary_depleted$cluster, rev(graphics$Cluster))
summary_depleted$phylogroup = factor(summary_depleted$phylogroup, c("B2", "F", "D","E","A","C","B1","U"))
write.table(summary_depleted, file = "../6_missing_specific_genes/missing_genes.csv", sep = ",", col.names = T, row.names = F, quote = F)

p1 = ggplot(summary_depleted, aes(x = cluster, y = num_genes)) + geom_bar(stat = "identity") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0))+
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


## genes which are missing only in one cluster
specific = gene_assignments$V1[which(gene_assignments$V2 == "Core and specific")]
specific = freqs[which(rownames(freqs) %in% specific),]
summary_specific = data.frame(cluster =  colnames(specific),
                              num_genes = rep(0, dim(specific)[2]),
                              specific = rep("", dim(specific)[2]), stringsAsFactors = F)
for (i in 1:dim(specific)[1]){
  core = which(specific[i,]>=0.95) ## which clusters are these genes in very low frequency
  for (j in core){
    summary_specific$specific[j] = paste(summary_specific$specific[j], rownames(specific)[i], sep = "/")
    summary_specific$num_genes[j] = summary_specific$num_genes[j] + 1
  }
}
summary_specific$cluster =  sapply(summary_specific$cluster, FUN = gsub, pattern = "X", replacement = "", fixed = T)
summary_specific = cbind(summary_specific, phylogroup = graphics$Phylogroup[match(summary_specific$cluster, graphics$Cluster)])
summary_specific$cluster = factor(summary_specific$cluster, rev(graphics$Cluster))
summary_specific$phylogroup = factor(summary_specific$phylogroup, c("B2", "F", "D","E","A","C","B1","U"))
write.table(summary_specific, file = "../6_missing_specific_genes//specific_genes.csv", sep = ",", col.names = T, row.names = F, quote = F)

p2 = ggplot(summary_specific, aes(x = cluster, y = num_genes)) + geom_bar(stat = "identity") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0))+
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p1
p2
## relationship between number of depleted genes vs number of enriched genes

df = data.frame(cluster = summary_depleted$cluster,
                missing = summary_depleted$num_genes, specific = summary_specific$num_genes, stringsAsFactors = F)
lin = ggplot(df, aes(x = missing, y = specific)) + geom_point(alpha = 0.6, size = 1.4) +
  theme_classic(base_size = 14) +geom_text(aes(label=cluster),hjust=0, vjust=0) +
  xlab("Missing genes") + ylab("Specific genes")+ 
  geom_smooth(method='lm', col = 'black') + scale_x_continuous(expand = c(0.01,0,0,4))
lin
linear_model = lm(formula = specific~missing, data = df)
summary(linear_model)

max(df$missing)
max(df$specific)
median(df$missing)
median(df$specific)
min(df$missing)
min(df$specific)

ubiq = gene_classes[which(gene_classes$total == 47,),]
specific = gene_classes[which(gene_classes$total == 1,),]

## see what the connection is between number of secondary variants and 
## missing / specific genes
secondary = gene_assignments$V1[which(grepl(gene_assignments$V2,pattern = "Secondary"))]
secondary = freqs[which(rownames(freqs) %in% secondary),]
summary_secondary = data.frame(cluster =  colnames(secondary),
                              num_genes = rep(0, dim(secondary)[2]), stringsAsFactors = F)
for (i in 1:dim(secondary)[1]){
  core = which(secondary[i,]>=0.5) ## which clusters have this gene in at least 50% of genomes 
  for (j in core){
    summary_secondary$num_genes[j] = summary_secondary$num_genes[j] + 1
  }
}
summary_secondary$cluster =  sapply(summary_secondary$cluster, FUN = gsub, pattern = "X", replacement = "", fixed = T)
summary_secondary = summary_secondary[match(summary_specific$cluster, summary_secondary$cluster),]
summary_depleted = summary_depleted[match(summary_specific$cluster, summary_depleted$cluster),]


plot(summary_depleted$num_genes, summary_secondary$num_genes)
plot(summary_specific$num_genes, summary_secondary$num_genes)

## generally speaking there's still a connection between all of these which to me suggests that a lot
## of the lineage specific genes are short variants of other genes that are no longer in the population









