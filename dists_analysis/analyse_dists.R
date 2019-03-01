library(ggplot2)
library(RColorBrewer)
library(reshape2)

setwd("/Users/gh11/poppunk_pangenome/dists_analysis/")


cluster_sizes = read.table("cluster_sizes.csv", sep = ",",
                           header = T, stringsAsFactors = F)

## draw a point plot with the cluster sizes
cluster_sizes = cluster_sizes[order(cluster_sizes$Size, decreasing = T),]
cluster_sizes = cbind(ID = 1:dim(cluster_sizes)[1], cluster_sizes, Cumsum = cumsum(cluster_sizes$Size))

ggplot(cluster_sizes, aes(x = ID ,y = cluster_sizes$Cumsum)) + geom_point(alpha = 0.5, size = 3) +
  xlab("poppunk cluster") + ylab("Total number of genomes") + geom_vline(xintercept=39, col = "red") +
  scale_y_continuous(limits = c(0, 13500)) + theme_bw(base_size = 16)

## this goes up to 519, because from 520 the clusters are of size 1 so there's no distance
within_distances = read.table("all/within_cluster_dist.csv", sep = ",",
                           header = T, stringsAsFactors = F, comment.char = "")

## because of an if statement in "calc_cluster_dists", it actually only took clusters
## of size 2 or more and ignore clusters of size 1
between_distances = read.table("all/between_cluster_dist.csv", sep = ",",
                               header = T, stringsAsFactors = F, comment.char = "")

chosen = cluster_sizes$Cluster[which(cluster_sizes$Size >= 30)]

## add a color specifying which clusters I kept in the analysis
chosen_within = rep("no",dim(within_distances)[1])
chosen_within[which(within_distances$Cluster %in% chosen)] = "yes"
within_distances = cbind(within_distances, chosen = chosen_within)

chosen_between = rep("no",dim(between_distances)[1])
chosen_between[which(between_distances$cluster1 %in% chosen & between_distances$cluster2 %in% chosen)] = "yes"
between_distances = cbind(between_distances, chosen = chosen_between)

### see the range of distance
ggplot(within_distances, aes(x = Core, fill = chosen)) + geom_density(alpha = 0.7) +
  theme_classic(base_size = 16) + scale_fill_manual(values = brewer.pal(n = 3,"Greys")[c(1,3)]) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Within cluster mean core distance")
ggplot(within_distances, aes(x = Core_max, fill = chosen)) + geom_density(alpha = 0.5)
ggplot(within_distances, aes(x = Acc_max,  fill = chosen)) + geom_density(alpha = 0.5)
ggplot(within_distances, aes(x = Acc, fill = chosen)) + geom_density(alpha = 0.7) +
  theme_classic(base_size = 16) + scale_fill_manual(values = brewer.pal(n = 3,"Greys")[c(1,3)]) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Within cluster mean accessory distance")


ggplot(between_distances, aes(x = core, fill = chosen)) + geom_density(alpha = 0.7) +
  theme_classic(base_size = 16) + scale_fill_manual(values = brewer.pal(n = 3,"Greys")[c(1,3)]) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Between cluster mean core distance")
ggplot(between_distances, aes(x = core_max, fill = chosen)) + geom_density(alpha = 0.5)
ggplot(between_distances, aes(x = accessory, fill = chosen)) + geom_density(alpha = 0.7) +
  theme_classic(base_size = 16) + scale_fill_manual(values = brewer.pal(n = 3,"Greys")[c(1,3)]) + 
  scale_y_continuous(expand = c(0,0)) + ylab("Density") + xlab("Between cluster mean accessory distance")
ggplot(between_distances, aes(x = accessory_max, fill = chosen)) + geom_density(alpha = 0.5)

### get general values
max(within_distances$Core_max)
max(within_distances$Acc_max)
max(between_distances$core_max)
max(between_distances$accessory_max)

### what's the connection between the size of the cluster and the core
## and accessory distances? There isn't a direct connection
max(roary_outputs$core_max_dist)
max(roary_outputs$acc_max_dist)

within_distances = cbind(within_distances, size = cluster_sizes$Size[match(within_distances$Cluster, cluster_sizes$Cluster)])
ggplot(within_distances, aes(y = Core, x = size)) + geom_point(alpha = 0.4, size = 2) +
  theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Within cluster core distance")

ggplot(within_distances, aes(y = Acc, x = size)) + geom_point(alpha = 0.4, size = 2) +
  theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Within cluster accessory distance")


#### see how the number of genes in each cluster compares to the 
# 1. Size of the cluster
# 2. Within cluster acc dist
# 3. Within cluster core dist

roary_outputs = read.table("roary_summary_file.csv", sep = ",", header = T, stringsAsFactors = F)
## add info to the roary outputs CSV
sizes = cluster_sizes$Size[match(roary_outputs$cluster, cluster_sizes$Cluster)]
core_dist = within_distances$Core[match(roary_outputs$cluster, within_distances$Cluster)]
core_max_dist = within_distances$Core_max[match(roary_outputs$cluster, within_distances$Cluster)]
acc_dist = within_distances$Acc[match(roary_outputs$cluster, within_distances$Cluster)]
acc_max_dist = within_distances$Acc_max[match(roary_outputs$cluster, within_distances$Cluster)]
roary_outputs = cbind(roary_outputs, sizes, core_dist, core_max_dist, acc_dist, acc_max_dist)


roary_outputs.m = melt(id.vars  = "cluster", roary_outputs[,1:5])
roary_outputs.m = cbind(roary_outputs.m,
                        size = roary_outputs$sizes[match(roary_outputs.m$cluster, roary_outputs$cluster)])

roary_outputs.m$variable = factor(roary_outputs.m$variable, c("rare_genes","inter_genes",  "soft_core_genes", "core_genes"))
labs = c("Rare", "Intermediate", "Soft core", "Core")
roary_outputs.m$cluster = factor(roary_outputs.m$cluster, sort(as.numeric(unique(roary_outputs$cluster))))

ggplot(roary_outputs.m, aes(x = cluster, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Blues") + theme_classic(base_size = 16) +
  scale_y_continuous(expand = c(0,0)) + xlab("Cluster") +
  ylab("Genes")

roary_outputs.m$variable = factor(roary_outputs.m$variable, c("core_genes", "soft_core_genes", "inter_genes","rare_genes"))
ggplot(roary_outputs.m, aes(x = variable, y = value, fill = variable))+ 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Blues", guide = F, direction = -1) +
  theme_bw(base_size = 16) + ylab("Genes") + xlab("") + 
  scale_x_discrete(labels = rev(labs))


median(roary_outputs$core_genes)
sd(roary_outputs$core_genes)
median(roary_outputs$soft_core_genes)
sd(roary_outputs$soft_core_genes)
median(roary_outputs$inter_genes)
sd(roary_outputs$inter_genes)

roary_outputs.m$variable = factor(roary_outputs.m$variable, c("rare_genes","inter_genes",  "soft_core_genes", "core_genes"))
ggplot(roary_outputs.m, aes(y = value, x = size, fill = variable)) + geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Cluster size") + ylab("Genes") +
  scale_fill_brewer(palette = "Blues", guide = F)+ 
  geom_smooth(method='lm',formula= y~x, se = T, aes(color = variable)) +
  scale_color_brewer(palette = "Blues", guide = F)

## add distances
roary_outputs.m = cbind(roary_outputs.m,
                        acc_dist = roary_outputs$acc_dist[match(roary_outputs.m$cluster, roary_outputs$cluster)])
ggplot(roary_outputs.m, aes(y = value, x = acc_dist, fill = variable)) + geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Within cluster accessory distance") + ylab("Genes") +
  scale_fill_brewer(palette = "Blues", guide = F) + 
  geom_smooth(method='lm',formula= y~x, se = T, aes(color = variable)) +
  scale_color_brewer(palette = "Blues", guide = F)


roary_outputs.m = cbind(roary_outputs.m,
                        core_dist = roary_outputs$core_dist[match(roary_outputs.m$cluster, roary_outputs$cluster)])
ggplot(roary_outputs.m, aes(y = value, x = core_dist, fill = variable)) + geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 16) + xlab("Within cluster core distance") + ylab("Genes") +
  scale_fill_brewer(palette = "Blues", guide = F) + 
  geom_smooth(method='lm',formula= y~x, se = T, aes(color = variable)) +
  scale_color_brewer(palette = "Blues", guide = F)



### look at the gene lengths within each type of gene (rare, core, soft_core ...)
gene_lengths = read.table("Gene_lengths.csv", sep = ",", header = T,
                          stringsAsFactors = F)
gene_lengths$Type = factor(gene_lengths$Type , rev(c("rare","inter","soft_core", "core")))
labs = c()
for (type in rev(c("rare","inter",  "soft_core", "core"))){
  labs = c(labs,paste(type,"\n(",length(which(gene_lengths$Type == type)), ")", sep = ""))
}
ggplot(gene_lengths, aes(x = Type, y = Length, fill = Type)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Blues", guide = F, direction = -1) +
  theme_bw(base_size = 16) + ylab("Gene length (aa)") + xlab("") + 
  scale_x_discrete(labels = labs)

median(gene_lengths$Length[which(gene_lengths$Type == "core")])
median(gene_lengths$Length[which(gene_lengths$Type == "soft_core")])
median(gene_lengths$Length[which(gene_lengths$Type == "inter")])
median(gene_lengths$Length[which(gene_lengths$Type == "rare")])

## stratify by cluster
gene_lengths$Cluster = factor(gene_lengths$Cluster, 1:39)
core_lengths = gene_lengths[which(gene_lengths$Type == "core"),]
ggplot(core_lengths, aes(x = Cluster, y = Length)) + geom_boxplot()
soft_core_lengths = gene_lengths[which(gene_lengths$Type == "soft_core"),]
ggplot(soft_core_lengths, aes(x = Cluster, y = Length)) + geom_boxplot()
inter = gene_lengths[which(gene_lengths$Type == "inter"),]
ggplot(inter, aes(x = Cluster, y = Length)) + geom_boxplot()
rare = gene_lengths[which(gene_lengths$Type == "rare"),]
ggplot(rare, aes(x = Cluster, y = Length)) + geom_boxplot()



## zoom in 
gene_lengths = gene_lengths[-which(gene_lengths$Length > 500), ]
labs = c()
for (type in rev(c("rare","inter",  "soft_core", "core"))){
  labs = c(labs,paste(type,"\n(",length(which(gene_lengths$Type == type)), ")", sep = ""))
}
ggplot(gene_lengths, aes(x = Type, y = Length, fill = Type)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  scale_fill_brewer(palette = "Blues", guide = F, direction = -1) +
  theme_bw(base_size = 16) + ylab("Gene length (aa)") + xlab("") + 
  scale_x_discrete(labels = labs)

ggplot(gene_lengths, aes(x = Length, fill = Type)) + geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "Blues", guide = F, direction = -1)  +
  theme_classic(base_size = 16) + xlab("Gene length (aa)") + ylab("Density")  



