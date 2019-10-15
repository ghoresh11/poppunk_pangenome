library(ggplot2)
library(RColorBrewer)
library(vegan)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(reshape2)
library(ape)

setwd("/Users/gh11/poppunk_pangenome/4_pairwise_roary/")

variable_order = c("rare","inter", "core")


## graphics
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "", stringsAsFactors = F)
freqs = read.table("071019_mode_rep/freqs.csv", header = T,row.names = 1,
                   stringsAsFactors = F, comment.char = "", quote = "", sep =",")
clusters = sapply(X = colnames(freqs), FUN = gsub, pattern = "X", replacement = "")
graphics = graphics[match(clusters, graphics$Cluster),]
o = clusters

for_pca = t(freqs)
remove = c()
for (i in 1:dim(for_pca)[2]) {
  if (length(unique(for_pca[,i])) == 1) { ## no variation in gene
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
B = ggplot(freqs.pca, aes(x = PC2, y = PC3, color = Cluster, shape = Cluster)) + geom_point(size = 3.5, stroke = 1, alpha = 0.7) +
  scale_color_manual(values = graphics$Color, guide = F) + theme_bw(base_size = 12) +
  geom_text(aes(label=Cluster),hjust=-0.3, vjust=-0.3) +
  scale_shape_manual(values =  graphics$Shape, guide = F) + ggtitle("C") + scale_x_continuous(expand = c(0.03,0.03))+
  scale_y_continuous(expand = c(0.02,0.02))   + 
   annotate("text", x = -0.05, y = 0.1 , label = "bold(A/B1)", size = 5, parse = T) + 
   annotate("text", x = -0.05, y = -0.12 , label = "bold(E)", size = 5,parse = T)+ 
   annotate("text", x = 0.05, y = -0.16 , label = "bold(F)", size = 5,parse = T)+ 
   annotate("text", x = 0.17, y = 0.08 , label = "bold(B2)", size = 5,parse = T)


## Trying to answer these questions:
# 1. How many rare/core genes are shared?
# 2. How many core genes are shared?
# gene_classes = data.frame(gene = rownames(freqs),
#                           core = rep(0, dim(freqs)[1]),
#                           inter = rep(0, dim(freqs)[1]),
#                           rare = rep(0, dim(freqs)[1]), stringsAsFactors = F)
# for (i in 1:dim(freqs)[1]){
#   curr_vec = freqs[i,]
#   gene_classes$core[i] = length(which(curr_vec >= 0.95))
#   gene_classes$inter[i] = length(which(curr_vec < 0.95 & curr_vec >= 0.15))
#   gene_classes$rare[i] =  length(which(curr_vec < 0.15 & curr_vec > 0))
# }
#write.table(x = gene_classes, file  = "071019_mode_rep//gene_classes.csv", sep = ",", col.names = T, row.names = F, quote = F)

gene_classes = read.table("071019_mode_rep/gene_classes.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "", row.names = 1)
gene_classes = cbind(gene_classes, total_presence = rowSums(gene_classes))
total_presence = data.frame(table(gene_classes$total_presence))

### get the tree object
tree = read.tree("../7_AMR_vir_plasmid/smaller_tree/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "NC_011740")
tree = drop.tip(tree, tip =  "NC_011740")
for (i in c(21,43,49)){
  tree = drop.tip(tree, tip =  as.character(i))
}
p = ggtree(tree)  +
  theme(legend.position="right") + geom_tiplab(align = T)

## examples of genes
gene_name = "intA_1" ## change here to check the frequency of a gene across the popPUNK clusters
curr = data.frame(name = o,
                  freq = unlist(freqs[which(rownames(freqs) == gene_name),]), stringsAsFactors = F)
curr = cbind(curr, phylo = graphics$Phylogroup[match(curr$name, graphics$Cluster)])
curr$name = factor(curr$name, o)
ggplot(curr, aes(x = name, y = freq)) + geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Frequency") + xlab("PopPUNK cluster") + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x")  +
  ggtitle("A - intA")

# facet_plot(p, panel = 'intA', data = curr, 
#            geom = geom_barh, 
#            mapping = aes(x = freq, fill = curr$name), 
#            stat='identity' ) +theme_minimal()


gene_name = "group_4083" ## change here to check the frequency of a gene across the popPUNK clusters
curr = data.frame(name = o,
                  freq = unlist(freqs[which(rownames(freqs) == gene_name),]), stringsAsFactors = F)
curr = cbind(curr, phylo = graphics$Phylogroup[match(curr$name, graphics$Cluster)])
curr$name = factor(curr$name, o)
ggplot(curr, aes(x = name, y = freq)) + geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Frequency") + xlab("PopPUNK cluster") + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x")  +
  ggtitle("C - group_4083")

# facet_plot(p1, panel = 'wzyE', data = curr, 
#                   geom = geom_hbar, 
#                   mapping = aes(x = as.numeric(freq), ), 
#                   stat='identity' ) +theme_minimal()


## plot presence and absence of two genes -> i.e. which clusters have one gene and not another - on a tree
gene_names = c("group_5441", "group_2329*************", "group_1059**")
curr = data.frame(name = rep(o, length(gene_names)),
                  gene = rep("", length.out = length(o)*length(gene_names)),
                  freq = rep(0, length.out = length(o)*length(gene_names)), stringsAsFactors = F)
for (i in 1:length(gene_names)) {
  curr_gene = gene_names[i]
  curr$gene[(length(o)*(i-1) + 1): (length(o)*i)] = rep(curr_gene, length(o))
  curr$freq[(length(o)*(i-1) + 1): (length(o)*i)] = unlist(freqs[which(rownames(freqs) == curr_gene),])
}
curr$name = factor(curr$name, o)
curr = cbind(curr, phylo = graphics$Phylogroup[match(curr$name, graphics$Cluster)])

res = dcast(curr, formula = name~gene, value.var = "freq")
rownames(res) = res[,1]
res = res[,-1]

gheatmap(p = p, data = res, offset = 0.1, width=5, color = "black", 
             high = "#023858", low = "#fff7fb", colnames_angle = 90, colnames = T, colnames_position = "top",
             font.size = 3, colnames_offset_y = 5)


ggplot(curr,aes(x = name, y = gene, fill = freq)) + geom_tile(stat = "identity", color = "black") + 
  theme_classic(base_size = 14) + scale_fill_gradient(low = "white", high = "black") +
  scale_y_discrete(labels = c("effector","hypothetical\nprotein","immunity")) +
  xlab("PopPUNK cluster") + 
  facet_grid(~phylo, scales = "free",  space = "free",switch = "x")


## this is what I would call the core -> these genes are found across all the cluster
## but there are a few clusters that are missing some key genes
df = data.frame(variable = character(0), type = character(0),value = numeric(0),stringsAsFactors = F)
for (size in 1:max(gene_classes$core)) {
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
p2 =  ggplot(df, aes(x = variable, y = value, fill = type)) + geom_bar(stat = "identity", color = "black") +
  theme_classic(base_size = 12) + scale_fill_manual(values = rev(brewer.pal(n=4,"Blues"))) + theme(legend.position = "bottom")+
  xlab("Number of clusters in which gene is present") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Fraction of genes")
p2
total_presence = cbind(total_presence, lab = paste(total_presence$Var1,"/47",sep=""))
total_presence = total_presence[-which(total_presence$Var1 == 0),]
p1 =  ggplot(total_presence, aes(x = total_presence$Var1, y = total_presence$Freq)) + geom_bar(stat = "identity") +
  xlab("Number of clusters in which gene is present") + ylab("Number of genes") +
  theme_bw(base_size = 12)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels = total_presence$lab) + ggtitle("D") +
  annotate("text", x = 4, y = 60000, label = format(total_presence$Freq[which(total_presence$Var1 == 1)], big.mark = ",", scientific = F)) + 
  annotate("text", x = 45, y = 4000, label = format(total_presence$Freq[which(total_presence$Var1 == 47)], big.mark = ",", scientific = F))

p1

grid.arrange(intA,wzyE, p1, B, layout_matrix = rbind(c(1,1,1,2,2,2),
                                                    c(1,1,1,2,2,2),
                                                    c(4,4,3,3,3,3),
                                                    c(4,4,3,3,3,3),
                                                    c(4,4,3,3,3,3)))

#### Ubiq Figure ####

plot_for_num <- function(num, p, lay){
  ubiq = gene_classes[which(gene_classes$total_presence == num),]
  ubiq = cbind(ubiq, pattern = apply(FUN = paste, X= ubiq[,1:3],1,  collapse = "-"))
  for_barplot = data.frame(table(ubiq$pattern), stringsAsFactors = F)
  for_barplot = for_barplot[order(for_barplot$Freq, decreasing = T),]
  ## first make the heatmap
  for_hp = ubiq[match(for_barplot$Var1, ubiq$pattern),]
  for_hp = melt(for_hp[,c(5,1,2,3)], id.vars = "pattern")
  for_hp$pattern = factor(for_hp$pattern, for_barplot$Var1)
  for_hp$variable = factor(for_hp$variable, c("rare","inter","core"))
  hp = ggplot(for_hp, aes(x = pattern, y = variable, fill = value)) + geom_tile(color = "black") +
    scale_fill_gradient(low = "white",high = "#023858", limits = c(0,num), name = "") + theme_minimal(base_size = 12) + ylab("") + xlab("")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), panel.grid = element_blank()) + 
    guides(fill = guide_colourbar(ticks = FALSE,  frame.colour = "black",  nbin = 47)) + theme(legend.position = "bottom")
  for_barplot$Var1 = factor(for_barplot$Var1, for_barplot$Var1) 
  legend = as_ggplot(get_legend(hp))
  hp = hp + theme(legend.position = "None")
  bp = ggplot(for_barplot, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity", fill = "#023858") +
    theme_classic(base_size = 12)+
    xlab("") + ylab("Genes")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), panel.grid = element_blank())+ 
    scale_y_continuous(expand = c(0.05,0.05)) + ggtitle("A")
  ## Analysing entirely rare variants
  p = p + ggtitle("B")
  grid.arrange(bp, hp,p, legend, layout_matrix = lay)
}


## genes which are missing only in one cluster
missing_in_one = gene_classes[which(gene_classes$total_presence == 46 & gene_classes$core == 46),]
missing_in_one = rownames(missing_in_one)
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

lay = rbind(c(1),
            c(1),
            c(1),
            c(2),
            3,
            3,
            3,
            3,
            4)
plot_for_num(47, p1, lay)



## this figure explains PC1
### find genes which are core and specific to a cluster
## genes which are missing only in one cluster
specific = gene_classes[which(gene_classes$total_presence == 1 & gene_classes$core == 1),]
specific = rownames(specific)
## for the ubiq genes, which clusters are they missing from
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
lay2 = rbind(c(1,3,3,3),
            c(1,3,3,3),
            c(1,3,3,3),
            c(2,3,3,3),
            c(4,4,4,NA))
plot_for_num(1, p2, lay2)



## relationship between number of depleted genes vs number of enriched genes

df = data.frame(cluster = summary_depleted$cluster,
                missing = summary_depleted$num_genes, specific = summary_specific$num_genes, stringsAsFactors = F)
lin = ggplot(df, aes(x = missing, y = specific)) + geom_point(alpha = 0.6, size = 1.4) +
  theme_classic(base_size = 14) +geom_text(aes(label=cluster),hjust=0, vjust=0) +
  xlab("Missing genes") + ylab("Specific genes")+ 
  geom_smooth(method='lm', col = 'black') + scale_x_continuous(expand = c(0.01,0,0,4))
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



