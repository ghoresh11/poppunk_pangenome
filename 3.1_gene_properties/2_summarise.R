library(ggplot2)
library(RColorBrewer)
library(reshape2)

setwd("/Users/gh11/poppunk_pangenome/3.1_gene_properties/")

variables = c("?","rare", "inter", "soft_core","core")
files = list.files(path = "results/", full.names = T)

cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes.csv", header = T,
                           comment.char = "", quote = "", stringsAsFactors = F, sep = ",")

## Read in and change to fit my needs
props = read.table(files[1], sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "")
cluster = strsplit(x = strsplit(files[1], split = "//", fixed = T)[[1]][2], split = "_", fixed = T)[[1]][1]
props = cbind(Cluster = rep(cluster, dim(props)[1]),
              Size = rep(cluster_sizes$Size[which(cluster_sizes$Cluster == cluster)], dim(props)[1]),
              props)
for (i in 2:length(files)){
  f = files[i]
  cluster = strsplit(x = strsplit(f, split = "//", fixed = T)[[1]][2], split = "_", fixed = T)[[1]][1]
  df = read.table(f, sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "")
  df = cbind(Cluster = rep(cluster, dim(df)[1]), 
             Size = rep(cluster_sizes$Size[which(cluster_sizes$Cluster == cluster)], dim(df)[1]),
             df, stringsAsFactors = F)
  props = rbind(props,df, stringsAsFactors = F)
}

props = props[-which(is.na(props$class)),] # not sure why this still happens
props$class = factor(props$class, variables)

### FUNCTIONS ###

plot_boxplot_per_prop <- function(name, props, column, path, xlabs, ylab, colors, log = F){
  colnames(props)[which(colnames(props) == name)] = "Class"
  colnames(props)[which(colnames(props) == column)] = "Var"
  if (log) {
    props$Var = log10(props$Var + 1)
  }
  p = ggplot(props, aes(x = Class, y = Var, fill = Class)) + geom_violin() +
    geom_boxplot(width = 0.1) + xlab("") +
    theme_classic(base_size = 16) + scale_fill_manual(values = colors, guide = F)+ 
    scale_x_discrete(labels = labs)+ theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    ylab(ylab)
  ggsave(p, file = paste(path, name, "_", column, ".pdf", sep = ""), width = 7, height = 4)
}

plot_size_to_genes_stratifies<- function(props, gene_type, property_column, property, probs, path){
  column = which(colnames(props) == property)
  quantiles = round(quantile(props[,column], probs = probs))
  vec = props[props[,property_column] == gene_type, column]
  curr_df = data.frame(Size = props$Size[props[,property_column] == gene_type], 
                       vals = rep("", length(vec)), stringsAsFactors = F)
  
  label = paste("<=", quantiles[1], sep = "")
  factor_vec = c(label)
  curr_df$vals[vec <= quantiles[1] ]= label
  
  for (i in 1:(length(quantiles)-1)) {
    label = paste(">", quantiles[i]," & <=", quantiles[i+1], sep = "")
    curr_df$vals[vec > quantiles[i] & vec <= quantiles[i+1] ] = label
    factor_vec = c(factor_vec, label)
  }
  
  label = paste(">", quantiles[length(quantiles)], sep = "")
  factor_vec = c(factor_vec, label)
  curr_df$vals[vec > quantiles[length(quantiles)] ]= label
  curr_df = data.frame(table(curr_df), stringsAsFactors = F)
  curr_df$Size = as.numeric(as.character(curr_df$Size))
  curr_df$vals = factor(curr_df$vals, factor_vec)
  p = ggplot(curr_df, aes(x = Size, y = Freq, fill = vals)) + geom_point(size = 3, shape = 21, color = "black") +
    theme_bw(base_size = 16) + scale_fill_brewer(palette = "Reds")  + xlab("Cluster Size") +
    ylab("Genes") + scale_x_continuous(limits = c(0, 1100)) +
    ggtitle(paste(gene_type, property, sep = "\n"))
  ggsave(plot =p,
         filename = paste(path, gene_type, "_", property, ".pdf", sep = ""), height = 4, width = 7)
}

### What's the difference between the different gene types for these properties
labs = c("Filtered by Roary", "Rare", "Intermediate", "Soft core", "Core")
colors = brewer.pal(n = 5, "Blues")
plot_boxplot_per_prop("class", props, "GC", "figures/", labs, "%GC", colors)
plot_boxplot_per_prop("class", props, "mean_length", "figures/", labs, "Gene length (log10(aa))", colors, log = T)
plot_boxplot_per_prop("class", props, "mean_pos", "figures/", labs, "Distance from edge (log10(bp))", colors, log = T)
plot_boxplot_per_prop("class", props, "mean_contig_length", "figures/",labs, "Contig length (log10(bp))", colors, log = T)

# ## How do these metrics correlate with each other? 
# ## DO smaller genes also have low GC content?NO
# for (i in c(7,8,10,12)) {
#   for (j in c(7,8,10,12)) {
#     if (i > j){
#       x = colnames(props)[i]
#       y = colnames(props)[j]
#       val = cor(props[,i], props[,j], method = "spearman")
#       print(paste(x, "-", y, ":", val, sep = " "))
#     }
#   }
# }
# ## they aren't correlated, that's why I see an increase in rare genes of size 100-200, but the contig lengths 
# ## are actually long and they're not at the edge of the contig


## Is the increase in rare genes driven by any of these properties?
column = which(colnames(props) == "class") 
path = "figures/props_per_class/"
for (var in c("rare", "inter", "soft_core", "core")){
  plot_size_to_genes_stratifies(props, var, column, "mean_length", seq(from = 0.1, to = 0.9, by = 0.2), path)
  plot_size_to_genes_stratifies(props, var, column,"mean_contig_length", seq(from = 0.05, to = 0.5, by = 0.1), path)
  plot_size_to_genes_stratifies(props, var, column, "mean_pos", seq(from = 0.05, to = 0.6, by = 0.15), path)
  plot_size_to_genes_stratifies(props, var, column,"GC", c(0.1, 0.9), path)
}

### Connection to COG category
cogs = read.table("../3_eggnog//cog_descs.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")
cogs = cogs[-which(cogs$Col == "black"),]
props = props[-which(props$class == "?" | props$COG %in% c("A","B")),]
props = props[which(unlist(props$COG) %in% cogs$COG),]
cogs = cogs[-which(!cogs$COG %in% unlist(props$COG) ),] # -> if the colours are wrong some COGs need to be removed
props$COG = unlist(as.character(props$COG))
props$COG = factor(props$COG, levels = cogs$COG)


labs = cogs$COG
colors = cogs$Col
plot_boxplot_per_prop("COG", props, "GC", "figures/", labs, "%GC", colors)
plot_boxplot_per_prop("COG", props, "mean_length", "figures/", labs, "Gene length (log10(aa))", colors, log = T)
plot_boxplot_per_prop("COG", props, "mean_pos", "figures/", labs, "Distance from edge (log10(bp))", colors, log = T)
plot_boxplot_per_prop("COG", props, "mean_contig_length", "figures/",labs, "Contig length (log10(bp))", colors, log = T)


## Is the increase in rare genes driven by any of these properties?
column = which(colnames(props) == "COG") 
for (cog in cogs$COG){
  path = "figures/props_per_COG/"
  plot_size_to_genes_stratifies(props, cog, column, "mean_length", seq(from = 0.1, to = 0.9, by = 0.2), path)
  plot_size_to_genes_stratifies(props, cog, column,"mean_contig_length", seq(from = 0.05, to = 0.5, by = 0.1), path)
  plot_size_to_genes_stratifies(props, cog, column, "mean_pos", seq(from = 0.05, to = 0.6, by = 0.15), path)
  plot_size_to_genes_stratifies(props, cog, column,"GC", c(0.1, 0.9), path)
}





