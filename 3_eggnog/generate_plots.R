library(ggplot2)
library(RColorBrewer)

setwd("/Users/gh11/poppunk_pangenome/3_eggnog/")


cogs = read.table("eggnog_cog_summary.csv", sep = ",", header = T,
                  stringsAsFactors = F, quote = "")

## connection between cluster size and category changing
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes.csv", sep = ",",
                           header = T, stringsAsFactors = F)

cogs = cbind(cogs, size = cluster_sizes$Size[match(cogs$cluster, cluster_sizes$Cluster)])

gene_type_cols = brewer.pal(n = 5, "Blues")[-1]
variable_order = c("rare", "inter", "soft_core", "core")
labs = c("Rare", "Intermediate", "Soft core", "Core")


CELLULAR_PROCESSES_AND_SIGNALING = c("D", "M", "N", "O", "T", "U", "V", "W", "Z")
INFORMATION_STORAGE_AND_PROCESSING = c("A", "B", "J", "K", "L")
METABOLISM = c("C", "E", "F", "G", "H", "I", "P", "Q")
POORLY_CHARACTERIZED = c("S", "?")

descs = read.table("cog_descs.csv", sep = ",", header = F, stringsAsFactors = F)

all =  c(CELLULAR_PROCESSES_AND_SIGNALING, INFORMATION_STORAGE_AND_PROCESSING, METABOLISM, "X", POORLY_CHARACTERIZED)
cols = c(
  brewer.pal(n = 8, "Blues")[c(7,8)],
  brewer.pal(n = 8 , "BuPu")[-1],
  brewer.pal(n = length(INFORMATION_STORAGE_AND_PROCESSING) + 1, "Greens")[-1],
  brewer.pal(n = length(METABOLISM) + 1, "Oranges")[-1],
  brewer.pal(n = 8, "Greys")[2:(length(POORLY_CHARACTERIZED)+2)]
)
names(cols) = all


### barplots of all rare, inter, core and soft-core genes
## for all the cluster seperately
for (gene_class in variable_order){
  rare = cogs[which(cogs$gene_type == gene_class),]
  for (c in unique(cogs$cog_cat)){
    if (! c %in% rare$cog_cat) {
      rare = rbind(rare, data.frame(cluster = 1,
                                    gene_type = gene_class,
                                    cog_cat = c,
                                    count = 0))
    }
  }
  rare$cluster = factor(rare$cluster, 1:39)
  rare$cog_cat = factor(rare$cog_cat, all)
  p = ggplot(rare, aes(x = cluster, y = count, fill = cog_cat)) + geom_bar(stat = "identity", color = "black", size = 0.1) +
    scale_fill_manual(values = cols) + theme_classic(base_size = 16) + ggtitle(gene_class)
  print(p)
}

figures_out = "/Users/gh11/Submissions/ecoli_mobilome/figures/functions/per_cog_category/"
## Look at how each COG category changes in count across all clusters
for (c in unique(cogs$cog_cat)){
  curr_cat = cogs[which(cogs$cog_cat == c),]
  curr_cat$gene_type = factor(curr_cat$gene_type, rev(variable_order))
  if (c != "?"){
    desc = paste(c, descs$V2[which(descs$V1 == c)], sep = ": ")
  } else {
    desc = "Unknown"
  }
  p = ggplot(curr_cat, aes(x = gene_type, y = count)) + geom_point(size = 3, alpha = 0.6) +
    theme_bw(base_size = 16) + ggtitle(desc) + xlab("")
  out = paste(figures_out, c, ".pdf", sep = "")
  ggsave(plot = p, filename = out, height = 4.5, width = 5)
}

barplot_df = data.frame(cog_cat = character(0), gene_type = character(0), count = numeric(0), sd = numeric(0), stringsAsFactors = F)
for (c in unique(cogs$cog_cat)) {
  for (gene_type in variable_order) {
    mean_genes = mean(cogs$count[which(cogs$gene_type == gene_type & cogs$cog_cat == c)])
    sd_genes = sd(cogs$count[which(cogs$gene_type == gene_type & cogs$cog_cat == c)])
    if (is.nan(mean_genes)) {
      mean_genes = 0
      sd_genes = 0
    }
    barplot_df = rbind(barplot_df,
                       data.frame(cog_cat = c, gene_type = gene_type, count = mean_genes, sd = sd_genes, stringsAsFactors = F))
  }
}

barplot_df$gene_type = factor(barplot_df$gene_type, variable_order)
barplot_df$cog_cat = factor(barplot_df$cog_cat, all)
barplot_df = barplot_df[-which(!barplot_df$cog_cat %in% c("S","?")),]
ggplot(barplot_df, aes(fill=gene_type, y=count, x=cog_cat)) + 
  geom_bar(stat="identity", color="black", size = 0.1,
           position=position_dodge()) + scale_fill_manual(values = gene_type_cols)+
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.1,color = "black",
                position=position_dodge(.9)) + theme_classic(base_size = 16) + ggtitle("Unknown functions") 



core_cats = c("K", "D", "M",  "N","U","L", "E","F","G","H","I","C","J","O","Q","T","P") ## Anything described
core_cats = c( "?", "S") # Mix
core_cogs = cogs[which(cogs$cog_cat %in% core_cats),]
curr_cols = cols[which(names(cols) %in% core_cats)]
core_cogs$cog_cat = factor(core_cogs$cog_cat, all)
core_cogs$cluster = factor(core_cogs$cluster, 1:39)
ggplot(core_cogs, aes(x = cluster, y = count, fill = cog_cat)) + 
  geom_bar(stat = "identity", color = NA) +
  ggtitle("Unknown functions") + scale_fill_manual(values = curr_cols, name = "COG category") +
  theme_classic(base_size = 16) + scale_y_continuous(expand = c(0,0)) + theme(legend.position="bottom") + 
  xlab("Cluster") + ylab("Genes")

ggplot(cogs, aes(x = size, y = count, fill = cog_cat)) + geom_point(size = 2, alpha = 0.7, pch = 21, color = "black") +
  scale_fill_manual(values = cols) + theme_bw(base_size = 16)
## zoom
zoom = cogs[-which(cogs$size>1000),]
ggplot(zoom, aes(x = size, y = count, fill = cog_cat)) + geom_point(size = 3, alpha = 0.7, pch = 21, color = "black") +
  scale_fill_manual(values = cols) + theme_bw(base_size = 16) 


