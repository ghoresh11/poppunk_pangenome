library(ggplot2)
library(RColorBrewer)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis/pre_comp_post/")


counts = read.table("pre_post_comparison.csv", sep = ",", header = T, stringsAsFactors = F)
counts= counts[-which(counts$cluster == 43),]
counts$type = factor(counts$type, c("core","soft_core", "inter",  "rare","Total"))
counts$pre_or_post = factor(counts$pre_or_post, c("pre", "post"))
## how does the number of genes change with this modification? 
## There's an increase in the number of genes, and some core genes are no longer core
ggplot(counts, aes( x = type, y = count, fill = pre_or_post)) +
  geom_boxplot() + theme_bw(base_size = 16) + xlab("") + ylab("Genes") +
  scale_fill_brewer(palette = "Set3", name = "", labels = c("No splitting", "Splitting"))

