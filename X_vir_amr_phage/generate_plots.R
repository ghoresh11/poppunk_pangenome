library(ggplot2)
library(RColorBrewer)


setwd("/Users/gh11/poppunk_pangenome/4_vir_amr_phage/")


cols = brewer.pal(n = 5, "Blues")[-1]
variable_order = c("rare", "inter", "soft_core", "core")
labs = c("Rare", "Intermediate", "Soft core", "Core")

file_names = dir(path = "gene_DB_hits/", pattern = "*.csv", full.names = T)
blast_hit_counts = do.call(rbind,lapply(file_names,read.table, sep = ",", header = T, 
                                        comment.char = "", stringsAsFactors = F, quote = ""))

blast_hit_counts = blast_hit_counts[-which(blast_hit_counts$DB %in% c("card", "phasta_bact")),]

databases = unique(blast_hit_counts$DB)

### see how this connects to the cluster size
## Does that number of Phage genes and AMR genes increase with the cluster size?
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes.csv", sep = ",",
                           header = T, stringsAsFactors = F)
blast_hit_counts = cbind(blast_hit_counts, 
                         size = cluster_sizes$Size[match(blast_hit_counts$Cluster, cluster_sizes$Cluster)])

outdir = "/Users/gh11/Submissions/ecoli_mobilome/figures/amr_vir_phage/"

for (db in databases) {
  curr_df = blast_hit_counts[which(blast_hit_counts$DB == db),]
  
  figure1 = paste(outdir, db, "_barplot.pdf", sep = "")
  figure2 = paste(outdir, db, "_size_genes.pdf", sep = "")
  figure3 = paste(outdir, db, "_size_genes_zoomed.pdf", sep = "")
  figure4 = paste(outdir, db, "_boxplot.pdf", sep = "")
  
  curr_df$Cluster = factor(curr_df$Cluster, 1:39)
  curr_df$Gene_class = factor(curr_df$Gene_class, variable_order)
  
  p = ggplot(curr_df, aes(x = Cluster, y = Count, fill = Gene_class)) + 
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = cols, labels = labs) + theme_classic(base_size = 16)  +
    scale_y_continuous(expand = c(0,0)) + ylab("Genes") + ggtitle(db)
  ggsave(plot = p, filename = figure1, width = 14, height = 4.5)
  
  p = ggplot(curr_df, aes( x = size, y = Count, fill = Gene_class)) + 
    geom_point(size = 3, alpha = 0.7, pch = 21, color = "black") +
  scale_fill_manual(values = cols) + ggtitle(db) + theme_bw(base_size = 16) +
    xlab("Cluster Size")
  ggsave(plot = p, filename = figure2, width = 6, height = 4)
  
  ## zoom_in
  zoomed = curr_df[which(curr_df$size<1000),]
  p = ggplot(zoomed, aes( x = size, y = Count, fill = Gene_class)) + 
    geom_point(size = 3, alpha = 0.7, pch = 21, color = "black") +
    scale_fill_manual(values = cols) + ggtitle(db) + theme_bw(base_size = 16) + 
    geom_smooth(method='lm',formula= y~x, se = T, aes(color = Gene_class)) +
    scale_color_manual(values = cols) +
    xlab("Cluster Size")
  ggsave(plot = p, filename = figure3, width = 6, height = 4)
  
  curr_df$Gene_class = factor(curr_df$Gene_class, rev(variable_order))
  p = ggplot(curr_df, aes( x= Gene_class, y = Count, fill = Gene_class)) + geom_violin()+
    geom_boxplot(width = 0.1) + scale_fill_manual(values = rev(cols)) +
    theme_bw(base_size = 16) + ggtitle(db)
  ggsave(plot = p, filename = figure4, width = 6, height = 4)
}





