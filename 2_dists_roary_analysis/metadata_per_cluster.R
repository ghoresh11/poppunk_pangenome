library(ggplot2)
library(RColorBrewer)


#setwd("/Users/gh11/poppunk_pangenome/dists_analysis/")

# metadata = read.table("metadata_per_cluster.csv", sep = "\t",
#                       header = T, stringsAsFactors = F, comment.char = "")
# 
# 
# clusters = unique(metadata$cluster)
# 
# ### ST
# ST = metadata[which(metadata$variable == "ST"),]
# ## count = data.frame(table(ST$value)) ## other than ST131,10,59 -> each ST is confined to one poppunk cluster
# for (c in clusters){
#   curr_df = ST[which(ST$cluster == c),]
#   curr_df$cluster = as.character(curr_df$cluster)
#   p = ggplot(curr_df, aes(x = cluster, y = count, fill = value)) + 
#     geom_bar(stat = "identity", width = 1)+ coord_polar(theta="y") +
#     theme_minimal(base_size = 16) +
#     xlab("") + ylab("")+
#     theme(axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           panel.grid  = element_blank()) + 
#     ggtitle(paste("Cluster: ",c, "\nVar: ", "ST", sep = ""))
#   ggsave(plot = p, 
#          filename = paste("/Users/gh11/Submissions/ecoli_mobilome/figures/clusters_metadata/ST/", c ,".png", sep = ""),
#          height = 4, width = 5)
# }
# 
# ### MASH
# # MASH = metadata[which(metadata$variable == "MASH"),]
# # MASH$cluster = factor(MASH$cluster, 1:39)
# # ggplot(MASH, aes( x = cluster, y = count, fill = value)) + geom_bar(stat = "identity", color = "black") + 
# #   theme_bw(base_size = 16) +
# #   xlab("Cluster") + ylab("% of isolates")
# 
# Pathotype = metadata[which(metadata$variable == "Pathotype"),]
# Pathotype$cluster = factor(Pathotype$cluster, 1:39)
# Pathotype$value = factor(Pathotype$value, c("nd", "epec", "etec","ehec","stec","expec", "eaec","epec/eaec","commensal"))
# cols = c("#dddddd", brewer.pal(n = 7, "Set2"), "#4c4cff")
# ggplot(Pathotype, aes( x = cluster, y = count, fill = value)) + geom_bar(stat = "identity", color = "black") +
#   scale_fill_manual(values = cols) + theme_bw(base_size = 16) +
#   xlab("Cluster") + ylab("% of isolates")
# 
# Continent = metadata[which(metadata$variable == "Continent"),]
# Continent$cluster = factor(Continent$cluster, 1:39)
# Continent$value = factor(Continent$value, c("Europe", "North America","Africa","Asia","Oceania", "South America","nd"))
# cols = c( brewer.pal(n = 6, "Set2"), "#dddddd")
# ggplot(Continent, aes( x = cluster, y = count, fill = value)) + geom_bar(stat = "identity", color = "black") +
#   theme_bw(base_size = 16) + scale_fill_manual(values = cols)+
#   xlab("Cluster") + ylab("% of isolates")
# 
# Isolation = metadata[which(metadata$variable == "Isolation"),]
# Isolation$cluster = factor(Isolation$cluster, 1:39)
# Isolation$value = factor(Isolation$value, c("feces","urine","blood","other/unknown"))
# cols = c( "#731d1d","#ffff66","#ff7f7f","#dddddd")
# ggplot(Isolation, aes( x = cluster, y = count, fill = value)) + 
#   geom_bar(stat = "identity", color = "black") +
#   theme_bw(base_size = 16) + scale_fill_manual(values = cols)+
#   xlab("Cluster") + ylab("% of isolates")
# 
# Publication = metadata[which(metadata$variable == "Publication"),]
# Publication$cluster = factor(Publication$cluster, 1:39)
# cols = c(brewer.pal(n=8, "Dark2"), brewer.pal(n = 8, "Set3"), brewer.pal(n = 8, "Set1"))
# ggplot(Publication, aes( x = cluster, y = count, fill = value)) + 
#   geom_bar(stat = "identity", color = "black") + scale_fill_manual(values = cols)+
#   theme_bw(base_size = 16) + 
#   xlab("Cluster") + ylab("% of isolates")
# 
# 
# Year = metadata[which(metadata$variable == "Year" & metadata$value != "nd"),]
# Year$cluster = factor(Year$cluster, 1:39)
# Year$value = as.numeric(Year$value)
# cols = rep(c(brewer.pal(n = 8, "Dark2")),5)
# shapes = c(0:25, 0:25)
# ggplot(Year, aes(x = value, y = count, color = cluster, shape = cluster)) + geom_point(size = 3, alpha = 0.7)+
#   scale_color_manual(values = cols) + scale_shape_manual(values = shapes) +
#   theme_bw(base_size = 16) + xlab("Year") + ylab("% of isolates")
# 

### Look at the metadata per cluster, stratified

loc = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/metadata_within_cluster"
plots = paste(loc, "plots/", sep = "/")
#loc = "/Users/gh11/poppunk_pangenome/dists_analysis/"
#plots = loc
setwd(loc)

num_tests = 546 ## I had a look how many files were created
required_pval = 0.05/ num_tests ## bonferroni corrected
## The hypothesis is that x<y (i.e. distance within a group is smaller than between groups)

for (cluster in 1:39){
  curr_cluster = read.table(paste(cluster,"_metadata_within_cluster.csv", sep = ""), sep = ",",
                            stringsAsFactors = F, header = T, comment.char = "", quote = "")
  same = rep("Different", dim(curr_cluster)[1])
  curr_cluster = cbind(curr_cluster, same)
  columns = c(4:10)
  for (c in columns){
    curr_cluster$same = rep("Different", dim(curr_cluster)[1])
    curr_cluster$same[which(curr_cluster[,c] == curr_cluster[,c+7])] = "Same"
    title = paste("Cluster:", cluster )
    name = strsplit(x = colnames(curr_cluster)[c], split = "_", fixed = T)[[1]][1]
    
    if (length(which(curr_cluster$same == "Different")) == 0  || length(which(curr_cluster$same == "Same")) == 0 ) {
      next
    }
    
    core_test = wilcox.test(x = curr_cluster$Core_dist[which(curr_cluster$same == "Same")],
                            y = curr_cluster$Core_dist[which(curr_cluster$same == "Different")], alternative = "less")
    
    acc_test = wilcox.test(x = curr_cluster$Acc_dist[which(curr_cluster$same == "Same")],
                           y = curr_cluster$Acc_dist[which(curr_cluster$same == "Different")], alternative = "less")
    
    if (core_test$p.value < required_pval || acc_test$p.value < required_pval) {
      p = ggplot(curr_cluster, aes(x = Core_dist, y = Acc_dist, color = same)) + 
        geom_point( size = 3, alpha = 0.7) +
        theme_bw(base_size = 16) + 
        scale_color_manual(values = c("#909090","#7c26cb")) +
        xlab("Core distance") + ylab("Accessory distance") + ggtitle(paste(title, name, sep = "\n")) 
      ggsave(p, height = 5, width = 6,
               filename = paste(plots, cluster, "_", name, ".pdf", sep = ""))
      
      p = ggplot(curr_cluster, aes(x = same, y = Core_dist)) + geom_violin() +
        geom_boxplot(width = 0.1)  +  ggtitle(paste(title, "\n",name, "\n",  core_test$p.value , sep = "")) +
        ylab("Core distance") + xlab("") + theme_classic(base_size = 16)
      ggsave(p, height = 4, width = 4, 
             filename = paste(plots, cluster, "_", name, "_core_dists.pdf", sep = ""))
      p = ggplot(curr_cluster, aes(x = same, y = Acc_dist)) + geom_violin() +
        geom_boxplot(width = 0.1)  + ggtitle(paste(title, "\n",name, "\n",  acc_test$p.value , sep = "")) +
        ylab("Accessory distance") + xlab("") + theme_classic(base_size = 16)
      ggsave(p, height = 4, width = 4, 
             filename = paste(plots, cluster, "_", name, "_acc_dists.pdf", sep = ""))
    }
  }
}


df = data.frame(gene_systems = c(1, 60000), genomes = c(300, 10000), type = c("ta", "ecoli"))
ggplot(df, aes(x = gene_systems, y = genomes, color = type)) + geom_point(size = 6, pch = 13, stroke = 2)  + 
  scale_y_continuous(limits = c(0, 10000)) + theme_bw(base_size = 20) + scale_alpha_continuous(limits = c(0, 60000)) + 
  xlab("Genes") + ylab("Genomes") + scale_color_manual(values = c("red", "#d3d3d3"), guide = F)
