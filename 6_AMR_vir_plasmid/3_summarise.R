library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(reshape2)
library(ape)

######## General things ###########
setwd("/Users/gh11/poppunk_pangenome/6_AMR_vir_plasmid/")

orig_md = read.table("/Users/gh11/e_colis/FILTERED_MD_FINAL_ALL.tab", sep = "\t",
                     header = T, comment.char = "", quote = "", stringsAsFactors = F)
md_names = as.character(sapply(sapply(sapply(sapply(orig_md$New_annot_loc, strsplit, split = "/", fixed = T),tail, n = 1), 
                                      strsplit, split = ".",fixed = T), head, n=1))
rownames(orig_md) = md_names

graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv",
                      sep = ",", comment.char = "", header = T, stringsAsFactors = F)


graphics = graphics[order(as.numeric(graphics$Cluster), decreasing = F),]

### plotting heatmaps of genes along the tree
tree = read.tree("smaller_tree/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "45")
tree_md = read.table("smaller_tree/smaller_tree.csv", sep = ",", comment.char = "", stringsAsFactors = F,
                     quote = "", header = T)

plot(tree)

############## FUNCTIONS #############

plot_boxplot <- function(md, vec, breaks, ylab){
  md = cbind(md, vec)
  pairwise.wilcox.test(x = md$vec, g = md$Poppunk_cluster, p.adjust.method = "fdr")
  md$Poppunk_cluster = factor(md$Poppunk_cluster, graphics$Cluster)
  p = ggplot(md, aes(x = Poppunk_cluster, y = vec, shape = Poppunk_cluster, color = Poppunk_cluster)) + 
    geom_jitter(width = 0.3, height = 0.2, size = 1.2, stroke = 1.2, alpha = 0.5) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA,) +
    scale_color_manual(values = graphics$Color) + scale_shape_manual(values = graphics$Shape) +
    theme_classic(base_size = 14) + theme(legend.position = "None") + 
    scale_y_continuous(breaks = breaks) + ylab(ylab) + xlab("Cluster") +
    coord_flip()
  return(p)
}

calculate_stats <- function(vec, md) {
  rs = sample(x  = unique(md$Poppunk_cluster),size = 10000, replace = T)
  res = rep(0, 10000)
  for (index in 1:10000) {
    cluster = rs[index]
    curr = vec[which(md$Poppunk_cluster == cluster)]
    res[index] = sample(size = 1, x = curr)
  } 
  signif = c()
  for (cluster in unique(md$Poppunk_cluster)) {
    curr_test = wilcox.test(res, vec[which(md$Poppunk_cluster == cluster)], alternative = "less")
    if (curr_test$p.value<(0.05/47)) {
      signif = c(signif, cluster)
    }
  }
  return(signif)
}

read_one_file <- function(filename, ylab, desc_file){
  df = read.table(filename, sep = ",", header = T, comment.char = "",
                  quote = "", stringsAsFactors = F, row.names = 1)
  df_no_trunc = df
  df_no_trunc[df_no_trunc=="1*"] = 0
  df_no_trunc = data.frame(lapply(df_no_trunc, as.numeric))
  rownames(df_no_trunc) = row.names(df)
  ## get the metadata file in the same order as the plasmids file
  md = orig_md[match(rownames(df_no_trunc), rownames(orig_md)),]
  num = rowSums(df_no_trunc)
  add_orig_md = num[match(rownames(orig_md), rownames(df_no_trunc))]
  orig_md = cbind(orig_md, add_orig_md)
  
  descs = read.table(desc_file, sep = ",", comment.char = "", quote = "", stringsAsFactors = F, header = T)

  sigs = calculate_stats(vec = num, md = md)
  
  empty_descs = descs[-which(is.na(descs)),]
  res = cbind(data.frame(cluster = character(0), freq = character(0), stringsAsFactors = F), empty_descs)
  
  for (s in sigs) {
    curr_df = df_no_trunc[which(md$Poppunk_cluster == s),]
    genes = colSums(curr_df)/dim(curr_df)[1]
    all_genes = sapply(aggregate(df_no_trunc, list(md$Poppunk_cluster), mean)[,-1], mean)
    genes = genes[which(genes > all_genes + 0.1)]
    curr = descs[match(names(genes), descs$identifier),]
    curr = cbind(cluster = rep(s, length(genes)), freq = genes, curr)
    res = rbind(res, curr)
  }
  p = plot_boxplot(md, num, seq(from=0, to=max(num),by=1), ylab)
  write.table(x = res, file = paste("results/",ylab, "_signif.csv", sep = ""), sep = ",", row.names = F, col.names = T, quote = F)
  return(list(orig_md, sigs, p))
}

get_gene_for_one_cluster <- function(cluster, md, df){
  df = df[which(md$Poppunk_cluster == cluster),]
  df[df=="1*"] = 0 ## here: decide what to do with proteins which appear to be truncated (maybe want to look at them)
  df = data.frame(lapply(df, as.numeric))
  freq = colSums(df, na.rm = T)/dim(df)[1]
  # df[df=="1"] = 0
  # df[is.na(df)] = 1
  # freq_trunc = colSums(df, na.rm = T)/dim(df)[1]
  res = data.frame(gene = rep(names(freq), 1),
                   cluster = rep(cluster,length(freq)),
                   #  type = c(rep("full", length(freq)), rep("trunc",length(freq))),
                   freq = c(freq), stringsAsFactors = F)
  return(res)
}

plot_on_tree <- function(filename, desc_file, tree, signif){
  df = read.table(filename, sep = ",", header = T, comment.char = "",
                  quote = "", stringsAsFactors = F, row.names = 1)
  md = orig_md[match(rownames(df), rownames(orig_md)),]
  dfs = lapply(X = unique(md$Poppunk_cluster), FUN = get_gene_for_one_cluster, md = md, df = df)
  res = do.call(rbind, dfs)
  res = dcast(res, cluster ~ gene, value.var = "freq")
  rownames(res) = res[,1]
  res = res[,-1]
  res = res[match(tree$tip.label, rownames(res)),]
  
  ## change the order of the columns
  res = res[,order(colSums(res), decreasing = T)]
  res = res[,-which(sapply(FUN = max, X = res)<0.1)] ## to remove very rare genes
  
  descs = read.table(desc_file, sep = ",", header = T, stringsAsFactors = F, quote = "")
  new_labs = descs$gene[match(colnames(res), descs$identifier)]
  labs = c()
  for (l in new_labs) {
    while (l %in% labs) {
      l = paste(l, "^", sep = "")
    }
    labs = c(labs,l)
  }
  colnames(res) = labs
  col.order <- rev(hclust(dist(t(res)))$order)
  res = res[,col.order]

  tip_labels = tree$tip.label
  for (i in 1:length(tip_labels)) {
    if (tip_labels[i] %in% signif){
      tip_labels[i] = paste(tip_labels[i], "*", sep = "")
    }
  }
  tree$tip.label = tip_labels
  rownames(res) = tip_labels
  p = ggtree(tree)  +
    theme(legend.position="right") + geom_tiplab(align = T)
  
  
  p2 = gheatmap(p = p, data = res, offset = 0.1, width=5, color = "black", 
                high = "#023858", low = "#fff7fb", colnames_angle = 90, colnames = T, colnames_position = "bottom",
                font.size = 3)
  return(p2)
}

get_median_per_cluster <- function(filename,lab){
  df = read.table(filename, sep = ",", header = T, comment.char = "",
                  quote = "", stringsAsFactors = F, row.names = 1)
  df[df == "1*"] = 0
  md = orig_md[match(rownames(df), rownames(orig_md)),]
  df = data.frame(lapply(df, as.numeric))
  names = unique(md$Poppunk_cluster)
  res = c()
  res2 = c()
  for (clstr in names) {
    num1 = median(rowSums(df[which(md$Poppunk_cluster == clstr),]))
    num2 = length(which(colSums(df[which(md$Poppunk_cluster == clstr),])/length(which(md$Poppunk_cluster == clstr)) > 0.1))
    res = c(res, num1)
    res2 = c(res2, num2)
  }
  df_for_plot = data.frame(median = res,
                           num_genes_total = res2,
                           name = names, stringsAsFactors = F)
  print(ggplot(df_for_plot, aes(x = median, y = num_genes_total, label = name)) + geom_text(position=position_jitter(width=0.4,height=0.4))+
    theme_bw(base_size = 14) + xlab(paste("Median", lab, "per isolate")) + ylab("Total genes in PopPUNK cluster"))
  names(res) = names
  return(res)
}


## without truncated

### AMR genes
amr_file = "results/resfinder.csv"
res = read_one_file(amr_file, "AMR_genes", "DBs/break_names/resfinder_genes_fixed.csv")
orig_md = res[[1]]
signif = res[[2]]
A = res[[3]]
colnames(orig_md)[dim(orig_md)[2]] = "amr"
B = plot_on_tree("results/resfinder.csv","DBs/break_names/resfinder_genes_fixed.csv", tree, signif )
med_amr = get_median_per_cluster("results/resfinder.csv", "resistance genes")
grid.arrange(A, B, layout_matrix = lay)


### virulence genes
vir_file = "results/virulence.csv"
res = read_one_file(vir_file, "virulence", "DBs/break_names/virulence.csv")
orig_md = res[[1]]
signif = res[[2]]
colnames(orig_md)[dim(orig_md)[2]] = "vir"
plot_on_tree(vir_file,"DBs/break_names/virulence.csv", tree, signif )
med_vir = get_median_per_cluster("results/virulence.csv", "virulence genes")

### plasmids
res = read_one_file("results/plasmid.csv", "plasmid", "DBs/break_names/virulence.csv")
orig_md = res[[1]]
signif = res[[2]]
colnames(orig_md)[dim(orig_md)[2]] = "plasmid"
plot_on_tree("results/plasmid.csv","DBs/break_names/plasmid.csv", tree, signif)
med_plasmid = get_median_per_cluster("results/plasmid.csv", "plasmid replicons")

### plot a 3-way scatter plot to see the relationships between having any of these genes



med = data.frame(cluster = names(med_amr),
                 amr = med_amr,
                 vir= med_vir, plasmid = med_amr, stringsAsFactors = F)
A = ggplot(med, aes(x = med_plasmid, y = med_amr, label = cluster)) + geom_text(position=position_jitter(width=0.4,height=0.4))+
  theme_bw(base_size = 14) + xlab("Median plasmid replicons per isolate") + ylab("Median resistance genes per isolate")

B = ggplot(med, aes(x = med_plasmid, y = med_vir, label = cluster)) + geom_text(position=position_jitter(width=0.4,height=0.4))+
  theme_bw(base_size = 14) + xlab("Median plasmid replicons per isolate") + ylab("Median virulence genes per isolate")

C = ggplot(med, aes(x = med_amr, y = med_vir, label = cluster)) + geom_text(position=position_jitter(width=0.4,height=0.4))+
  theme_bw(base_size = 14) + xlab("Median resistance genes per isolate") + ylab("Median virulence genes per isolate")


## rewriting the metadata nicely:
orig_md = cbind(orig_md, phylogroup = graphics$Phylogroup[match(orig_md$Poppunk_cluster, graphics$Cluster)])
orig_md = orig_md[c(1:6, 9, 21:23, 34:36, 38:42)]
write.table(orig_md, file = "../SELECTED_METADATA.csv", sep = "\t", col.names = T, row.names = F, quote = F)



### understanding the which genes are stable in which clusters
## and how many genes are stable in each cluster and which


# orig_md$Poppunk_cluster = factor(orig_md$Poppunk_cluster, graphics$Cluster)
# 
# 
# ## change Y in this plot to look at plasmids/amr/virulence/stx
# orig_md$Isolation[which(!orig_md$Isolation %in% c("urine","feces","blood"))] = "other"
# ggplot(orig_md, aes(x = phylogroup, y = plasmids, color = Poppunk_cluster, shape = Poppunk_cluster)) + 
#   geom_jitter(width = 0.3, height = 0.2, size = 1.2, stroke = 1.2, alpha = 0.5) +
#   scale_color_manual(values = graphics$Color) + scale_shape_manual(values = graphics$Shape) +
#   theme_classic(base_size = 14) + theme(legend.position = "None") +  
#   scale_y_continuous(breaks = seq(from = 0, to = 28, by = 1))
