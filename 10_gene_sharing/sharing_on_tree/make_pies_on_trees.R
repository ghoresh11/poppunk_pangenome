library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(ape)


### all functions required to draw pie plots on a phylogenetic tree to see which clades share genes
## FUNCTIONS

setwd("/Users/gh11/poppunk_pangenome/10_gene_sharing/sharing_on_tree/")

find_specific_core_genes <- function(curr_clusters){
  ## find all the genes that are specific to this clade
  opposite_clade = genes_of_interest_freqs[,which(!clusters %in% curr_clusters)]
  missing_in_rest = row.names(opposite_clade)[which(apply(opposite_clade, 1, FUN = sum) == 0 )]
  return(missing_in_rest)
}

### MAIN

### step 1: load all relevant files
freqs = read.table("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv",
                   sep = ",",comment.char = "", stringsAsFactors = F, header = F, quote = "", row.names = 1)
freqs = freqs[,-47]
clusters = as.character(freqs[1,])
freqs = freqs[-1,]

cluster_graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv", sep =",",
                              comment.char = "", stringsAsFactors = F, header = T)

tree = read.tree("/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/tree_for_treeseg.nwk")
all_subtrees = subtrees(tree)

classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep = "\t",
                            header = T, stringsAsFactors = F, comment.char = "", quote = "")

colours = read.table("../../5_classify_genes/colours_v2.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")
### run on each gene class separately

gene_classes = unique(classification$fill)

## for each gene class save the percentage of genes which are
## present at the root of the tree against the mean number of clusters in which they were found
root_df = data.frame(gene_classes,
                     percentage = rep(0, length(gene_classes)),
                     mean_presence = rep(0, length(gene_classes)),
                     sd_presence = rep(0, length(gene_classes)),
                     stringsAsFactors = F)

for (gene_class in gene_classes) {
  col1 = colours$Colour[which(colours$Class == gene_class)]
  col = "#d3d3d3"
  
  genes_of_interest_freqs = freqs[rownames(freqs) %in% classification$gene[classification$fill == gene_class],]
  
  fileConn = paste(gene_class, "_itol_pie.txt",sep = "")
  file_open =  file(fileConn)
  writeLines(c("DATASET_PIECHART","SEPARATOR COMMA",
               "DATASET_LABEL,pies",paste("COLOR,",col1,sep = ""),paste("FIELD_COLORS,",col, sep =""),"FIELD_LABELS,f1","BORDER_WIDTH,2",
               "DATA"), file_open)
  close(file_open)
  
  label_file = paste(gene_class, "_itol_labels.txt",sep = "")
  file_open =  file(label_file)
  writeLines(c("DATASET_TEXT","SEPARATOR COMMA","DATASET_LABEL,example text dataset",paste("COLOR,",col,sep=""),
               "DATA"), file_open)
  close(file_open)
  complete = c()
  for (i in length(all_subtrees):2) {
    curr = all_subtrees[[i]]
    curr_specifics = find_specific_core_genes(curr$tip.label)
    if (length(complete) > 0){
      curr_specifics = curr_specifics[which(!curr_specifics %in% complete)]
    }
    complete = c(complete, curr_specifics)
    
    curr_specifics_count = dim(genes_of_interest_freqs)[1]
    write(paste("INT",i+length(all_subtrees)+1, ",0.3," ,curr_specifics_count, ",",curr_specifics_count,sep = ""), file = fileConn, append = T)
    if (curr_specifics_count > 10){
      write(paste("INT",i+length(all_subtrees)+1,",", curr_specifics_count, ",0.3,black,bold,2,270",sep = ""), file = label_file,append = T)
    }
  }
  
  root_df$percentage[which(root_df$gene_classes == gene_class)] = round( (dim(genes_of_interest_freqs)[1] - length(complete)) / dim(genes_of_interest_freqs)[1]*100, digits = 0)
  curr = genes_of_interest_freqs[which(!rownames(genes_of_interest_freqs) %in% complete),]
  curr = data.frame(sapply(curr, FUN = ceiling), stringsAsFactors = F)
  root_df$mean_presence[which(root_df$gene_classes == gene_class)]  = mean(apply(curr, 1, FUN = sum))
  root_df$sd_presence[which(root_df$gene_classes == gene_class)]  = sd(apply(curr, 1, FUN = sum))
}


## view the results of root df
root_df$gene_classes = factor(root_df$gene_classes, colours$Class)
ggplot(root_df, aes(y = mean_presence, x = percentage, fill = gene_classes)) + geom_point(pch = 21, size = 4) +
  scale_fill_manual(values = colours$Colour, guide = F) + theme_bw(base_size = 14)


