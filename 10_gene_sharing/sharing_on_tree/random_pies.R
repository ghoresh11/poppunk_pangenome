library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ape)


### all functions required to draw pie plots on a phylogenetic tree to see which clades share genes
## FUNCTIONS


#setwd("/Users/gh11/poppunk_pangenome/10_gene_sharing/sharing_on_tree/")
setwd("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/random_sharing")

find_specific_core_genes <- function(curr_clusters){
  ## find all the genes that are specific to this clade
  opposite_clade = genes_of_interest_freqs[,which(!clusters %in% curr_clusters)]
  missing_in_rest = row.names(opposite_clade)[which(apply(opposite_clade, 1, FUN = sum) == 0 )]
  return(missing_in_rest)
}


### MAIN

### step 1: load all relevant files
# freqs = read.table("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv",
#                    sep = ",",comment.char = "", stringsAsFactors = F, header = F, quote = "", row.names = 1)
freqs = read.table("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/correction/corrected_freqs.csv",
                   sep = ",",comment.char = "", stringsAsFactors = F, header = F, quote = "", row.names = 1)


freqs = freqs[,-47]
clusters = as.character(freqs[1,])
freqs = freqs[-1,]

# cluster_graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv", sep =",",
#                               comment.char = "", stringsAsFactors = F, header = T)
cluster_graphics = read.table("cluster_graphics.csv", sep =",",
                              comment.char = "", stringsAsFactors = F, header = T)

#tree = read.tree("/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/tree_for_treeseg.nwk")
tree = read.tree("tree_for_treeseg.nwk")
all_subtrees = subtrees(tree)
# 
# classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep = "\t",
#                             header = T, stringsAsFactors = F, comment.char = "", quote = "")
classification = read.table("classification_v2.csv", sep = "\t",
                            header = T, stringsAsFactors = F, comment.char = "", quote = "")


# colours = read.table("../../5_classify_genes/colours_v2.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")
colours = read.table("colours_v2.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")

### run on each gene class separately

gene_classes = unique(classification$fill)


precentage_explained = data.frame(gene_classes,
                                  num_genes = rep(0, length(gene_classes)),
                                  num_explained = rep(0, length(gene_classes)),
                                  stringsAsFactors = F)
gene_classes = c(gene_classes)
num_reps = 1000

for (gene_class in gene_classes) {
  if (gene_class %in% c("Population core", "Cluster specific rare", "Cluster specific core","Cluster specific intermediate")) {next}
  print(gene_class)
  col1 = colours$Colour[which(colours$Class == gene_class)]
  genes_of_interest_freqs = freqs[rownames(freqs) %in% classification$gene[
    which(classification$fill == gene_class & classification$total_presence < 41)],]
  
  df_out_spec = data.frame(matrix(0, nrow = length(all_subtrees), ncol = num_reps+1))
  df_out_missing = data.frame(matrix(0, nrow = length(all_subtrees), ncol = num_reps+1))
  rownames(df_out_spec) = paste("INT", 1:length(all_subtrees) + length(all_subtrees) + 1, sep = "")
  rownames(df_out_missing) = paste("INT", 1:length(all_subtrees) + length(all_subtrees) + 1, sep = "")
  
  for (i in length(all_subtrees):1) {
    curr = all_subtrees[[i]]
    df_out_spec[,1][i] = length(curr$tip.label)
    df_out_missing[,1][i] = length(curr$tip.label)
  }
  df_out_spec = df_out_spec[order(df_out_spec[,1], decreasing = T),]
  df_out_missing = df_out_missing[order(df_out_missing[,1], decreasing = T),]
  
  for (rep in 1:num_reps) {
    genes_of_interest_freqs = genes_of_interest_freqs[,sample(1:47,47, replace = F)]
    complete = c()
    for (i in dim(df_out_spec)[1]:2) {
      curr_subtree = as.numeric(gsub(pattern = "INT", x = row.names(df_out_spec)[i], replacement = "")) - length(all_subtrees) - 1
      curr = all_subtrees[[curr_subtree]]
      curr_specifics = find_specific_core_genes(curr$tip.label)
      if (length(complete) > 0){
        curr_specifics = curr_specifics[which(!curr_specifics %in% complete)]
      }
      complete = c(complete, curr_specifics)
      curr_specifics_count = length(curr_specifics)
      df_out_spec[i,rep+1] = curr_specifics_count
    }
    
    for (i in 1:dim(df_out_missing)[1]) {
      curr_subtree = as.numeric(gsub(pattern = "INT", x = row.names(df_out_missing)[i], replacement = "")) - length(all_subtrees) - 1
      curr = all_subtrees[[curr_subtree]]
      ## find genes which are entirely missing in this subclade
      curr_clade = genes_of_interest_freqs[,which(clusters %in% curr$tip.label)]
      missing_in_clade = row.names(curr_clade)[which(apply(curr_clade, 1, FUN = sum) == 0 )]
      
      opposite_clade = genes_of_interest_freqs[,which(!clusters %in% curr$tip.label)]
      common_in_opposite = c()
      if (dim(opposite_clade)[1] > 0){
        common_in_opposite = row.names(opposite_clade)[which((rowSums(opposite_clade != 0) / dim(opposite_clade)[2])> 0.8)] #
      }
      missing_in_clade = missing_in_clade[which(missing_in_clade %in% common_in_opposite)]
      
      if (length(complete) > 0){
        missing_in_clade = missing_in_clade[which(!missing_in_clade %in% complete)] ## already missing in a higher up ancestor
      }
      complete = c(complete, missing_in_clade)
      curr_missing_count = length(missing_in_clade)
      df_out_missing[i,rep+1]  = curr_missing_count
    }
    
  }
  
  df_out = data.frame(id = row.names(df_out_missing),
                      position = rep(0.3, length(all_subtrees)),
                      size = rep(0, length(all_subtrees)),
                      specific = round(apply(X = df_out_spec[,-1], 1, FUN = mean)),
                      missing =  round(apply(X = df_out_missing[,-1], 1, FUN = mean)), 
                      stringsAsFactors = F)
  df_out$size = df_out$missing + df_out$specific
  
  ## what percentage of genes can be explained by the phylogeny?
  explained = sum(df_out$size)/dim(genes_of_interest_freqs)[1] ## compare to measured
  print(explained)
  ## rewrite the labels so they have asterisks if it's significant
  label_file = paste(gene_class, "_itol_labels_signif.txt",sep = "")
  file_open =  file(label_file)
  writeLines(c("DATASET_TEXT","SEPARATOR COMMA","DATASET_LABEL,labels",paste("COLOR,",col1,sep=""),
               "DATA"), file_open)
  close(file_open)
  
  labels = read.table(paste(gene_class, "_itol_labels.txt",sep=""), sep = ",", skip=5, comment.char = "", stringsAsFactors = F)
  labels$V2 = as.character(labels$V2)
  ## what if the p-value for each internal node?
  measured = read.table(paste(gene_class, "_itol_pie.txt",sep=""), sep = ",", skip=8, comment.char = "", stringsAsFactors = F)
  measured_spec = measured$V4
  measured_missing = measured$V5
  
  df_out_spec = df_out_spec[,-1]
  df_out_missing = df_out_missing[,-1]
  
  write.table(df_out_spec, paste("random_results/",gene_class,"_specific.csv", sep = ""), col.names = F, row.names = T, sep = ",", quote = F)
  write.table(df_out_missing, paste("random_results/",gene_class,"_missing.csv", sep = ""), col.names = F, row.names = T, sep = ",", quote = F)
  
  for (i in 1:length(measured_spec)) {
    pval = 1-length(which(df_out_spec[i,] < measured_spec[i]))/num_reps
    if (pval < 0.05/47){
      index = which(labels$V1 == rownames(df_out_spec)[i])
      if (length(index) == 0) { next }
      res = strsplit(x = labels$V2[index], split = " ", fixed = T)[[1]]
      if (length(res)>1) {
        labels$V2[index] = paste(res[1],"* ",res[2], sep = "")
      } else {
        labels$V2[index] =  paste(res[1],"*", sep = "")
      }
    }
    pval = 1-length(which(df_out_missing[i,] < measured_missing[i]))/num_reps
    if (pval < 0.05/47){
      index = which(labels$V1 == rownames(df_out_spec)[i])
      if (length(index) == 0) { next }
      res = strsplit(x = labels$V2[index], split = " ", fixed = T)[[1]]
      if (length(res)>1) {
        labels$V2[index] = paste(res[1]," ", res[2],"*" ,sep = "")
      }
    }
  }
  write.table(labels, file =  paste(gene_class, "_itol_labels_signif.txt",sep = ""), col.names = F, row.names = F, quote = F, sep = ",", append = T)
}

### compare the random values to the measured values to see if some subtrees really share more 
## genes than expected given the gene counts
## actually rather than saving to a file I need to load the existing results and I need to add an asterisk 
## where the number of genes I'm seeing is lower or higher than expected by random sampling





