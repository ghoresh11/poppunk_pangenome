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


precentage_explained = data.frame(gene_classes,
                                  num_genes = rep(0, length(gene_classes)),
                                  num_explained = rep(0, length(gene_classes)),
                                  stringsAsFactors = F)
gene_classes = c(gene_classes)

for (gene_class in gene_classes) {
  
  col1 = colours$Colour[which(colours$Class == gene_class)]
  
  # genes_of_interest_freqs = freqs[rownames(freqs) %in% classification$gene[classification$fill == gene_class],]
  
  ## do the same for all of them, remove anything that was present in more than 40
  genes_of_interest_freqs = freqs[rownames(freqs) %in% classification$gene[
    which(classification$fill == gene_class & classification$total_presence < 41)],]
  
  precentage_explained$num_genes[precentage_explained$gene_classes == gene_class] = dim(genes_of_interest_freqs)[1]
  
  fileConn = paste(gene_class, "_itol_pie.txt",sep = "")
  file_open =  file(fileConn)
  writeLines(c("DATASET_PIECHART","SEPARATOR COMMA",
               "DATASET_LABEL,pies",paste("COLOR,",col1,sep = ""),paste("FIELD_COLORS,#1b9e77,#d95f02", sep =""),"FIELD_LABELS,specific,missing","BORDER_WIDTH,2",
               "DATA"), file_open)
  close(file_open)
  
  label_file = paste(gene_class, "_itol_labels.txt",sep = "")
  file_open =  file(label_file)
  writeLines(c("DATASET_TEXT","SEPARATOR COMMA","DATASET_LABEL,labels",paste("COLOR,",col1,sep=""),
               "DATA"), file_open)
  close(file_open)
  complete = c()
  df_out = data.frame(id = paste("INT", 1:length(all_subtrees) + length(all_subtrees) + 1, sep = ""),
                      position = rep(0.3, length(all_subtrees)),
                      size = rep(0, length(all_subtrees)),
                      specific = rep(0, length(all_subtrees)),
                      missing = rep(0, length(all_subtrees)), stringsAsFactors = F)
  
  for (i in length(all_subtrees):1) {
    curr = all_subtrees[[i]]
    df_out$size[i] = length(curr$tip.label)
  }
  df_out = df_out[order(df_out$size, decreasing = T),]
  
  
  for (i in dim(df_out)[1]:2) {
    curr_subtree = as.numeric(gsub(pattern = "INT", x = df_out$id[i], replacement = "")) - length(all_subtrees) - 1
    curr = all_subtrees[[curr_subtree]]
    curr_specifics = find_specific_core_genes(curr$tip.label)
    if (length(complete) > 0){
      curr_specifics = curr_specifics[which(!curr_specifics %in% complete)]
    }
    complete = c(complete, curr_specifics)
    curr_specifics_count = length(curr_specifics)
    df_out$specific[i] = curr_specifics_count
    df_out$size[i] = length(curr$tip.label)
  }
  
  
  for (i in 1:dim(df_out)[1]) {
    curr_subtree = as.numeric(gsub(pattern = "INT", x = df_out$id[i], replacement = "")) - length(all_subtrees) - 1
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
    df_out$missing[i] = curr_missing_count
  }
  df_out$size = df_out$missing + df_out$specific
  precentage_explained$num_explained[precentage_explained$gene_classes == gene_class] = length(complete)
  
  ## write df-out to the file
  write.table(df_out, file = fileConn, col.names = F, row.names = F, quote = F, sep = ",", append = T)
  min_num = max(quantile(probs = 0.7,df_out$specific), quantile(probs = 0.7,c(df_out$missing)))
  for (i in 1:dim(df_out)[1]) {
    if (df_out$specific[i] > min_num | df_out$missing[i]> min_num) {
      if (df_out$specific[i] == 0) {
        write(paste(df_out$id[i],",",df_out$missing[i], ",0.3,black,bold,2,270",sep = ""), file = label_file,append = T)
      } else if (df_out$missing[i] == 0) {
        write(paste(df_out$id[i],",", df_out$specific[i], ",0.3,black,bold,2,270",sep = ""), file = label_file,append = T)
      } else  {
        write(paste(df_out$id[i],",", df_out$specific[i]," ",df_out$missing[i], ",0.3,black,bold,2,270",sep = ""), file = label_file,append = T)
      }
    }
  }
}






## MAKE PIE OF ALL LOW FREQUENCY GENES
## combine the three varied genes into one tree
varied_gene_classes = c("Multi-cluster rare", "Intermediate and rare", "Core, intermediate and rare")

fileConn = "all_mobile_itol_pie.txt"
file_open =  file(fileConn)
writeLines(c("DATASET_PIECHART","SEPARATOR COMMA",
             "DATASET_LABEL,varied_genes","COLOR,#458B00",paste("FIELD_COLORS","#d95f0e","#edf8e9","#bae4b3", sep =","),
             "FIELD_LABELS,Multirare,Inter and rare, Core inter and rare","BORDER_WIDTH,1",
             "DATA"), file_open)
close(file_open)
label_file = "all_mobile_itol_labels.txt"
file_open =  file(label_file)
writeLines(c("DATASET_TEXT","SEPARATOR COMMA","DATASET_LABEL,labels",paste("COLOR,","#458B00",sep=""),
             "DATA"), file_open)
close(file_open)

df_out = data.frame(id = paste("INT", 1:length(all_subtrees) + length(all_subtrees) + 1, sep = ""),
                    position = rep(0.3, length(all_subtrees)),
                    size = rep(0, length(all_subtrees)),
                    multirare = rep(0, length(all_subtrees)),
                    inter_rare = rep(0, length(all_subtrees)), 
                    core_inter_rare = rep(0, length(all_subtrees)), stringsAsFactors = F)

for (j in 1:3) {
  gene_class = varied_gene_classes[j]
  col1 = colours$Colour[which(colours$Class == gene_class)]
  genes_of_interest_freqs = freqs[rownames(freqs) %in% classification$gene[classification$fill == gene_class],]
  
  for (i in length(all_subtrees):1) {
    curr = all_subtrees[[i]]
    df_out$size[i] = length(curr$tip.label)
  }
  df_out = df_out[order(df_out$size, decreasing = T),]
  complete = c()
  for (i in dim(df_out)[1]:2) {
    curr_subtree = as.numeric(gsub(pattern = "INT", x = df_out$id[i], replacement = "")) - length(all_subtrees) - 1
    curr = all_subtrees[[curr_subtree]]
    curr_specifics = find_specific_core_genes(curr$tip.label)
    if (length(complete) > 0){
      curr_specifics = curr_specifics[which(!curr_specifics %in% complete)]
    }
    complete = c(complete, curr_specifics)
    curr_specifics_count = length(curr_specifics)
    df_out[i,j+3] = curr_specifics_count
  }
}
df_out$size = df_out$multirare + df_out$inter_rare + df_out$core_inter_rare


## write df-out to the file
write.table(df_out, file = fileConn, col.names = F, row.names = F, quote = F, sep = ",", append = T)
min_num = quantile(probs = 0.7,df_out$size)
for (i in 1:dim(df_out)[1]) {
  if (df_out$size[i] > min_num) {
    write(paste(df_out$id[i],",", df_out$multirare[i]," ",df_out$inter_rare[i]," ",df_out$core_inter_rare[i], ",0.3,black,bold,2,270",sep = ""), file = label_file,append = T)
  }
}

### plot the percentage explained for the "varied genes"
precentage_explained$percent = precentage_explained$num_explained/precentage_explained$num_genes
precentage_explained = precentage_explained[-which(precentage_explained$percent == 1 | is.nan(precentage_explained$percent)),]
precentage_explained$gene_classes = factor(precentage_explained$gene_classes, colours$Class)
ggplot(precentage_explained, aes(x = gene_classes, y = percent, fill = gene_classes)) + geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = colours$Colour, drop = F, guide = F) + theme_bw(base_size = 14) +
  xlab("") + ylab("Fraction of clade specific genes") + scale_y_continuous(expand = c(0,0,0.1,0))+ 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 


