library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(scatterplot3d)

setwd("/Users/gh11/poppunk_pangenome/9_gene_properties/")

classification = read.table("../5_classify_genes/classification_v2.csv", sep = "\t", header = T,
                            stringsAsFactors = F, comment.char = "", quote = "")
classification$gene = gsub("[(/)/']","_" , classification$gene ,ignore.case = TRUE)
colours = read.table("../5_classify_genes/colours_v2.csv", sep = ",", stringsAsFactors = F, comment.char = "", header = T)


props = read.table("easy_gene_props.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T, fill = NA)
harder_props =  read.table("harder_gene_props.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
harder_props$Median = rep(NA, dim(harder_props)[1])
harder_props = harder_props[,c(1,2,3,5,4)]
props = rbind(props,harder_props)


get_value <- function(x) { return(curr_props$Mean[which(curr_props$Gene == x)])}
prop_names = c("length","gravy","aromaticity", "distance","contig_length", "GC")

for (p in prop_names){
  print(p)
  curr_props = props[which(props$Property == p),]
  classification = cbind(classification, as.numeric(curr_props$Mean[match(classification$gene, curr_props$Gene)]))
  colnames(classification)[dim(classification)[2]] = p
}
## save this because it takes quite a long time to aggeragate this data
write.table(classification, "classification_w_props_v2.csv", sep = "\t", col.names = T, row.names = F, quote = F)

classification = read.table("classification_w_props_v2.csv", sep = "\t", header = T,
                            stringsAsFactors = F, comment.char = "", quote = "")

swap_values_with_cols <- function(vec){
  num_cols = length(unique(vec))
  cols = brewer.pal(8, "Set2")
  vec_to_return = as.character(vec)
  i = 1
  for (value in unique(vec)){
    vec_to_return[which(vec == value)] = cols[i]
    i = i + 1
  }
  return(vec_to_return)
}

classification$label = factor(classification$label, paste(1:47, "/47",sep = ""))
ggplot(classification, aes(x = label, y = length, color = gene_class)) + geom_point(alpha = 0.4, size = 3) + 
  geom_boxplot(width = 0.1, colour = "black", fill = NA, outlier.shape = NA)+
  scale_y_continuous(trans='log10')

classification$colour = rep("#d3d3d3",dim(classification)[1]) 
classification$colour[which(classification$total_presence >= 45)] = "green"
classification$colour[which(classification$total_presence <= 3)] = "red"

classification.m = melt(classification[,c(1,18,7,11,13,17)], id.vars = c("gene","colour"))

gc = classification.m[classification.m$variable == "GC",]
ggplot(classification, aes(x = GC, fill = colour)) + geom_density(alpha = 0.4)





plot_prop <- function(p) {
  std_content = aggregate(x = classification[,which(colnames(classification) == p)], by = list(classification$total_presence, classification$gene_class), FUN = sd)
  mean_content = aggregate(x = classification[,which(colnames(classification) == p)], by = list(classification$total_presence, classification$gene_class), FUN = median)
  p1 = ggplot(mean_content, aes(x = Group.1, y = x, fill = Group.2, color = Group.2)) + geom_point(pch = 21, size = 3) + 
    geom_smooth(method='lm', formula= y~x, se = T) + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
    ylab("Mean") + xlab("PopPUNK clusters in which gene was present")+ scale_y_continuous(trans='log10')
  
  p2 = ggplot(std_content, aes(x = Group.1, y = x, fill = Group.2, color = Group.2)) + geom_point(pch = 21, size = 3)+ 
    geom_smooth(method='lm', formula= y~x, se = T) + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
    ylab("Std") + xlab("PopPUNK clusters in which gene was present") +  scale_y_continuous(trans='log10')
  return(list(p1,p2))  
}

grid.arrange(plot_prop("length")[[1]],plot_prop("length")[[2]])





