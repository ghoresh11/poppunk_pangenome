library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(reshape2)
library(scatterplot3d)

setwd("/Users/gh11/poppunk_pangenome/9_gene_properties/")

## What about length?? (that's the most basic one)

gene_classes = read.table("../5_classify_genes/gene_classification.csv", sep = ",", stringsAsFactors = F,
                          comment.char = "", quote = "")
gene_classes$V1 = gsub("[^0-9A-Za-z\\*_-]","_" , gene_classes$V1 ,ignore.case = TRUE)
gene_cats = read.table("../5_classify_genes/descs_template.csv", sep = ",", stringsAsFactors = F,
                       comment.char = "", quote = "", header = T)

harder_props = read.table("harder_gene_props.csv", sep = ",", stringsAsFactors = F,
                          comment.char = "", quote = "", header = T)
easier_props = read.table("easy_gene_props.csv",  sep = ",", stringsAsFactors = F,
                          comment.char = "", quote = "", header = T)
props = rbind(easier_props, harder_props)

props$class = gene_classes$V2[match(props$Gene, gene_classes$V1)]
props = props[-which(is.na(props$class)),]
props$Cat = gene_cats$Main.Label[match(props$class, gene_cats$Var1)]
props$Mean = as.numeric(props$Mean)
props = props[-(22331:22337),]


plot_one_property <- function(curr_prop, ylabel, title, log = F, hline = T, blank = T) {
  curr_df = props[props$Property == curr_prop,]
  if (log) {curr_df$Mean = curr_df$Mean + 2}
  curr_df$Cat = factor(curr_df$Cat, unique(gene_cats$Main.Label))
  curr_df$class = factor(curr_df$class, unique(gene_cats$Var1))
  p = ggplot(curr_df, aes(x = class, y = Mean, fill = class, color = class)) + geom_violin(alpha = 0.8, lwd = 0.2) +
    geom_boxplot(width = 0.2, lwd = 0.2, outlier.size = 0.2)+
    scale_fill_manual(values = gene_cats$Color, guide = F) + scale_color_manual(values = gene_cats$Border, guide = F) +
    theme_classic(base_size = 12)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("") + ylab(ylabel) + scale_x_discrete(labels = rep("",21)) + 
    theme(axis.ticks.x=element_blank())+ggtitle(title)
  if (log) {
    p = p + scale_y_continuous(trans='log10', 
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
      coord_trans( y="log10") }
  p = p +  facet_grid(~Cat, scales = "free",  space = "free",switch = "x") 
  if (hline) {
    p = p + geom_hline(yintercept = mean(props$Mean[props$Property == curr_prop], na.rm = T), lty = 2,  color = "blue", lwd = 0.2)
  }
  if (curr_prop == "GC") {
    p = p + geom_hline(yintercept = 51.6, lty = 2,  color = "blue", lwd = 0.2)
  }
  
  p = p + theme(axis.line.x=element_line(color="black", size = 0.5),
                axis.text.x=element_blank(),
                axis.title.x=element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_blank()
  )
  
  # ggsave(p, filename = paste("/Users/gh11/Submissions/my_thesis/Chapter4/prep/genes_properties/", curr_prop, ".pdf", sep = ""),
  #        width = 11, height = 4)
  
  return(p)
}

plot_property_by_num_clusters <- function(curr_prop, ylabel, title, log = F, hline = T, blank = F) {
  curr_df = props[props$Property == curr_prop,]
  if (log) {curr_df$Mean = curr_df$Mean + 2}
  curr_df$New_class = factor(curr_df$New_class, c("=>10","<10","Secondary"))
  p = ggplot(curr_df, aes(x = New_class, y = Mean)) + geom_violin(alpha = 0.8, lwd = 0.2, fill = "#d3d3d3") +
    geom_boxplot(width = 0.2, lwd = 0.2, outlier.size = 0.2, fill = NA)+
    theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("") + ylab("")  + ggtitle(title)
  if (log) {
    p = p + scale_y_continuous(trans='log10', 
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) +
      coord_trans( y="log10") }
  if (hline) {
    p = p + geom_hline(yintercept = mean(props$Mean[props$Property == curr_prop], na.rm = T), lty = 2, color = "blue", lwd = 0.2)
  }
  if (curr_prop == "GC") {
    p = p + geom_hline(yintercept = 51.6, lty = 2,  color = "blue", lwd = 0.2)
  }
  
  p = p + theme(axis.line.x=element_line(color="black", size = 0.5),
                axis.text.x=element_blank(),
                axis.title.x=element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_blank()) + scale_x_discrete(labels = rep("",21)) + 
    theme(axis.ticks.x=element_blank())
# ggsave(p, filename = paste("/Users/gh11/Submissions/my_thesis/Chapter4/prep/genes_properties/", curr_prop, ".pdf", sep = ""),
#        width = 11, height = 4)

return(p)
}

new_class = rep("Secondary",dim(props)[1])
new_class[which(props$class %in% c( "Ubiquitous (soft-core)", "40-45 varied"  ,   
                                    "Real core", "25-39 varied" ,"10-24 specific",
                                    "40-45 specific"  ,  "Missing in one" , "10-24 varied"  , "25-39 specific" ))] = "=>10"
new_class[which(props$class %in% c("Rare and specific","2-9 varied",   "Intermediate and specific","Multicluster rare",
                                   "Intermediate and rare","2-9 specific","Multicluster intermediate",
                                   "Core and specific"))] = "<10"
props$New_class = new_class



A = plot_one_property("length", "Length (aa)", "A", log = T) 
B = plot_one_property("contig_length", "Contig length (bp)","C", log = T)
C = plot_one_property("distance", "Distance from edge (bp)", "E",log = T)
D = plot_one_property("GC", "%GC", "G", log = F, hline = F)
E = plot_one_property("gravy", "GRAVY score", "I", log = F) 
F_plt = plot_one_property("aromaticity", "fraction aromatic aa", "K", log = F, blank = T) 
G = plot_one_property("start_codon", "%ATG", "M", log = F, hline = F) 

## to run a test:
## pairwise.wilcox.test(g = curr_df$class, x = curr_df$Mean, p.adjust.method = "fdr")

Ax = plot_property_by_num_clusters("length", "Length (aa)", "B", log = T) 
Bx=plot_property_by_num_clusters("contig_length", "Contig length (bp)","D", log = T)
Cx=plot_property_by_num_clusters("distance", "Distance from edge (bp)", "F",log = T)
Dx=plot_property_by_num_clusters("GC", "%GC", "H", log = F, hline = F)
Ex=plot_property_by_num_clusters("gravy", "GRAVY score", "J", log = F) 
Fx=plot_property_by_num_clusters("aromaticity", "fraction aromatic aa", "L", log = F, blank = F) 
Gx=plot_property_by_num_clusters("start_codon", "%ATG", "N", log = F, hline = F) 

layout_matrix = rbind(c(1,1,1,2),
                      c(3,3,3,4),
                      c(5,5,5,6),
                      c(7,7,7,8),
                      c(9,9,9,10),
                      c(11,11,11,12),
                      c(13,13,13,14)
)

grid.arrange(A,Ax,B,Bx,C,Cx,D,Dx,E,Ex,F_plt,Fx,G,Gx, layout_matrix = layout_matrix)




## For the text
mean(props$Mean[which(props$New_class == "=>10" & props$Property == "aromaticity")])
mean(props$Mean[which(props$New_class == "<10" & props$Property == "aromaticity")], na.rm = T)


curr = props[which(props$New_class == "<10" & props$Property == "length"),]

length(which(curr$Mean > 40 & curr$Mean < 60)) / dim(curr)[1]
  length(which(curr$Mean > 150))/ dim(curr)[1]

curr = props[which(props$New_class != "Secondary" & props$Property == "start_codon"),]
length(which(curr$Mean>0.9))/dim(curr)

## see if using a PCA on the matrix of properties, if it seperates the classes well

#$Mean[props$Property %in% c("coting_length","distance","length")] = log10(props$Mean[props$Property %in% c("coting_length","distance","length")] + 2)

props_matrix = dcast(data = props, formula = Gene ~ Property, fill = NA, value.var = "Mean")
#props_matrix = props_matrix[-which(is.na(props_matrix$GC)),]
props_matrix = props_matrix[-which(is.na(props_matrix$gravy)),]
rownames(props_matrix) = props_matrix[,1]
props_matrix = props_matrix[,-1]
props_matrix = props_matrix[,-which( colnames(props_matrix) %in% c("num_members"))] 


#tsne_res = tsne(props_matrix)

pca_res = prcomp(x = props_matrix, center = T, scale. = T)
summary(pca_res)
pca_res = as.data.frame(pca_res$x)
pca_res= cbind(pca_res, Cat = props$Cat[match(rownames(pca_res), props$Gene)])
pca_res= cbind(pca_res, Class = props$class[match(rownames(pca_res), props$Gene)])
pca_res = cbind(pca_res, prevalence = props$New_class[match(rownames(pca_res), props$Gene)])

pca_res$Cat = factor(pca_res$Cat, unique(gene_cats$Main.Label))
pca_res$Class = factor(pca_res$Class, gene_cats$Var1)

#pca_res = pca_res[-which(pca_res$Cat == "Rare"),]

ggplot(pca_res, aes(x = PC3, y = PC2, color = prevalence)) + geom_point() 

ggplot(pca_res, aes(x = PC1, y = PC2, color = Cat)) + geom_point() + 
  scale_color_manual(values = c("Blue","Purple","Green","Red","Gray"))

ggplot(pca_res, aes(x = PC1, y = PC2, color = Class, fill = Class)) + geom_point(pch = 21, size = 2) +
  scale_color_manual(values = gene_cats$Border) + scale_fill_manual(values = gene_cats$Color)

cols = as.character(unlist(pca_res$prevalence))
cols[pca_res$prevalence == "Secondary"] = "yellow"
cols[pca_res$prevalence == "=>10"] = "blue"
cols[pca_res$prevalence == "<10"] = "green"
scatterplot3d(y = pca_res$PC1, z = pca_res$PC2, x = pca_res$PC3, color = cols)
