library(ggplot2)
library(RColorBrewer)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/4_pairwise_roary/classify_genes/")

## classify genes according to their distribution pattern and see if there
## is a composition of genes which are "universal" to E. coli, or specific to 
## a PopPUNK cluster
gene_classes = read.table("../040919_eps0.5_nodupl/gene_classes.csv", sep = ",", comment.char = "", header = T, 
                          stringsAsFactors = F, quote = "", row.names = 1)
gene_classes = cbind(gene_classes, total = rowSums(gene_classes))

remove_from_gene_classes <- function(vec){
  gene_classes = gene_classes[-which(row.names(gene_classes) %in% vec),]
  return (gene_classes)
}


## writing exactly according to the decision tree:
one = gene_classes[which(gene_classes$total == 1),]
core_and_specific = row.names(one)[which(one$core == 1)]
intermediate_and_specific =  row.names(one)[which(one$inter == 1)]
rare_and_specific = row.names(one)[which(one$rare == 1)]

fortyseven = gene_classes[which(gene_classes$total == 47),]
complete_core = row.names(fortyseven)[which(fortyseven$core == 47)]
ubiq_almost_core = row.names(fortyseven)[which(fortyseven$core != 47)]

other = gene_classes[which(gene_classes$total>1 & gene_classes$total < 47),]

## 0
no_core = other[which(other$core == 0),]
rare_but_multiclade = row.names(no_core)[which(no_core$inter == 0)]
inter_multiclade = row.names(no_core)[which(no_core$rare == 0)]
inter_and_rare = row.names(no_core)[which(no_core$rare != 0 & no_core$inter != 0)]

### Other
get_values_in_range <- function(minimum, maximum) {
  range = other[which(other$core >= minimum & other$core<= maximum),]
  specific = row.names(range)[which(range$core == range$total)]
  varied = row.names(range)[which(range$core != range$total)]
  return(list(specific, varied))
}

rarely = get_values_in_range(1,9)
rarely_specific = rarely[[1]]
rarely_varied = rarely[[2]]

occasional = get_values_in_range(10,24)
occasional_specific = occasional[[1]]
occasional_varied = occasional[[2]]

sometimes = get_values_in_range(25,39)
sometimes_specific = sometimes[[1]]
sometimes_varied = sometimes[[2]]

mostly = get_values_in_range(40,45)
mostly_specific = mostly[[1]]
mostly_varied = mostly[[2]]


missing_in_one = rownames(other)[which(other$core == 46)]


## create classification table for all the genes
all = list(complete_core, ubiq_almost_core,
           missing_in_one, mostly_specific, sometimes_specific, occasional_specific, rarely_specific, core_and_specific,
           intermediate_and_specific, inter_multiclade, rare_and_specific, rare_but_multiclade,
           inter_and_rare, rarely_varied, occasional_varied, sometimes_varied,mostly_varied)

df = read.table("descs_template.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F)
core_cols = rev(brewer.pal(n=4, "Blues")[-1])
core_specific = brewer.pal(n = 6, "Purples")[-1]
inter_cols = brewer.pal(n = 3, "Reds")[-1]
rare_cols = brewer.pal(n = 3, "Oranges")[-1]
varied = brewer.pal(n = 5, "Greens")
df = cbind(df, Color = c(core_cols, core_specific, inter_cols, rare_cols, varied))

df$Cat = factor(df$Cat, rev(c("Core","Core Specific","Intermediate","Rare","Varied")))
df$Name = factor(df$Name, df$Name)
df$Count = as.numeric(df$Count)
df$Color = as.character(df$Color)
A = ggplot(df, aes(x = Cat, fill = Name, y= Count)) + geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = df$Color, name = "") +
  theme_classic(base_size = 12) + ylab("Genes") + coord_flip() + scale_y_continuous(expand = c(0.01,0.01))+
  theme(legend.position = "bottom") + xlab("")+guides(fill=guide_legend(ncol=2,byrow=F))

A

names = as.character(df$Name)
cols = df$Color
file.remove("gene_classification.csv")
for (i in 1:length(all)){
  write.table(file = "gene_classification.csv",
              data.frame(all[[i]], 
                         rep(names[i], length(all[[i]])),
                         rep(cols[i], length(all[[i]]))), col.names = F, row.names = F, quote = F, sep = "\t", append = T)
}

