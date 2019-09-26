library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(plyr)
library(readr)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

all_props_files = list.files(path = "props/", full.names = T)
all_props = ldply(all_props_files, read_csv)


ggplot(all_props, aes(x = type, y = instability)) + geom_boxplot()

ggplot(all_props, aes(x = type, y = GC)) + geom_boxplot()

start = data.frame(table(all_props$type, all_props$start))
ggplot(start, aes(x = Var2, y = Freq, fill = Var1)) + geom_bar(stat = "identity") +scale_fill_gradient(hi)


t = read.csv("/Users/gh11/poppunk_pangenome/6_missing_specific_genes/DBs/ecosyc_All_pathways_of_E._coli_K-12_substr._MG1655_170919.txt",
           sep = "\t", comment.char = "", quote = )
write.table(x = t$Sequence...coordinates.of.DNA.region, file = "/Users/gh11/poppunk_pangenome/6_missing_specific_genes/DBs/gene_locs.txt",
            row.names =  F, col.names = F, quote = F)
