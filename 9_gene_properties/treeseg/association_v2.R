library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(stringr)

setwd("/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/")

classification = read.table("../classification_w_props_v2.csv", sep = "\t", header = T,
                            stringsAsFactors = F, comment.char = "", quote = "")

association = read.table("associations_final_fix_forR.csv", sep = ",", header = T,
                         stringsAsFactors = F, comment.char = "", quote = "")[,1:4]
association = association[-which(association$gene == ""),]
association$gene = gsub("[(/)/']","_" , association$gene ,ignore.case = TRUE)


association$score = rep(NA, length(association$gene))
association$score[association$combined == "All"] = 1
association$score[association$combined %in% c("treeseg+smallest_subtree","min_subtrees+smallest_subtree")] = 2/3
association$score[association$combined %in% c("treeseg+min_subtrees")] = 1.5/3
association$score[association$combined %in% c("treeseg","smallest_subtree")] = 1/3
association$score[association$combined %in% c("min_subtrees")] = 0.5/3

classification$association = association$score[match(classification$gene, association$gene )]
classification$type = association$combined[match(classification$gene, association$gene )]

test = data.frame(table(classification$fill, classification$type), stringsAsFactors = F)


ggplot(test, aes(x = Var1, fill = Var2, y = Freq)) + geom_bar(stat = "identity")



ggplot(classification, aes(y = means, x = total_presence, fill = association)) + 
  geom_jitter(height = 0, alpha = 0.9, width = 0.1, pch = 21, color = "black", stroke = 0.1, size = 2)  + theme_classic(base_size = 14) +
  ylab("Mean frequency when present") + xlab("PopPUNK clusters in which gene is present") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient(low = "#d3d3d3",high = "black")
