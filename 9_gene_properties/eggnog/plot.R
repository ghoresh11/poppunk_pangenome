library(ggplot2)
library(dplyr)
library(stringr)

setwd("/Users/gh11/poppunk_pangenome/9_gene_properties/eggnog/")

cogs = read.csv("eggnog_results.emapper.annotations", sep = "\t",stringsAsFactors = F, header = F, comment.char = "#")

cogs_desc = read.csv("cog_descs.csv", sep = ",", comment.char = "", stringsAsFactors = F, quote = "", header = T)
cogs_desc$Title = str_wrap(cogs_desc$Title, width = 10)

class_desc = read.csv(file = "../../5_classify_genes/descs_template.csv", sep = ",", comment.char = "",
                      stringsAsFactors = F, header = T)
class_desc$Main.Label[which(class_desc$Main.Label == "Intermediate")] = "Inter."

classification = read.table("../../5_classify_genes/gene_classification.csv", sep = ",", comment.char = "",
                            stringsAsFactors = F, header = F, quote = "")

classification$COG = cogs$V12[match(classification$V1, cogs$V1)]
classification$COG[is.na(classification$COG) | classification$COG == ""] = "?"
classification = cbind(classification, cogs_desc[match(classification$COG, cogs_desc$COG),])
classification = classification[,-4]


## by fraction of genes from each class rather than pure counts
counts  = table(classification$V2, classification$COG)
counts = counts/ rowSums(counts)
counts = data.frame(counts, stringsAsFactors = F)
counts$Title = cogs_desc$Title[match(counts$Var2, cogs_desc$COG)]
counts = cbind(counts, class_desc[match(counts$Var1, class_desc$Var1),])
counts = counts[,-(5:7)]

counts$Var1 = factor(counts$Var1, class_desc$Var1)
counts$Main.Label = factor(counts$Main.Label, unique(class_desc$Main.Label))
counts$Var2 = factor(counts$Var2, cogs_desc$COG)
counts$Title = factor(counts$Title, unique(cogs_desc$Title))

ggplot(counts, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity", color = "black", lwd = 0.2) + 
  facet_grid(Title~Main.Label, scales = "free_x",  space = "free_x",switch = "both") +
  scale_fill_manual(values = cogs_desc$Col, name = "") + theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("Fraction of genes")



### to plot the words in particular gene glasses for some COG cats

word_files = dir("words/", full.names = T) #where you have your files
words = do.call(rbind,lapply(word_files,read.csv, sep = ",", comment.char = "", quote = "", stringsAsFactors = F, header = T))
words$Title =  cogs_desc$Title[match(words$COG, cogs_desc$COG)]
words$Class_cat = class_desc$Main.Label[match(words$Class, class_desc$Var1)]


## counts per class

make_one_cog <- function(cog){
  test = words[words$COG == cog,]
  
  curr_count_per_word = aggregate(x = test$Count, by = list(test$Label, test$Class_cat), FUN = sum)
  curr_count_per_word = curr_count_per_word[order(curr_count_per_word$Group.2,curr_count_per_word$x, decreasing = T),]
  
  words_to_keep = c()
  for (gene_class in unique(curr_count_per_word$Group.2)) {
    words_to_add = as.character(curr_count_per_word$Group.1[which(curr_count_per_word$Group.2 == gene_class)])
    if (length(words_to_add) > 5) { words_to_add = words_to_add[1:5]}
    words_to_keep = c(words_to_keep, words_to_add)
  }
  words_to_keep = unique(words_to_keep)
  
  count_per_class = aggregate(test$Count, list(test$Class_cat), FUN = sum)
  test$fraction = test$Count/ count_per_class$x[match(test$Class_cat, count_per_class$Group.1)]
  
  test = test[which(test$Label %in% words_to_keep),]
  
  test = rbind(test,
               data.frame(
                 Class = class_desc$Var1, 
                 COG = rep("XX", dim(class_desc)[1]),
                 Word = rep(test$Word[1], dim(class_desc)[1]),
                 Label = rep(test$Label[1], dim(class_desc)[1]),
                 Count = rep(0, dim(class_desc)[1]),
                 Title = rep("XX", dim(class_desc)[1]),
                 Class_cat = class_desc$Main.Label,
                 fraction = rep(0, dim(class_desc)[1])
               ))
  
  test$Label = factor(test$Label, words_to_keep)
  test = test[order(test$Label),]
  
  test$Class = factor(test$Class, class_desc$Var1)
  test$Class_cat = factor(test$Class_cat, unique(class_desc$Main.Label))
  test$COG = factor(test$COG, cogs_desc$COG[which(cogs_desc$COG %in% test$COG)])
  
  test = test[-which(test$Class_cat == "Secondary"),]
  
  p = ggplot(test, aes(x = Label, y = fraction, fill = Class, color = Class)) + geom_bar(stat = "identity", lwd = 0.2) + 
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = class_desc$Color, guide = F) +
    scale_color_manual(values = class_desc$Border, guide = F) +
    facet_grid(Class_cat~., scales = "free_x",  space = "free",switch = "both") 
 ggsave(plot = p, filename = paste("plots/", cog, ".pdf", sep = ""), height = 10, width = 15)
  means = aggregate(test$fraction, by = list(test$Class_cat, test$Label), FUN = mean)
  return(means)
}

###
gene_types = class_desc$Var1
for (curr_gene_type in gene_types){
  word_counts = read.table(paste("words/", gsub(x = curr_gene_type, pattern = " ", replacement = "_"), "_words.csv", sep = ""),
                           sep = ",", comment.char = "", quote = "", stringsAsFactors = F, header = T)
  word_counts$Title =  cogs_desc$Title[match(word_counts$COG_Cat, cogs_desc$COG)]
  word_counts$Color = cogs_desc$Col[match(word_counts$COG_Cat, cogs_desc$COG)]
  word_counts = word_counts[-which(word_counts$Count < 4), ]
  p = ggplot(word_counts, aes(x = Word, y = Count, fill = COG_Cat))+ geom_bar(stat = "identity", color = "black", lwd = 0.2) + 
    facet_grid(.~Title, scales = "free_x",  space = "free_x",switch = "both") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(labels = word_counts$Label) + ggtitle(curr_gene_type)
  print(p)
}


K = make_one_cog("K")
L = make_one_cog("L")
M = make_one_cog("M")
T_cog = make_one_cog("T")

## for text
summary =  aggregate(counts$Freq, by = list(counts$Var1, counts$Title), FUN = sum)
summary_2 = aggregate(counts$Freq, by = list(counts$Var1, counts$Var2), FUN = sum)
summary_2$Group.1 = factor(summary_2$Group.1, class_desc$Var1)
summary_2$title = class_desc$Main.Label[match(summary_2$Group.1, class_desc$Var1)]
summary_2$title = factor(summary_2$title, unique(class_desc$Main.Label))
summary_2$Group.2 = factor(summary_2$Group.2 , unique(cogs_desc$COG))

chosen = c("J","K","L","N","M","V","U","O","T")
summary_2 = summary_2[which(summary_2$Group.2 %in% chosen),]
ggplot(summary_2, aes(x = Group.1, y = x,fill = Group.2, group = Group.2))  +
  geom_point(size = 2, pch = 21, color = "black")+
  scale_fill_manual(values = cogs_desc$Col[cogs_desc$COG %in% chosen], guide = F)+
  facet_grid(Group.2~title, scales = "free",  space = "free_x",switch = "both") + theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0,0.01,0,0.03)) +
  xlab("") + ylab("Fraction of genes")



