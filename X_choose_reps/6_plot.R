library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)


setwd("/Users/gh11/poppunk_pangenome/X_choose_reps/")


pairwise_dists = read.table("pairwise_dists.csv", sep = ",", header = T, stringsAsFactors = F)

type = rep("between", dim(pairwise_dists)[1])
type[pairwise_dists$Cluster1 == pairwise_dists$Cluster2] = "within"

pairwise_dists = cbind(pairwise_dists, type)


ggplot(pairwise_dists, aes(x = type, y = Min)) + geom_violin(fill = "#d2d2d2") + geom_boxplot(width = 0.1) +
  theme_bw(base_size = 16) + xlab("") + ylab("Mean core distance identity") 


pairwise_dists$Cluster1 = factor(pairwise_dists$Cluster1, 1:51)
pairwise_dists$Cluster2 = factor(pairwise_dists$Cluster2, 1:51)
ggplot(pairwise_dists, aes(x = Cluster1, y = Cluster2, fill = Mean)) + geom_tile(color = "black")+
  scale_fill_gradient(low = "white", high = "black") + theme_classic(base_size = 12)

ggplot(pairwise_dists, aes(x = Min, y = Mean)) + geom_point(size = 3, alpha = 0.7) +
  theme_bw(base_size = 16)

for (curr_cluster in 1:51){
  if (curr_cluster == 50) {next}
  curr_dists = pairwise_dists[which(pairwise_dists$Cluster1 == curr_cluster | pairwise_dists$Cluster2 == curr_cluster),]
  curr_dists = cbind(curr_dists, max = curr_dists$Mean + curr_dists$SD, min = curr_dists$Mean - curr_dists$SD)
  curr_dists = curr_dists[order(curr_dists$Mean, decreasing = T),]
  curr_dists = cbind(curr_dists, label = paste(curr_dists$Cluster1,"_",curr_dists$Cluster2, sep = ""))
  curr_dists$label = factor(curr_dists$label, curr_dists$label)
  p = ggplot(curr_dists, aes(x = label, y = Mean)) + geom_point() + 
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,color = "black", position=position_dodge(0.05))  +
    theme_bw(base_size = 13) + scale_y_continuous(limits = c(0.85, 1))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
}
