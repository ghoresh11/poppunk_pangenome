library(ggplot2)
library(RColorBrewer)


setwd("/Users/gh11/e_colis/poppunk/rare_genes/")


df = read.table("rare_filtering_summary.csv", sep = ",",
                header = T, stringsAsFactors = F, comment.char = "")
df$Cluster = factor(df$Cluster, 1:max(unique(df$Cluster)))


ggplot(df, aes(x = Threshold, y = Decrease, col = Cluster)) + geom_line()
ggplot(df, aes(x = Threshold, y = Gene_count, col = Cluster)) + geom_line()


df$Threshold = as.character(df$Threshold)
df$Threshold = factor(df$Threshold, c(50,60,70,80,90,100))
ggplot(df, aes(x = Cluster, y = Gene_count, fill = Threshold)) + geom_bar(stat = "identity", position = "dodge")

ggplot(df, aes(x = Cluster, y = Decrease, fill = Threshold)) + geom_bar(stat = "identity", position = "dodge")

df$Threshold = as.character(df$Threshold)
ggplot(df, aes(x = Threshold, y = Gene_count)) +
  geom_violin(trim=FALSE, fill="gray")  +
  geom_boxplot(width=0.1)+
  theme_classic(base_size = 16)

cluster_sizes = read.table("../dists_analysis/cluster_sizes.csv", sep = ",",
                           header = T, comment.char = "", stringsAsFactors = F)
df = cbind(df, Size = cluster_sizes$Size[match(df$Cluster, cluster_sizes$Cluster)])

fifty = df[which(df$Threshold == 50),]
ggplot(fifty, aes(x = Size, y = Decrease)) + geom_point()

eighty = df[which(df$Threshold == 80),]
ggplot(eighty, aes(x = Size, y = Decrease)) + geom_point()
