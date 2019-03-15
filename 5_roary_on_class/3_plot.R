library(ggplot2)
library(RColorBrewer)
library(reshape)

gene_type = "soft_core"
setwd(paste("/Users/gh11/poppunk_pangenome/5_roary_on_class/", gene_type, "/", sep = ""))

### 1. How many genes are shared between how many clusters?
num_members = read.table("membership_count.csv", sep = ",",
                         header = T, stringsAsFactors = F)
num_members$Num_members = factor(num_members$Num_members,1:39 )
ggplot(num_members, aes(x = Num_members, y = Count)) + geom_bar(stat = "identity") +
  theme_classic(base_size = 16) + xlab("Number of clusters") 

## 2. Look at the most common patterns
patterns = read.table("pattern_count.csv", sep = ",",
                      header = T, stringsAsFactors = F)
patterns = patterns[-which(patterns$Num_members == 1),]
## most combinations are only seen once
ggplot(patterns, aes(x = Count)) + geom_histogram(bins = 50)

patterns = patterns[which(patterns$Num_members > 1 & patterns$Count > 10),]

patterns$Pattern = factor(patterns$Pattern, patterns$Pattern[order(patterns$Count)])
ggplot(patterns, aes(x = Pattern, y = Count)) + geom_bar(stat = "identity") +
  coord_flip() + theme_classic(base_size = 16)


## plot the matrix
get_lower_tri<-function(m){
  m[lower.tri(m)] <- NA
  diag(m) <- NA
  return(m)
}

reorder_matrix <- function(m){
  hc <- hclust(dist(m))
  m <-m[hc$order, hc$order]
  return(m)
}


m = as.matrix(read.table("matrix.csv", header = T, row.names = 1, stringsAsFactors = F, sep = ","))
colnames(m) = 1:39
rownames(m) = 1:39
m <- reorder_matrix(m)
lower_tri = get_lower_tri(m)
melted_m <- melt(lower_tri)
melted_m$X1 = factor(melted_m$X1, melted_m$X1[1:length(unique(melted_m$X1))])
melted_m$X2 = factor(melted_m$X2, unique(melted_m$X2))
ggplot(data = melted_m, aes(X2, X1, fill = value)) +
  geom_tile(color = "black") +
  theme_bw(base_size = 16)+
  coord_fixed() + xlab("Cluster") + ylab("Cluster") + theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_fill_gradient(low = "white", high = "blue", space = "Lab", name = "") 

