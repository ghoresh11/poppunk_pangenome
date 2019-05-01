library(ggplot2)
library(RColorBrewer)

setwd("/Users/gh11/poppunk_pangenome/7_tree/")

### metadata_file for tree
cluster_colors = c(
  brewer.pal(n = 8, "Blues")[-1],
  brewer.pal(n = 8, "Reds")[-1],
  brewer.pal(n = 8, "Greens")[-1],
  brewer.pal(n = 8, "Greys")[-1],
  brewer.pal(n = 8, "Purples")[-1],
  brewer.pal(n = 8, "PuOr")[1:4]
)
legend = data.frame(cluster = 1:39, col = cluster_colors, stringsAsFactors = F)
legend$cluster = factor(legend$cluster, 1:39)
ggplot(legend, aes(x = cluster,y = cluster, color = cluster)) + geom_point(size = 5) +
  theme_classic(base_size = 16) + scale_color_manual(values = legend$col)


tree = read.table("chosen_for_tree.csv", sep = ",", header = T, comment.char = "",
                  quote = "", stringsAsFactors = F)
tree = data.frame(name = tree$Annotation, shape = rep(2,  length(tree$popppunk_cluster)),
                  size = rep(10,  length(tree$popppunk_cluster)),
                  cluster = tree$popppunk_cluster,
                  col = rep("", length(tree$popppunk_cluster)),
                  fill = rep(3,  length(tree$popppunk_cluster)),
                  loc = rep(-1,  length(tree$popppunk_cluster)),
                  font = rep("bold",  length(tree$popppunk_cluster)),
                  stringsAsFactors = F)
for (i in 1:dim(tree)[1]) {
  tree$col[i] = cluster_colors[tree$cluster[i]]
  name = strsplit(x = tree$name[i], split = "/", fixed = T)[[1]]
  name = gsub(pattern = ".velvet.gff", replacement = "",x = name[length(name)], fixed = T)
  name = gsub(pattern = ".gff", replacement = "",x = name, fixed = T)
  tree$name[i] = name
}

tree = data.frame(tree$name, tree$cluster, tree$loc, tree$col, tree$font, tree$fill, rep(0, dim(tree)[1]),stringsAsFactors = F)
## save to file
write.table(tree, file = "tree_labels.txt", sep = ",", col.names = F, row.names = F, quote = F)
