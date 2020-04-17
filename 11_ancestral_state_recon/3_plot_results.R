library(ggplot2)
library(reshape2)
library(ggimage)
library(ape)
library(ggtree)
library(gridExtra)

setwd("/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/ancestral_state_recon/")

classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep = "\t",
                            header = T, comment.char = "", stringsAsFactors = F, quote = "")
colours = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/colours_v2.csv", header = T, comment.char = "",
                     stringsAsFactors = F, sep = ",")


## Figure for gain/loss per gene

### Calculations looking at the median/mean gain/loss events across gene classes
gain_loss_per_gene = read.table(file = "gain_loss_per_gene.csv", sep = "\t", header = T, quote = "", stringsAsFactors = F, comment.char = "")

gain_loss_per_gene = read.table(file = "all_genes.tab", sep = "\t", header = T, quote = "", stringsAsFactors = F, comment.char = "")
gain_loss_per_gene$class = classification$fill[match(gain_loss_per_gene$gene, classification$gene)]

gain_loss_per_gene$total_presence = classification$total_presence[match(gain_loss_per_gene$gene, classification$gene)]

gain_loss_per_gene.m = melt(gain_loss_per_gene[,-c(1,5)], id.vars = "class") 
o=   c("Multi-cluster core",
       "Core and intermediate","Multi-cluster intermediate","Core and rare",
       "Core, intermediate and rare", "Intermediate and rare",
       "Multi-cluster rare")
gain_loss_per_gene.m$class = factor(gain_loss_per_gene.m$class,
                                  o)
plots = list()
plots[["A"]] = ggplot(gain_loss_per_gene.m, aes(x = class, fill = variable, y = value )) + geom_boxplot() +
  xlab("Occurence class") + ylab("Gain/loss events per gene") + theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Dark2", name = "Event type", labels = c("Gain","Loss"))+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom") + 
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2)) + ggtitle("A") 


#pairwise.wilcox.test(gain_loss_per_gene$gains,gain_loss_per_gene$class, method = "fdr")
gene_class_order = c("Multi-cluster core", 
                     "Core and intermediate", 
                     "Multi-cluster intermediate",
                     "Core and rare",
                     "Core, intermediate and rare", 
                     "Intermediate and rare",
                     "Multi-cluster rare")
letters = c("B","C","D","E","F","G","H")
i = 0
for (gene_class in gene_class_order) {
  i = i + 1
  test = gain_loss_per_gene[gain_loss_per_gene$class == gene_class,]
  test$gains[test$gains >= 5] = ">4"
  test$losses[test$losses >= 5] = ">4"
  curr_summary = data.frame(table(test$gains, test$losses), stringsAsFactors = F)
  colnames(curr_summary) = c("gains","losses","count")
  curr_summary$count = curr_summary$count / dim(test)[1]
  for (gain in c(0:4, ">4")) {
    for (loss in c(0:4, ">4")) {
      index = which(curr_summary$gains == gain & curr_summary$losses == loss)
      if (length(index) == 0) {
        curr_summary = rbind(curr_summary,
                             data.frame( 
                               gains = gain,
                               losses = loss,
                               count = 0, stringsAsFactors = F))
      }
    }
  }
  curr_summary$gains = factor(curr_summary$gains, c(0:4, ">4"))
  curr_summary$losses = factor(curr_summary$losses, c(0:4, ">4"))
  
  plots[[letters[i]]] = ggplot(curr_summary, aes(x = gains, y = losses, fill = count)) + geom_tile(color = "black") +
    theme_classic(base_size = 14) + 
    scale_fill_gradient(low = "white", high = "black", name = "Fraction of genes", breaks = seq(from = 0, to = 0.5, by = 0.1), limits = c(0,0.5)) +
    xlab("Gain events") + ylab("Loss events") +
    guides(ticks = F) + labs(title = letters[i], subtitle = gene_class) +
    scale_y_discrete(expand = c(0,0)) + theme(axis.line.x = element_blank(), axis.line.y = element_blank()) +
    scale_x_discrete(na.value = "black")
}
legend = as_ggplot(get_legend(plots[[i]]))
for (letter in letters){
  plots[[letter]] = plots[[letter]] + theme(legend.position = "None")
}

## I saved is as 975x770
grid.arrange(grobs = plots, layout_matrix = rbind(c(1,2,3),
                                                  c(1,4,5),
                                                  c(6,7,8)))
legend
## for the text, it is useful to have the median for each occurrence class
medians = aggregate(gain_loss_per_gene.m$value, by = list(gain_loss_per_gene.m$class, gain_loss_per_gene.m$variable), FUN = median)
median_total_presence = aggregate(gain_loss_per_gene$total_presence, by = list(gain_loss_per_gene$class), FUN = median)
medians$total_presence = median_total_presence$x[match(medians$Group.1, median_total_presence$Group.1)]

### FIgure for gain/loss pair subtree

## Look at each subtree to count how many gains/loss events on each internal node
tree = read.tree("/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/tree_for_treeseg.nwk")
clades = read.table("phylogroup_clades.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
all_subtrees = read.table("gain_loss_per_subtree.csv", sep = "\t", comment.char = "", stringsAsFactors = F, header = T)
all_subtrees = read.table("all_subtrees.tab", sep = "\t", comment.char = "", stringsAsFactors = F, header = T)



## make a new plot where the x axis is the size of the node (from root to tip), the y axis is the 
## relative number of gain/loss events happening on the branch and that should be coloured by occurrence class
all_subtrees.m = melt(all_subtrees[,c(1,2,6,3,4,5)], id.vars = c("id","name","size", "class"))

all_subtrees.m$label = clades$Name[match(all_subtrees.m$id, clades$Node)]
all_subtrees.m$label[all_subtrees.m$label %in% c("C","U")] = "Other"
all_subtrees.m$label[is.na(all_subtrees.m$label)] = "Other"
all_subtrees.m$label[all_subtrees.m$size == 1] = "Tip"

all_subtrees.m$class = factor(all_subtrees.m$class, o)
all_subtrees.m$label = factor(all_subtrees.m$label, c(clades$Name,"Tip","Other"))

aggregated = aggregate(all_subtrees.m$value, by = list(all_subtrees.m$class, all_subtrees.m$variable, all_subtrees.m$label), FUN = sum)
colnames(aggregated) = c("class", "event","clade","value")
aggregated$value[aggregated$clade == "Tip"] = aggregated$value[aggregated$clade == "Tip"]/47
## how many internal nodes are there that aren't the phylogroups
num_other = 46 - (length(unique(aggregated$clade))-2)
aggregated$value[aggregated$clade == "Other"] = aggregated$value[aggregated$clade == "Other"]/num_other

clade_sizes = data.frame(clade = unique(aggregated$clade), size = rep(0, length(unique(aggregated$clade))), stringsAsFactors = F)
clade_sizes$size = all_subtrees.m$size[match(clade_sizes$clade, all_subtrees.m$label)]
clade_sizes$size[clade_sizes$clade == "Other"] = 0
clade_sizes = clade_sizes[order(clade_sizes$size, decreasing = T),]


aggregated$clade = factor(aggregated$clade, clade_sizes$clade)
xlabs = paste(clade_sizes$clade, "\n(n=",clade_sizes$size, ")", sep="")
xlabs[length(xlabs)] = "Mean across other branches\n(n=2-47)"
xlabs[length(xlabs)-1] = "Mean across all tips\n(n=1)"
aggregated$class = factor(aggregated$class, o)
C= ggplot(aggregated, aes(x = clade, y = value, fill = event)) + geom_bar(width = 0.5,color = "black",stat = "identity", position = "dodge") +
  facet_grid(class~., scales = "free", labeller = label_wrap_gen(width=7)) +
  theme_bw(base_size = 14) + ylab("Gain/loss events") + xlab("Clade") +
  scale_y_continuous(expand = c(0,0,0.1,0)) + 
  scale_fill_brewer(palette = "Dark2", labels = c("Gain","Loss"), name = "Event type") +
  scale_x_discrete(labels = xlabs)+ 
  theme( panel.grid.minor = element_blank()) + ggtitle("C")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


## Build the basic tree plot
gene_class = "Multi-cluster core"

curr = all_subtrees[all_subtrees$class %in% gene_class,]
curr$label = clades$Name[match(curr$id, clades$Node)]
curr$fill = clades$Colour[match(curr$id, clades$Node)]

A =  ggtree(tree, aes(color=curr$gains), size = 1.2) + geom_tiplab(size = 4) + geom_tippoint()+geom_treescale(x = -0.0001, y = 44, width = 0.05)+
    scale_color_viridis_c(name = "Gain events", direction = -1) + theme(legend.position="bottom")+
    geom_label(aes(x=branch, label=curr$label, fill = curr$label), color = "black") +
    scale_fill_manual(values = clades$Colour, guide = F) + ggtitle("A")


B=  ggtree(tree, aes(color=curr$losses), size = 1.2) + geom_tiplab(size = 4) + geom_tippoint()+geom_treescale(x = -0.0001, y = 44, width = 0.05)+
    scale_color_viridis_c(name = "Loss events", direction = -1) + theme(legend.position="bottom")+
    geom_label(aes(x=branch, label=curr$label, fill = curr$label), color = "black") +
    scale_fill_manual(values = clades$Colour, guide = F) + ggtitle("B")



grid.arrange(A, B, C, layout_matrix=rbind(c(1,2,3,3)))



## Plot the number of events per tip-> might be useful for me to talk about but too much to present
 ### tips plot -> sack it
tips = all_subtrees.m[all_subtrees.m$size == 1,]
tips$normalised_value = rep(NA, dim(tips)[1])

normalise <- function(df, type, gene_class) {
  filtered_df = df[df$variable == type & df$class == gene_class, ]
  normalised_value = (filtered_df$value - min(filtered_df$value))/(max(filtered_df$value) - min(filtered_df$value))
  df$normalised_value[df$variable == type & df$class == gene_class] = normalised_value
  return(df)
}
for (gene_class in o) {
  for (type in c("gains", "losses")) {
    tips = normalise(tips, type, gene_class)
  }
}

cluster_sizes = read.table("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",",
                              comment.char = "", header = T, stringsAsFactors = F)
tips$name = factor(tips$name, cluster_sizes$Cluster)
tips$class = factor(tips$class, o)
ggplot(tips, aes(x = name, y = value, label = name)) + geom_bar(stat = "identity") +
  facet_grid(class~variable, scales = "free") +
  theme_bw(base_size = 14) + ylab("Gain/loss events") + xlab("PopPUNK Cluster")+
  scale_y_continuous(expand = c(0,0,0.1,0))  + geom_text(color = "red") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = rep("",47))


## to plot barplots on internal nodes, really not easy to use, gave up on it
# bars = list()
# for (id in curr$id) {
#   df.m = melt(curr[curr$id == id, c(1,4,5)], id.vars = "id")
#   df.m$label = df.m$value
#   df.m$label[df.m$label<med_val] = ""
#   bars[[id]] = ggplot(df.m, aes(x = variable, y = value, fill = variable, label = label)) + geom_bar(stat = "identity") + 
#     scale_fill_brewer(palette = "Dark2", guide = F) +
#     scale_y_continuous(limits = c(0, max_val), expand = c(0,0,0.1,0)) + 
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(), 
#           axis.title.y=element_blank(),
#           axis.text.y=element_blank(),
#           axis.ticks.y=element_blank(), 
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#           panel.border = element_blank(), 
#           panel.background = element_rect(fill = "transparent",colour = NA),
#           plot.background = element_rect(fill = "transparent",colour = NA))  + geom_text(nudge_y  = -1)
#     
# }
# names(bars) = curr$id               
# 
# inset(tree_view = p,insets =  bars, height = 4, width = 0.018, vjust = 0.5)  ## if you can't see anything, the height + width are essential)
# 
# ## multi-cluster core was saved as 1460x668
# 
