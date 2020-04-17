library(ape)
library(treeSeg)
library(ggplot2)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/9_gene_properties/treeseg/")

### LOAD THE TREE
tree = read.tree("/Users/gh11/poppunk_pangenome/7_AMR_vir_plasmid/smaller_tree/raxml_tree_mod.nwk")
tree = root(tree,outgroup = "NC_011740")
tree = drop.tip(tree, tip =  "NC_011740")
for (i in c(21,43,49)){
  tree = drop.tip(tree, tip =  as.character(i))
}
## need to fix the rooting... I think I should save the tree then reload it!!
write.tree(tree, file = "tree_for_treeseg.nwk")
tree = read.tree("tree_for_treeseg.nwk")
tree_labels = read.table("subtree_labels.csv", sep = ",", header=T, comment.char = "", stringsAsFactors = F)
all_subtrees = subtrees(tree = tree)

### LOAD THE FREQUENCY TABLE
freqs = read.table("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv", sep = ",",
                   comment.char = "", stringsAsFactors = F, quote = "", header = F, row.names = 1)
freqs[1,] = as.character(freqs[1,])
colnames(freqs) = freqs[1,]
freqs = freqs[-1,]

### LOAD THE CLASSIFICATION TABLE
classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/gene_classification.csv", sep = ",",
                            comment.char = "", stringsAsFactors = F, quote = "", header = F)
all_classes = unique(classification$V2)

class_desc = read.table("../../5_classify_genes/descs_template.csv", sep = ",", comment.char = "",
                        stringsAsFactors = F, header = T)

### FUNCTIONS ####

## Run a single treeSeg test and plot it ##
run_one_test <- function(gene_name, plot = T){
  test = freqs[row.names(freqs) == gene_name,]
  test = t(test)[,1][match(tree$tip.label, colnames(test))]
  names(test) = tree$tip.label
  test[which(test>0)] = 1
  test = as.numeric(test)
  seg = treeSeg(as.numeric(unlist(test)), tree, alpha = 0.05)
  if (plot) { plot_seg_result(seg, list(test)) }
  return(list(seg, test))
}

run_multiple_tests <- function(genes) {
  tests = list()
  segs = list()
  for (g in genes) {
    res = run_one_test(g, plot = F)
    segs[[g]] = res[[1]]
    tests[[g]] = res[[2]]
  }
  plot_seg_result(segs, tests)
  
}

## Plot results of treeSeg ##
plot_seg_result <- function(seg, tests){
  ### plot
  n = length(tree$tip.label)
  lwdEdge <- rep(1.5, dim(tree$edge)[1])
  if (length(tests) == 1) {
    layout(mat=matrix(1:3,ncol = 1),heights = rep(c(2,0.1, 1),6))
    test = tests[[1]]
  } else {
    layout(mat=matrix(1:(2+length(tests)),ncol = 1),heights = c(2,rep(0.1,length(tests))))
    seg = seg[[1]]
  }
  par(mar = c(0.1,1,0.1,0.1))
  
  plot(tree, type = "phylogram", show.tip.label = T , use.edge.length = F, cex = 1.5,
       node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 
  
  # plot maximum likelihood estimate for the nodes where the distribution of phenotype changes.
  ## to add the red dot on the clade...
  if (length(tests) == 1 && length(seg$mlAN) > 0 ) {
    for(i in 1:length(seg$mlAN)){
      nodelabels("", seg$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
    }
  }
  # plot tip phenos
  for (test in tests) {
    plot(1:n, rep(0,n), axes = F, bg = c("white", "black")[test+1],
         col = c("black", "black")[test+1], pch = 22, cex = 4,xlab = '',ylab = '') 
  }
  if (length(tests) > 1 ) { return() }
  # plot the maximum likelihood estimates for the phenotype distributions for each segment.
  par(mar = c(1,2,0.1,0.1))
  plot(1:n, seg$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')
  
  # plot the confidence interval for the phenotype distributions.
  polygon(c(1:n, rev(1:n)), c(seg$confBandP[,2], rev(seg$confBandP[,1])),
          col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)
  
  #axis(1, at = c(1,50,100,150,200), cex = 2)
  axis(2, at = c(0,0.5,1),labels = c('0','0.5','1') )
  mtext(side=2,'Prob.',line=2.5,cex=1)
}



## method to find the smallest subtree (whether for random or real data)
find_smallest_subtree <- function(gene_name, real = T){
  if (real) {
    test = freqs[row.names(freqs) == gene_name,]
    test = t(test)[,1][match(tree$tip.label, colnames(test))]
    test = as.numeric(test)
    names(test) = tree$tip.label
    members = names(test)[which(test>0)]
  } else {
    test = rep(0, 47)
    test[sample(1:47,size = gene_name,replace = F)] = 1 ## assume uniform dist on all tips
    names(test) = tree$tip.label
    members = names(test)[which(test>0)]
  }
  subtree_counts = data.frame(gene = rep(gene_name, length(all_subtrees)),
                              subtree = (1+47):(47 + length(all_subtrees)),
                              size = rep(0, length(all_subtrees)),
                              gene_cluster_size = rep(0, length(all_subtrees)))
  for (i in 1:length(all_subtrees)) {
    curr_tips = all_subtrees[[i]]$tip.label
    subtree_counts$gene_cluster_size[i] = length(which(members %in% curr_tips))
    subtree_counts$size[i] = length(curr_tips)
    
  }
  subtree_counts = subtree_counts[which(subtree_counts$gene_cluster_size == length(members)),]
  subtree_counts = subtree_counts[which(subtree_counts$size == min(subtree_counts$size)),]
  subtree_counts$score = subtree_counts$gene_cluster_size/subtree_counts$size ## the larger the better
  subtree_counts$misfits = subtree_counts$size - subtree_counts$gene_cluster_size 
  return (subtree_counts)
}

## Method to find the minimum number of subtrees that contain leaves with the gene with no errors ##
count_minimal_subtrees <- function(gene_name, real = T) {
  if (real) {
    test = freqs[row.names(freqs) == gene_name,]
    test = t(test)[,1][match(tree$tip.label, colnames(test))]
    test = as.numeric(test)
    names(test) = tree$tip.label
    members = names(test)[which(test>0)]
    num_members = length(members)
  } else {
    num_members = gene_name
    test = rep(0, 47)
    test[sample(1:47,size = gene_name,replace = F)] = 1 ## assume uniform dist on all tips
    names(test) = tree$tip.label
    members = names(test)[which(test>0)]
  }
  
  subtrees = c()
  while(length(members) > 1){
    subtree_counts = data.frame(gene = rep(gene_name, length(all_subtrees)),
                                subtree = 48:(length(all_subtrees)+47),
                                size = rep(0, length(all_subtrees)),
                                num_correct = rep(0, length(all_subtrees)),
                                num_incorrect = rep(0, length(all_subtrees)), stringsAsFactors = F)
    for (i in 1:length(all_subtrees)) {
      curr_tips = all_subtrees[[i]]$tip.label
      subtree_counts$num_correct[i] = length(which(members %in% curr_tips))
      subtree_counts$num_incorrect[i] = length(which(!curr_tips %in% members))
      subtree_counts$size[i] = length(curr_tips)
    }
    subtree_counts = subtree_counts[which(subtree_counts$num_incorrect == 0),]
    if (dim(subtree_counts)[1] == 0) {break }
    subtree_counts = subtree_counts[which(subtree_counts$num_correct == max(subtree_counts$num_correct)),]
    chosen = subtree_counts$subtree[1]
    members = members[-which(members %in% all_subtrees[[chosen-47]]$tip.label)]
    subtrees = c(subtrees, chosen)
  }
  subtrees = c(subtrees, members)
  for (i in 1:length(subtrees)){
    if (subtrees[i] %in% tree_labels$Node) {
      subtrees[i] = tree_labels$Label[which(tree_labels$Node == subtrees[i])]
    }
  }
  num_subtrees = length(subtrees)
  num_singletons = length(members)
  res = data.frame(gene_name = gene_name, 
                   gene_cluster_size = num_members,
                   num_subtrees = num_subtrees,
                   num_singletons = num_singletons,
                   score = num_members / num_subtrees,  ## the larger the better
                   subtrees = paste(subtrees, collapse = ","), stringsAsFactors = F)
  return (res)
}

## Method to get random distributions of the smallest subtree for all presence values ##
generate_random_values_smallest_subtree <- function(n){
  random_samples = data.frame(sample = 1:n)
  for (i in 1:47) {
    minimal_subtree = rbind.fill(lapply(X = rep(i,n), FUN = find_smallest_subtree, real = F))
    random_samples = cbind(random_samples, minimal_subtree$score)
  }
  random_samples = random_samples[,-1]
  colnames(random_samples) = 1:47
  write.table(random_samples, "smallest_subtree_random_samples.csv", quote = F, col.names = T, row.names = F, sep = ",")
  return (random_samples)
}

## Method to get random distributions of the minimum number of subtrees for all presence values ##
generate_random_values_min_num_subtrees <- function(n) {
  random_samples = data.frame(sample = 1:n)
  for (i in 1:47) {
    min_num_subtrees = rbind.fill(lapply(X = rep(i,n), FUN = count_minimal_subtrees, real = F))
    random_samples = cbind(random_samples, min_num_subtrees$score)
  }
  random_samples = random_samples[,-1]
  colnames(random_samples) = 1:47
  write.table(random_samples, "min_num_trees_random_samples.csv", quote = F, col.names = T, row.names = F, sep = ",")
  return(random_samples)
}


## check if a number if signifactnly associated with a lineage/multiple lineages ##
check_significance <- function(df, random_samples, level = 0.95) {
  signif_vals = rep(0, 47)
  for (i in 1:47){
    signif_vals[i] = quantile(random_samples[,i], level)
  } ## anything higher than that value is significant, not seen by random (or seen randomly 5% of times)
  df$signif = rep("No",dim(df)[1])
  for (i in 1:dim(df)[1]){
    curr_val = signif_vals[ as.numeric(df$gene_cluster_size[i])]
    if (round(curr_val, 5) == round(df$score[i],5)) {
      next
    }
    if (df$score[i] > curr_val) {
      df$signif[i] = "Yes"
    }
  }
  return (df)
}


## Method to plot all possible lineage associations for one gene class ##
generate_plot_for_one_class <- function(curr_class){
  curr_genes = classification$V1[which(classification$V2 == curr_class)]
  
  minimal_subtree = rbind.fill(lapply(X = curr_genes, FUN = find_smallest_subtree))
  minimal_subtree$misfits = minimal_subtree$size - minimal_subtree$gene_cluster_size
  minimal_subtree$score = minimal_subtree$gene_cluster_size/minimal_subtree$size
  
  minimal_subtree$gene_cluster_size = as.character(minimal_subtree$gene_cluster_size)
  
  minimal_subtree$subtree= minimal_subtree$subtree + 47
  xlab_order =  unique(minimal_subtree$subtree[order(minimal_subtree$size, decreasing = T)])
  minimal_subtree$subtree = factor(minimal_subtree$subtree,xlab_order)
  
  x_labels = tree_labels$Label[match(xlab_order, tree_labels$Node)]
  minimal_subtree$label = tree_labels$Label[match(minimal_subtree$subtree, tree_labels$Node)]
  
  A = ggplot(minimal_subtree, aes(x = subtree)) + geom_bar(color = "black") + 
    theme_bw(base_size = 14)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_discrete(labels = x_labels)
  B = ggtree(tree) + geom_tiplab() + geom_text2(aes(subset=!isTip, label=c(1:47,tree_labels$Node)), hjust=-.1, color = "blue")
  return(list(minimal_subtree,A,B))
}


#### MAIN ####

RANDOMISE = F
CALCULATE = F

## Step 1: generate random samples for both type of scores
if (RANDOMISE) {
  smallest_subtree_random_samples = generate_random_values_smallest_subtree(n = 1000)
  min_subtrees_random_samples = generate_random_values_min_num_subtrees(n = 1000)
} else {
  smallest_subtree_random_samples = read.table("smallest_subtree_random_samples.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
  min_subtrees_random_samples = read.table("min_num_trees_random_samples.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
}

if (CALCULATE) {
  ## Step 2: generate real numbers for all genes
  smallest_subtree_genes = rbind.fill(lapply(X = classification$V1, FUN = find_smallest_subtree))
  min_subtrees_genes = rbind.fill(lapply(X = classification$V1, FUN = count_minimal_subtrees))
} else {
  smallest_subtree_genes = read.table("smallest_subtree.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
  min_subtrees_genes = read.table("all_minimal_subtrees.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
}

## Step 3: Modify dataframes
smallest_subtree_genes$class = classification$V2[match(smallest_subtree_genes$gene, classification$V1)] 
smallest_subtree_genes$cat = class_desc$Main.Label[match(smallest_subtree_genes$class, class_desc$Var1)]
min_subtrees_genes$class = classification$V2[match(min_subtrees_genes$gene, classification$V1)] 
min_subtrees_genes$cat = class_desc$Main.Label[match(min_subtrees_genes$class, class_desc$Var1)]

## reorder so they match
min_subtrees_genes = min_subtrees_genes[match(smallest_subtree_genes$gene, min_subtrees_genes$gene),]

## Step 4: write to file
write.table(smallest_subtree_genes, "smallest_subtree.csv", sep = ",", col.names = T, row.names = F, quote = T)
write.table(min_subtrees_genes, "all_minimal_subtrees.csv", sep = ",", quote = T, col.names = T, row.names = F)

## Step 5: Check significance of association
smallest_subtree_genes = check_significance(smallest_subtree_genes, smallest_subtree_random_samples, level = 0.95)
min_subtrees_genes = check_significance(min_subtrees_genes, min_subtrees_random_samples, level = 0.99)

## Step 6: Read treeSeg results
tree_seg_results = read.table("classification_treeseg.csv", sep = ",", stringsAsFactors = F, header = F, quote = "", comment.char = "")
tree_seg_results = tree_seg_results[match(min_subtrees_genes$gene, tree_seg_results$V1),]
tree_seg_results$signif = rep("Yes", dim(tree_seg_results)[1])
tree_seg_results$signif[which(is.na(tree_seg_results$V4))] = "No"


## Step 7 summarise results -> three scores, which are significant?
all_results = data.frame(gene = min_subtrees_genes$gene,
                         class = min_subtrees_genes$class,
                         size =  min_subtrees_genes$gene_cluster_size,
                         treeseg = tree_seg_results$signif,
                         smallest_subtree = smallest_subtree_genes$signif,
                         min_subtrees = min_subtrees_genes$signif,
                         treeseg_node = tree_seg_results$V4,
                         smallest_subtree_node = smallest_subtree_genes$subtree,
                         stringsAsFactors = F)

all_results$combined = rep("No", dim(all_results)[1])
all_results$combined[all_results$smallest_subtree == "Yes" & all_results$treeseg == "Yes" & all_results$min_subtrees == "Yes"] = "All"
all_results$combined[all_results$smallest_subtree != "Yes" & all_results$treeseg == "Yes" & all_results$min_subtrees == "Yes"] = "treeseg+min_subtrees"
all_results$combined[all_results$smallest_subtree == "Yes" & all_results$treeseg != "Yes" & all_results$min_subtrees == "Yes"] = "min_subtrees+smallest_subtree"
all_results$combined[all_results$smallest_subtree == "Yes" & all_results$treeseg == "Yes" & all_results$min_subtrees != "Yes"] = "treeseg+smallest_subtree"
all_results$combined[all_results$smallest_subtree != "Yes" & all_results$treeseg != "Yes" & all_results$min_subtrees == "Yes"] = "min_subtrees"
all_results$combined[all_results$smallest_subtree == "Yes" & all_results$treeseg != "Yes" & all_results$min_subtrees != "Yes"] = "smallest_subtree"
all_results$combined[all_results$smallest_subtree != "Yes" & all_results$treeseg == "Yes" & all_results$min_subtrees != "Yes"] = "treeseg"

all_results$subtree_lineage = tree_labels$Label[match(all_results$smallest_subtree_node, tree_labels$Node)]
all_results$treeseg_lineage = tree_labels$Label[match(all_results$treeseg_node, tree_labels$Node)]


all_genes_summary = table(all_results$class, all_results$combined)
all_genes_summary = all_genes_summary/ rowSums(all_genes_summary)
all_genes_summary = data.frame(all_genes_summary, stringsAsFactors = F)
all_genes_summary$Cat = class_desc$Main.Label[match(all_genes_summary$Var1, class_desc$Var1)] 

## factorise
all_genes_summary$Cat = factor(all_genes_summary$Cat, unique(class_desc$Main.Label))
all_genes_summary$Var1 = factor(all_genes_summary$Var1, class_desc$Var1)

all_genes_summary = all_genes_summary[-which(all_genes_summary$Var2 == "No"),]
all_genes_summary = all_genes_summary[-which(all_genes_summary$Cat == "Core"),]
all_genes_summary$Var2 = factor(all_genes_summary$Var2, c("treeseg", "min_subtrees", "smallest_subtree",
                                                          "treeseg+min_subtrees","min_subtrees+smallest_subtree","treeseg+smallest_subtree",
                                                          "All"))
## remove trivials
trivials = c("Rare and specific","Intermediate and specific","Core and specific")
all_genes_summary = all_genes_summary[-which(all_genes_summary$Var1 %in% trivials),]

## counts of genes associated with a lineage (there's no correction for multiple testing, so can assume that about 5% is wrong...)
ggplot(all_genes_summary, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity", color = "black") +
  facet_grid(.~Cat, scales = "free_x",  space = "free",switch = "x") + 
  scale_fill_brewer(palette = "YlGnBu", name = "", labels = 
                      c("TreeSeg", "minimum trees","smallest subtree", "TreeSeg + minimum trees", "minimum trees + smallest subtree",
                        "TreeSeg + smallest subtree", "All")) +
  theme_classic(base_size = 14) + xlab("") + ylab("Fraction of genes")+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + 
  theme(legend.position = "bottom")


## for the text
sum(all_genes_summary$Freq[which(all_genes_summary$Var1  == "25-39 specific")])



## plot a tree with the genes that have different associations
all = "dctB"
only_treeseg = "higA_3*"
only_min_tree = "bhsA_4"
only_smallest_tree = "pikAV_3"
smallest_min_subtrees = "group_3567"
smallest_treeseg = "group_4943"
treeseg_min_subtrees = "tsaR_2"
none = "soj"
run_multiple_tests(c(all, smallest_treeseg, smallest_min_subtrees, treeseg_min_subtrees,
                     only_smallest_tree, only_treeseg, only_min_tree, none))




## the associations mean different things (i don't trust only treeSeg or only min_trees)
tree_seg_min_subtree = all_results[all_results$combined %in% c("treeseg+min_subtrees"),] ## indicates more lineages/deletion in a lineage
all = all_results[all_results$combined %in% c("All"),] ## strong association
smallest_subtree = all_results[all_results$combined %in% c("smallest_subtree"),]  ## confinement to a clade

plot_one_group <- function(df, column) {
  colnames(df)[which(colnames(df) == column)] = "lineage"
  curr_class_desc = class_desc[which(class_desc$Var1 %in% df$class),]
  df$class = factor(df$class, curr_class_desc$Var1)
  df$lineage = factor(df$lineage, names(sort(table(df$lineage), decreasing = F)))
  p = ggplot(df, aes(x = lineage,fill = class, color = class)) + geom_bar(lwd = 0.3)+  
    theme_classic(base_size = 14)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = curr_class_desc$Color) +
    scale_color_manual(values = curr_class_desc$Border) +
    xlab("Clade") + ylab("Genes") + theme(legend.position = "None") + coord_flip()
  return (p)
}

signif = all_results[all_results$combined != "No",]
signif = signif[signif$combined != "min_subtrees",]
signif$lineage = signif$subtree_lineage 
signif$lineage[which(!is.na(signif$treeseg_lineage))] = signif$treeseg_lineage[which(!is.na(signif$treeseg_lineage))]
signif$Cat = class_desc$Main.Label[match(signif$class, class_desc$Var1)]
curr_class_desc = class_desc[which(class_desc$Var1 %in% signif$class),]
signif$class = factor(signif$class, curr_class_desc$Var1)
signif$lineage = factor(signif$lineage, names(sort(table(signif$lineage), decreasing = T)))
signif$Cat = factor(signif$Cat, unique(curr_class_desc$Main.Label))


ggplot(signif, aes(x = lineage, fill = class, color = class))  + geom_bar(lwd = 0.3)+  
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(Cat~., scales = "free_y",space = "free_x",switch = "x") +
  scale_fill_manual(values = curr_class_desc$Color) +
  scale_color_manual(values = curr_class_desc$Border) + xlab("Clade") +
  ylab("Genes")


## if plotting per significance type
A = plot_one_group(all, "subtree_lineage") + ggtitle("C")
tree_seg_min_subtree = tree_seg_min_subtree[-which(is.na(tree_seg_min_subtree$treeseg_lineage)),]
B = plot_one_group(tree_seg_min_subtree, "treeseg_lineage") + ggtitle("D")
C = plot_one_group(smallest_subtree, "subtree_lineage") + ggtitle("E")

grid.arrange(A,B,C, nrow = 1)

## write an output file that can be used in downstream analysis to distinguish genes which are associated 
## with a clade compared to those that aren't
results_output = all_results[,c(1,2,3,9,10,11)]
write.table(results_output, file = "associations_final.csv", sep = "\t", col.names = T, row.names = F, quote = F)





