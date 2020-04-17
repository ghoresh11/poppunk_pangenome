library(ggtree)
library(ggplot2)
library(gridExtra)
library(ape)
library(phangorn)

local = T

## How many methods were used in pastml? Using 3 maximum parsimony methods and 3 ML methods
num_methods = 6

setwd("/lustre/scratch118/infgen/team216/gh11/pastml_dir")

if (local) {
  setwd("/Users/gh11/Submissions/my_thesis/Chapter5/prep_v2/corrections/ancestral_state_recon/")
}
### FUNCTIONS

count_gain_loss <- function(res_vector, method, gene_class) {
  counts = list()
  counts[["gain"]] = 0
  counts[["loss"]] = 0
  ## find the node parent
  for (curr_node in names(res_vector)) {
    if (curr_node == "root") {next}
    parent = Ancestors(tree, node = which(o == curr_node), type = "parent")
    parent = o[parent]
    curr_val = res_vector[names(res_vector) == curr_node]
    parent_val = res_vector[names(res_vector) == parent]
    row = which(gain_loss_per_subtree$subtree == curr_node & gain_loss_per_subtree$gene_class == gene_class)
    if (curr_val == "0" && parent_val == "1") { ## loss
      counts[["loss"]] = counts[["loss"]] + 1
      gain_loss_per_subtree[row, method*2 + 1] = 
        gain_loss_per_subtree[row, method*2 + 1] + 1
    } else if (curr_val == "1" && parent_val == "0") { ## gain
      counts[["gain"]] = counts[["gain"]] + 1
      gain_loss_per_subtree[row, (method+1)*2] = 
        gain_loss_per_subtree[row, (method+1)*2] + 1
    }
  }
  return(list(counts,gain_loss_per_subtree))
}




## plot the ancestral state reconstruction outputed from pastml
gene_classification = read.table("classification_v2.csv", sep = "\t",
                                 header = T, stringsAsFactors = F, quote = "")

if (local) {
  gene_classification = read.table("/Users/gh11/poppunk_pangenome/5_classify_genes/classification_v2.csv", sep = "\t",
                                   header = T, stringsAsFactors = F, quote = "")
}

gene_order = read.table("order.txt", header = F, stringsAsFactors = F)[,1]
gene_classification = gene_classification[gene_order,] ## only use the genes that were the input to the ancestral state reconstruction

tree = read.tree("named.tree_tree_for_treeseg.nwk") ## use the output tree

name_corrections = read.table("corrections.csv", sep = ",", header = T, stringsAsFactors = F)

o = c(tree$tip.label, tree$node.label)

## count gain/loss per subtree
gene_classes = unique(gene_classification$fill)
out = c()
for (gene_class in gene_classes) {
  out = c(out, rep(gene_class, length(o)))
}

gain_loss_per_subtree = data.frame(subtree = rep(o,length(gene_classes)),
                                   gene_class = out, 
                                   gains_ACCTRAN = rep(0, length(out)),
                                   losses_ACCTRAN = rep(0, length(out)),
                                   gains_DELTRAN = rep(0, length(out)),
                                   losses_DELTRAN = rep(0, length(out)),
                                   gains_DOWNPASS = rep(0, length(out)),
                                   losses_DOWNPASS = rep(0, length(out)),
                                   gains_JOINT = rep(0, length(out)),
                                   losses_JOINT = rep(0, length(out)),
                                   gains_MAP = rep(0, length(out)),
                                   losses_MAP = rep(0, length(out)),
                                   gains_MPPA = rep(0, length(out)),
                                   losses_MPPA = rep(0, length(out)),stringsAsFactors = F)


gain_loss_per_gene = data.frame(gene = gene_classification$gene,
                                class = gene_classification$fill,
                                gains = rep(0, dim(gene_classification)[1]),
                                losses = rep(0, dim(gene_classification)[1]), stringsAsFactors = F)

## run the analysis on a single gene

for (filecount in seq(from = 1, to = 25049, by = 1000)) {
  print(filecount)
  asr_results = read.table(paste("results/",filecount,".tab.out", sep = ""), sep= "\t", header = T, comment.char = "", stringsAsFactors = F)
  for (j in 1:((dim(asr_results)[2]-1)/num_methods)) {
    curr_index = (j-1)*(num_methods)+2
    gains = c()
    losses = c()
    method = 0
    
    gene_name = colnames(asr_results)[curr_index]
    ## count gains/losses along tree branches
    gene_name = gsub(x = gene_name, pattern = "_ACCTRAN", replacement = "")
    
    
    ## change the name to be as it was for the input of pastml
    gene_name = gsub(x = gene_name, pattern = "XX", replacement = "*", fixed = T)
    gene_name = gsub(x = gene_name, pattern = "ZZ", replacement = "'", fixed = T)
    gene_name = gsub(x = gene_name, pattern = "YY", replacement = "(", fixed = T)
    gene_name = gsub(x = gene_name, pattern = "UU", replacement = ")", fixed = T)
    gene_name = gsub(x = gene_name, pattern = "VV", replacement = "/", fixed = T)
    
    if (gene_name %in% name_corrections$Gene.name) {
      gene_name = name_corrections$Actual.name[name_corrections$Gene.name == gene_name]
    }
    
    gene_index = which(gene_classification$gene == gene_name)
    
    if (length(gene_index)!= 1) { ## for some reason the gene wasn't found 
      print(gene_name)  ## --> this will be saved to the .o file and I will need to check these
      print(gene_index) ## to see if there were no matches or multiple matches
      next
    }
    
    for (i in curr_index:(curr_index+(num_methods-1))) {
      method = method + 1
      curr_gene_results = asr_results[,i]
      names(curr_gene_results) = asr_results$node
      curr_gene_results = curr_gene_results[!is.na(curr_gene_results)]
      
      ## find cases where the ancestral state is ambigous
      both = names(which(table(names(curr_gene_results))>1))
      curr_gene_results[names(curr_gene_results) %in% both] = "ambiguous"
      for (name in both) {
        if (length(which(names(curr_gene_results) == name)) > 1) {
          curr_gene_results = curr_gene_results[-which(names(curr_gene_results) == name)[-1]] ## remove duplicate entries
        }
      }
      
      curr_counts = count_gain_loss(curr_gene_results, method, gene_classification$fill[gene_index])
      gain_loss_per_subtree = curr_counts[[2]]
      gains = c(gains, curr_counts[[1]][["gain"]])
      losses = c(losses, curr_counts[[1]][["loss"]])
      ## Plotting
      ## order the results to match the tip order on the tree
      if (local) {
        curr_gene_results = curr_gene_results[match(o,names(curr_gene_results))]
        curr_gene_results = factor(curr_gene_results, c("0","1","ambiguous"))
        print(ggtree(tree) + geom_label(aes(label = c(tree$tip.label,tree$node.label), fill = curr_gene_results)) +
                ggtitle(method) + scale_fill_manual(values = c("#d2d2d2","#78AB46","#f8d591"), drop = F))
      }
    }
    methods_used = c(1,2,4,5) # acctran, deltran, join, map
    gain_loss_per_gene$gains[gene_index] = median(gains[methods_used])
    gain_loss_per_gene$losses[gene_index] = median(losses[methods_used])
  }
}
write.table(gain_loss_per_gene, file = "gain_loss_per_gene.csv", col.names = T, row.names = F, quote = F, sep= "\t")
write.table(gain_loss_per_subtree, file = "gain_loss_per_subtree.csv", col.names = T, row.names = F, quote = F, sep= "\t")
