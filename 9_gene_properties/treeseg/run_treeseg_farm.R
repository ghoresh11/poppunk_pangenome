library(ape)
library(treeSeg)

setwd("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/properties/treeseg/")

## LOAD TREE
tree = read.tree("tree_for_treeseg.nwk")

### LOAD THE FREQUENCY TABLE
freqs = read.table("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/correction/corrected_freqs.csv", sep = ",",
                   comment.char = "", stringsAsFactors = F, quote = "", header = F, row.names = 1)
freqs[1,] = as.character(freqs[1,])
colnames(freqs) = freqs[1,]
freqs = freqs[-1,]

## LOAD CLASSIFICTION TABLE
classification = read.table("gene_classification.csv", sep = ",",
                            comment.char = "", stringsAsFactors = F, quote = "", header = F)

run_one_test <- function(gene_name){
  test = freqs[row.names(freqs) == gene_name,]
  test = t(test)[,1][match(tree$tip.label, colnames(test))]
  names(test) = tree$tip.label
  test[which(test>0)] = 1
  test = as.numeric(test)
  seg = treeSeg(as.numeric(unlist(test)), tree, alpha = 0.05)
  return(seg$mlAN)
}

## run tree seg on all
classification$V4 = rep(NA, dim(classification)[1])
for (i in 1:dim(classification)[1]){
  res = run_one_test(classification$V1[i])
  if (length(res)>0){
    classification$V4[i] = paste(res, collapse = "+")
  }
}
write.table(classification,file = "classification_treeseg.csv", sep = ",",quote = F,row.names = F, col.names = F)
