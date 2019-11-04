library(ggplot2)
library(RColorBrewer)
library(reshape2)


setwd("/Users/gh11/poppunk_pangenome/4.1_correct_pan_genome/")


## see how much the mode and longest differ (currently choosing mode...)
mode_v_longest = read.table("mode_v_longest.csv", sep = ",", comment.char = "", header = T, stringsAsFactors = F, quote = "")
mode_v_longest = cbind(mode_v_longest, diff = mode_v_longest$Longest - mode_v_longest$Mode)
ggplot(mode_v_longest, aes(x = Mode, y = Longest)) + geom_point()

diff_in_length = mode_v_longest[which(mode_v_longest$diff > 50),]
ggplot(diff_in_length, aes(x = diff)) + geom_histogram(binwidth = 10)


### investigate the truncated genes
# 1. How groups do they belong to?
# 2. Are there paritcular clusters that obviosly have one gene truncated in many different ways? (a specific group...)
# 3. Can I understand how it's being truncated (N, C, middle, combination)?

relationships = read.table("trunc_relationships.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T, quote = "")


ggplot(relationships, aes(x = Location, y = Ratio)) + geom_violin() + geom_boxplot(width = 0.5)

## to plot against a particular long gene (doesn't plot ALL the genes in a connected component, only all the genes
## which are connected to a particular other gene)
plot_one_gene<-function(gene){
  check_one = relationships[relationships$GeneA == gene,]
  the_seq = check_one[1,]
  the_seq$Start = 0
  the_seq$Stop = the_seq$LengthA
  the_seq$Location = "Ref"
  the_seq$GeneB = the_seq$GeneA
  the_seq$StdB = the_seq$StdA
  check_one = rbind(the_seq, check_one)
  
  check_one$Location = factor(check_one$Location, c("Ref", "N-Terminus","Middle", "C-Terminus"))
  check_one = check_one[
    order( check_one$Location, -check_one$LengthB ),
    ]
  check_one$GeneB = factor(check_one$GeneB, rev(check_one$GeneB))
  check_one$STD_high = check_one$Stop + check_one$StdB
  check_one$STD_low = check_one$Start - check_one$StdB
  check_one$STD_low[check_one$STD_low < 0] = 0
  p = ggplot() +
    #geom_segment(data = check_one, aes(x = STD_low, y = GeneB, xend = STD_high, yend = GeneB),color = "blue", size = 0.5) +
    geom_segment(data = check_one, aes(x = Start, y = GeneB, xend = Stop, yend = GeneB, col= Location), size = 10) +
    scale_colour_manual(values = rev(brewer.pal(8,"Greys"))) + theme_bw(base_size = 14) +
    ylab("Gene") + xlab("Location (aa)") + ggtitle(gene)
  return(p)
}


only_long = unique(relationships$GeneA[which(!relationships$GeneA %in% relationships$GeneB)])
in_both = intersect(relationships$GeneA, relationships$GeneB) ## for the "in both", I'll consider them truncs as well (just an extra truncation...)
only_short = unique(relationships$GeneB[which(!relationships$GeneB %in% relationships$GeneA)])

long_genes = data.frame(ID = rep("", length(only_long)),
                        Name = only_long,
                        Num_shorter_genes = rep(0, length(only_long)),
                        Num_C_terminus = rep(0, length(only_long)),
                        Num_middle = rep(0, length(only_long)),
                        Num_N_terminus = rep(0, length(only_long)), stringsAsFactors = F)
for (i in 1:length(only_long)){
  name = long_genes$Name[i]
  curr = relationships[which(relationships$GeneA == name),]
  long_genes$ID[i] = curr$ID[1]
  long_genes$Num_shorter_genes[i] = dim(curr)[1]
  long_genes$Num_C_terminus[i] = length(which(curr$Location == "C-Terminus"))
  long_genes$Num_N_terminus[i] = length(which(curr$Location == "N-Terminus"))
  long_genes$Num_middle[i] = length(which(curr$Location == "Middle"))
}

long_genes$N_plus_middle = long_genes$Num_N_terminus + long_genes$Num_middle
ggplot(long_genes, aes(x = Num_C_terminus, y = Num_N_terminus)) + geom_jitter(alpha = 0.3, size = 1) + theme_classic(base_size = 14)
ggplot(long_genes, aes(x = Num_C_terminus, y = N_plus_middle)) + geom_jitter(alpha = 0.3, size = 1) + theme_classic(base_size = 14)

### 
components = read.table("trunc_components.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F, quote = "")
components = components[-which(!components$Gene %in% c(relationships$GeneA, relationships$GeneB)),]
type_per_gene = rep("",dim(components)[1])
for (i in 1:length(type_per_gene)) {
  if (components$Gene[i] %in% only_long){
    type_per_gene[i] = "only long"
  } else if (components$Gene[i] %in% in_both){
    type_per_gene[i] = "both"
  } else {
    type_per_gene[i] = "only short"
  }
}
type_per_gene = data.frame(Gene = components$Gene, type_per_gene, 
                           Num_clusters = components$Num_clusters, stringsAsFactors = F)
ggplot(type_per_gene, aes(x = type_per_gene, y = Num_clusters)) + geom_boxplot() 

## To decide if a gene is a truncated gene or fused, have a look at all the short genes of each long gene
## if the shorter version is found in more clusters, the long is the anomaly 
## if the longer version is found in more cluster, the short is the anomaly
## This is based on this dataset and does not say which is the functional gene, simply which gene is "more different" from the rest
gene_assignments = data.frame(
  Gene = components$Gene,
  Assignment = rep("", dim(components)[1]),
  Partner = rep("", dim(components)[1]),
  Length = components$Length,
  Length_partner =  rep(0, dim(components)[1]),
  Std = components$Std,
  Std_partner = rep(0, dim(components)[1]),
  Num_clusters = components$Num_clusters,
  Num_clusters_partner = rep(0, dim(components)[1]),
  stringsAsFactors = F
)
for (i in 1:dim(gene_assignments)[1]) {
  curr_gene = gene_assignments$Gene[i]
  num_clusters = components$Num_clusters[i]
  
  partner_longer_genes = data.frame(Gene = relationships$GeneA[which(relationships$GeneB == curr_gene)], 
                                    Length = relationships$LengthA[which(relationships$GeneB == curr_gene)],
                                    Std = relationships$StdA[which(relationships$GeneB == curr_gene)],
                                    stringsAsFactors = F)
  partner_longer_genes$Num_cluster = components$Num_clusters[match(partner_longer_genes$Gene, components$Gene)]
  if (dim(partner_longer_genes)[1] > 0) {
    partner_longer_genes = partner_longer_genes[which(partner_longer_genes$Num_cluster == max(partner_longer_genes$Num_cluster)),]
    partner = partner_longer_genes[which(partner_longer_genes$Length == max(partner_longer_genes$Length)),]
    gene_assignments$Partner[i] = partner$Gene[1]
    gene_assignments$Num_clusters_partner[i] = partner$Num_cluster[1]
    gene_assignments$Length_partner[i] = partner$Length[1]
    gene_assignments$Std_partner[i] = partner$Std[1]
    if (num_clusters > partner$Num_cluster[1]){
      gene_assignments$Assignment[i] = "common variant (short)"
    } else if (num_clusters < partner$Num_cluster[1]){
      gene_assignments$Assignment[i] = "secondary variant (short)"
    } else {
      gene_assignments$Assignment[i] = "short variant"
    }
    next
  } 
  partner_shorter_genes = data.frame(Gene = relationships$GeneB[which(relationships$GeneA == curr_gene)], 
                                    Length = relationships$LengthB[which(relationships$GeneA == curr_gene)],
                                    Std = relationships$StdB[which(relationships$GeneA == curr_gene)],
                                    stringsAsFactors = F)
  partner_shorter_genes$Num_cluster = components$Num_clusters[match(partner_shorter_genes$Gene, components$Gene)]
  partner_shorter_genes = partner_shorter_genes[which(partner_shorter_genes$Num_cluster == max(partner_shorter_genes$Num_cluster)),]
  partner = partner_shorter_genes[which(partner_shorter_genes$Length == max(partner_shorter_genes$Length)),]
  gene_assignments$Partner[i] = partner$Gene[1]
  gene_assignments$Num_clusters_partner[i] = partner$Num_cluster[1]
  gene_assignments$Length_partner[i] = partner$Length[1]
  gene_assignments$Std_partner[i] = partner$Std[1]
  if (num_clusters > partner$Num_cluster[1]){
    gene_assignments$Assignment[i] = "common variant (long)"
  } else if (num_clusters < partner$Num_cluster[1]){
    gene_assignments$Assignment[i] = "secondary variant (long)"
  } else {
    gene_assignments$Assignment[i] = "long variant"
  }
}

write.table(gene_assignments, file = "truncated_genes.csv", sep = ",", quote = F, row.names = F, col.names = T)


## median number of short genes per longer gene:
short_per_long = data.frame(table(relationships$GeneA[which(!relationships$GeneA %in% relationships$GeneB)]), stringsAsFactors = F)

## change here to look at any gene, can be quite useful
A = ggplot(short_per_long, aes(x = Freq)) + geom_histogram(binwidth = 0.5) +
  theme_bw(base_size = 14) + xlab("Short variants per long variant") + ylab("Long variants") +
  ggtitle("A") + scale_x_continuous(breaks = 1:15)

B = plot_one_gene("icsA*") + ggtitle("B") + theme(legend.position = "None")
C = plot_one_gene("group_5407**") + ggtitle("C")+ theme(legend.position = "None")
C

relationships$Location = factor(relationships$Location, c("N-Terminus","Middle", "C-Terminus"))
D = ggplot(relationships, aes(x = Location, fill = Location)) + geom_bar() + ylab("Count") + theme_bw(base_size = 14) +
  scale_fill_manual(values = rev(brewer.pal(8,"Greys")[-8]), guide = F) + ggtitle("D") + xlab("Location (on longer variant)")

table(relationships$Location)

grid.arrange(A,B,C,D, layout_matrix = rbind(c(1,2,3),
                                          c(4,2,3)) )



## connect lots of short variants and function
interpro_res = read.table("../5_classify_genes/interproscan_results.csv", sep = "\t", header = T, 
                          stringsAsFactors = F, quote = "", comment.char = "")
short_per_long = cbind(short_per_long, interpro_res[match(short_per_long$Var1, interpro_res$name),])
