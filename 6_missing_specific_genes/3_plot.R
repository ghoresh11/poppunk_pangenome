library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)
library(ape)
library(ggtree)
library(dplyr)

setwd("/Users/gh11/poppunk_pangenome/6_missing_specific_genes/")


#### MAIN ####

phylogroup_order =  c("B2", "F", "D","E","A","C","B1","U") ## define
type_order = c("wrong", "truncated","longer","same functional group","true")
cols = c("black", "#CC6677","#DDCC77","#88CCEE","#44AA99")

## 1. Using analysis on merged genes to remove wrong ones
specific = read.csv("specific_genes_detailed.csv", sep ="\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
specific_details = read.table("specific_gene_types.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F, quote = "")
#specific_details$Type[which(specific_details$Type == "wrong")] = "truncated"
specific_details = specific_details[match(specific$gene, specific_details$Gene),]
specific = cbind(specific, specific_details)


specific$cluster = as.character(specific$cluster)
specific$Type = factor(specific$Type, type_order)
specific$phylogroup = factor(specific$phylogroup,phylogroup_order)
ps1 = ggplot(specific, aes(x = cluster, fill = Type)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols)
ps = ggplot(specific, aes(x = cluster)) + geom_bar(color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ps1
ps

## plot it proportionally? -> doesn't tell me much
specific_prop = data.frame(prop.table(table(specific$cluster, specific$Type),1))
specific_prop = cbind(specific_prop, phylogroup = specific$phylogroup[match(specific_prop$Var1, specific$cluster)])
ggplot(specific_prop, aes(x = Var1, fill = Var2, y = Freq)) + geom_bar(stat = "identity", color = "black") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = cols)
counts = data.frame(table(specific$cluster)) ## for text, counts


missing = read.csv("missing_genes_detailed.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
missing_details = read.table("missing_gene_types.csv", sep = ",", header = T, comment.char = "", stringsAsFactors = F, quote = "")
missing_details = missing_details[match(missing$gene, missing_details$Gene),]
missing = cbind(missing, missing_details)

missing$cluster = as.character(missing$cluster)
missing$Type = factor(missing$Type, type_order)
missing$phylogroup = factor(missing$phylogroup,phylogroup_order)


for (cluster in unique(specific$cluster)) { ## add clusters with no missing genes
  if (!cluster %in% missing$cluster){
    add = missing[1,]
    add$cluster[1] = cluster
    add$num_genes[1] = 0
    add$phylogroup[1] = specific$phylogroup[specific$cluster == cluster][1]
    for (i in 4:length(add)) {
      add[1,i] = NA
    }
    missing = rbind(missing, add)
  }
}
pm1 = ggplot(missing, aes(x = cluster, fill = Type, color = Type)) + geom_bar() + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(cols,NA),
                    labels = c("wrong","truncated","longer","same functional group","truly missing in cluster","")) + 
  scale_color_manual(values = c(rep("black", 6),NA), guide = F) 

missing_counts = data.frame(cluster = unique(missing$cluster), counts = rep(0,47), 
                            phylogroup = rep("",47),stringsAsFactors = F)
for (i in 1:dim(missing_counts)[1]){
  missing_counts$counts[i] = missing$num_genes[which(missing$cluster == missing_counts$cluster[i])][1]
  missing_counts$phylogroup[i] = missing$phylogroup[which(missing$cluster == missing_counts$cluster[i])][1]
}

pm = ggplot(missing_counts, aes(x = cluster, y = counts)) + geom_bar(color = "black", stat = "identity") + theme_classic(base_size = 12)  + 
  ylab("Genes") + xlab("PopPUNK cluster") + scale_y_continuous(expand = c(0.01,0,0.1,0)) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
pm
pm1

ps = ps + ggtitle("Specific genes")
pm = pm + ggtitle("Missing genes")
grid.arrange(ps, pm, ncol = 1)+ theme(legend.position = "None")

legend = as_ggplot(get_legend(ps1))
pm1 = pm1 + ggtitle("Missing genes") + theme(legend.position = "None")
ps1 = ps1 + ggtitle("Specific genes")+ theme(legend.position = "None")

grid.arrange(ps1, pm1, ncol = 1)



## should have the linear regression between missing and specific
grid.arrange(pm1, ps1, legend, layout_matrix = rbind(c(1,1,NA),
                                                     c(2,2,3)))

## look at connection between number of fused genes/truncated genes and cluster size
## more genomes = more fused genes, which means it's a probability thing rather
## than an actual thing, alot of the fused genes are called in correctly
size_num_genes = data.frame(table(specific$Type, specific$cluster), stringsAsFactors = F)
colnames(size_num_genes) = c("type","cluster","count")
cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",",
                           header = T, stringsAsFactors = F)
size_num_genes = cbind(size_num_genes, size = cluster_sizes$Size[match(size_num_genes$cluster, cluster_sizes$Cluster)])
for (kind in c("longer","truncated","wrong")) {
  curr = size_num_genes[which(size_num_genes$type == kind),]
  print(ggplot(curr, aes(x = size, y = count)) + geom_point(alpha = 0.7, size = 2) + ggtitle(kind) + 
          geom_smooth(method='lm', col = "black")+ theme_bw(base_size = 14) + xlab("Cluster size") + ylab("Wrong gene clusters"))
}


## create a summary for each lineage after the correction where I have
## how many truncated genes (based on partner), how many fused genes?, how many truly specific?, how many truly missing?
## continue here

## add the partnered genes to the missing genes 
all_partnered_genes = as.character(specific$Partner)
all_partnered_genes = all_partnered_genes[-which(is.na(all_partnered_genes))]
all_partnered_genes = as.character(unlist(sapply(FUN = strsplit, X = all_partnered_genes, split = "/", fixed = T)))
for (gene in all_partnered_genes){
  index = which(missing$gene == gene)
  if (length(index) == 0){next}
  missing$Partner[index] = gene
}


### Add information of pathotype and drug resistance to see if it's connected??
pathytope_mdr = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_sts.csv",
                           sep = ",", comment.char = "", stringsAsFactors = F, header = T)

num_clusters = length(unique(specific$cluster))
final_summary = data.frame(clusters = unique(specific$cluster),
                           phylogroup = rep("", num_clusters),
                           pathotype = rep("", num_clusters),
                           mdr = rep("", num_clusters),
                           truncated = rep(0, num_clusters),
                           fused = rep(0, num_clusters),
                           truly_specific = rep(0, num_clusters),
                           truly_missing = rep(0, num_clusters),
                           stringsAsFactors = F)

for (i in 1:dim(final_summary)[1]){
  cluster = final_summary$clusters[i]
  final_summary$phylogroup[i] = specific$phylogroup[which(specific$cluster == cluster)][1]
  final_summary$pathotype[i] = pathytope_mdr$Pathotype[which(pathytope_mdr$cluster == cluster)][1]
  final_summary$mdr[i] = pathytope_mdr$MDR[which(pathytope_mdr$cluster == cluster)][1]
  curr = specific[which(specific$cluster == cluster),]
  final_summary$truncated[i] = length(unique(curr$Partner[which(curr$Type == "truncated")]))
  final_summary$fused[i] = length(curr$Partner[curr$Type == "longer"])
  final_summary$truly_specific[i] = length(which(curr$Type %in% c("true",  "same functional group")))
  final_summary$truly_missing[i] = length(which(missing$cluster == cluster & is.na(missing$Partner)))
}
final_summary.m = melt(final_summary, id.vars = c("clusters","phylogroup", "pathotype","mdr"))
final_summary.m$clusters = as.character(final_summary.m$clusters)
final_summary.m$variable = factor(final_summary.m$variable, c("truly_specific","truly_missing","truncated","fused"))
final_summary.m$phylogroup = factor(final_summary.m$phylogroup, phylogroup_order)

final_missing.m = final_summary.m[which(final_summary.m$variable != "truly_specific"),]
final_missing = ggplot(final_missing.m, aes(x = clusters, y = value, fill = variable)) + geom_bar(stat = "identity", color = "black") + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") +
  scale_fill_manual(values = rev(brewer.pal(3, "Reds")), labels = c("Truly missing","Truncated","Fused"), name = "") + 
  theme_classic(base_size = 14) + xlab("PopPUNK Cluster") + ylab("Genes") + ggtitle("Missing")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
legend = as_ggplot(get_legend(final_missing))
final_missing = final_missing + theme(legend.position = "None")

final_specific.m = final_summary.m[which(final_summary.m$variable == "truly_specific"),]
final_specific = ggplot(final_specific.m, aes(x = clusters, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black", fill = brewer.pal(3,"Blues")[3]) + 
  facet_grid(~phylogroup, scales = "free",  space = "free",switch = "x") + 
  theme_classic(base_size = 14) + xlab("PopPUNK Cluster") + ylab("Genes") + ggtitle("Specific")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
final_specific

## doesn't look like there's any particular connection between pathotype, being MDR to the number of cluster specific or
## missing genes
grid.arrange(final_missing, final_specific, nrow = 2)
legend


## find a story to tell
## cluster 10 is a non-O157 EHEC which has a lot of specific and missing genes, it's not really drug resistant
curr_cluster = 45
cluster10 = specific[which(specific$cluster == curr_cluster),]
#cluster10 = specific
cluster10 = data.frame(apply(cluster10, 2, gsub, patt=",", replace="."), stringsAsFactors = F)
truly_missing = missing[which(missing$cluster == curr_cluster & is.na(missing$Partner)),]
#truly_missing = missing[which(is.na(missing$Partner)),]
truly_missing  = data.frame(apply(truly_missing, 2, gsub, patt=",", replace="."), stringsAsFactors = F)
truly_missing$Type = rep("truly missing", length(truly_missing$Type))

cluster10 = rbind(cluster10, truly_missing)
cluster10 = cluster10[, c(1,14:16,4:13)]
cluster10 = cluster10[order(cluster10$Type, cluster10$Partner),]
#cluster10 = cluster10[-which(cluster10$Type == "wrong"),]
# write.table(cluster10, file = paste("stories/cluster","_all",".csv", sep = ""),
#             sep = ",", col.names = T, row.names = F, quote = F)
write.table(cluster10, file = paste("stories/cluster",curr_cluster,".csv", sep = ""),
            sep = ",", col.names = T, row.names = F, quote = F)


### Lengths of truly specific genes -> I suspect a lot of them are short (maybe add info on function, how many unknown??)
specific_lengths = read.table("stories/all_specific_gene_lengths.csv", sep = ",",
                              header = T, stringsAsFactors = F, comment.char = "")
mean(specific_lengths$Length)
A = ggplot(specific_lengths, aes(x = Length)) + geom_histogram(binwidth = 20) +
  scale_x_continuous(limits = c(0,2000)) + theme_bw(base_size = 14) + 
  geom_vline(xintercept = 100, color = "red") + xlab("Length(aa)")
num_not = length(which(specific_lengths$Functionally_annotated == "No"))
num_yes = length(which(specific_lengths$Functionally_annotated == "Yes"))

## 30% of specific genes have no functional annotation and these tend to be much shorter...
## they're likely not real or truncated versions of other genes, so they're not really cluster specific
## but are truncated genes
B = ggplot(specific_lengths, aes(y = Length, x = Functionally_annotated)) + geom_violin(fill = "#d3d3d3") +
  geom_boxplot(width = 0.1)+ theme_classic(base_size = 14) +
  scale_x_discrete(labels = 
                     c(paste("No\n(n=",num_not, ")",sep = ""),
                       paste("Yes\n(n=",num_yes, ")",sep = ""))) + ylab("Length(aa)")

short_genes_compared_to_trunc = data.frame(cluster = unique(specific_lengths$Cluster),
                                           num_trunc = rep(0, num_clusters),
                                           num_under_100 = rep(0, num_clusters), stringsAsFactors = F)
for (i in 1:dim(short_genes_compared_to_trunc)[1]){
  curr_cluster = short_genes_compared_to_trunc$cluster[i]
  short_genes_compared_to_trunc$num_trunc[i] = specific_lengths$Num_trunc_genes[which(specific_lengths$Cluster == curr_cluster)][1]
  short_genes_compared_to_trunc$num_under_100[i] = length(which(specific_lengths$Cluster == curr_cluster & specific_lengths$Length <= 100))
}
C =  ggplot(short_genes_compared_to_trunc, aes(x = num_trunc, y = num_under_100, label = cluster)) + geom_text() +
  xlab("Number of truncated genes") + ylab("Number of specific genes<100aa") +
  theme_bw(base_size = 14) + 
  geom_smooth(method='lm',formula=y~x, colour = "black")

grid.arrange(A,B,C, ncol = 3)



### Gene context -> when a gene is a truncated or longer variant of another gene
## it's important to check it if it's in the same genetic location, I did that by looking at the genes around it
curr_context = read.table("stories/contexts/prn.csv", sep = ",", comment.char = "",
                        stringsAsFactors = F, header = T)






