library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(ggpubr)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

complete_presence_absence = fread("../4_pairwise_roary/231019_corrected//complete_presence_absence.csv", sep = ",", header = T, stringsAsFactors = F)
classification = read.table("gene_classification.csv", sep = ",", header = F, stringsAsFactors = F, comment.char = "", quote = "")

strains = colnames(complete_presence_absence)[-1]
genes = complete_presence_absence[,1][-1]


clusters = unlist(complete_presence_absence[1,])[-1]
random_samples = c()
for (curr in unique(clusters)) {
  print(curr)
  samples = which(clusters == curr)
  samples = sample(x = samples ,size = 15, replace = T)
  random_samples = c(random_samples, samples)
}

# res =  data.frame( cluster = numeric(0), desc = character(0), count = numeric(0), stringsAsFactors = F)
# for (i in random_samples){
#   i = i+1
#   cluster = unlist(complete_presence_absence[1,..i])
#   vec = complete_presence_absence[,..i][-1]
#   vec = unlist(genes[which(vec>0)])
#   for (class in unique(classification$V2)) {
#     curr_count = length(intersect(classification$V1[which(classification$V2 == class)], vec))
#     res = rbind(res, data.frame( cluster = cluster, desc = class, count = curr_count, stringsAsFactors = F))
#   }
# }
# write.table(res, file = "class_per_isolate.csv", sep = ",", row.names = F, col.names = T, quote = F) ## randomly using only 15 per cluster
res = read.table(file = "class_per_isolate.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "")
summarised = aggregate(.~cluster+desc, res, mean)

## add the phylogroup information from the graphics file
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "")

summarised = cbind(summarised, phylogroup= graphics$Phylogroup[match(summarised$cluster, graphics$Cluster)])
summarised$cluster = as.character(unique(summarised$cluster))
all_order = c("Real core", "Ubiquitous (soft-core)","Missing in one","40-45 specific",
              "25-39 specific","10-24 specific","2-9 specific","Core and specific",
              "Multicluster intermediate", "Intermediate and specific",
              "Multicluster rare","Rare and specific",
              "40-45 varied","25-39 varied", "10-24 varied","2-9 varied",
              "Intermediate and rare","Secondary (specific)","Secondary (intermediate)",
              "Secondary (rare)","Secondary (varied)")

summarised$desc = factor(summarised$desc, rev(all_order))
cols = read.table("descs_template.csv", sep = ",", comment.char = "",
                  stringsAsFactors = F, header = T)


B = ggplot(summarised, aes(x = cluster, y = count, fill = desc)) + geom_bar(stat = "identity", color= "black", size = 0.2) +
  theme_classic(base_size = 12) + scale_fill_manual(values = rev(cols$Color), guide = F) + 
  facet_grid(. ~ phylogroup, scales='free',switch = "x", space = "free_x")
B

summarised$desc = factor(summarised$desc, all_order)
B = ggplot(summarised, aes(x = desc, y = count, label = cluster, color = phylogroup)) + geom_text(size = 5, hjust = 0, nudge_x = 0.05) + geom_point()+
  scale_color_brewer(palette = "Dark2") + theme_bw(base_size = 14)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Gene class") + ylab("Count")
B


summarised_filtered = summarised[which(summarised$desc %in% c("Real core","Missing in one","25-39 specific","10-24 specific","Core and specific",
                                                             "Intermediate and specific", "Multicluster rare","Rare and specific","2-9 varied","Intermediate and rare",
                                                             "Secondary (specific)")),]
ggplot(summarised_filtered, aes( x = phylogroup, y = count, label = cluster, color = phylogroup)) + geom_text(size = 5, nudge_x = -0.08) +
  facet_grid(~desc, scales = "free",  space = "free",switch = "x",  labeller = label_wrap_gen(width=5)) + scale_color_brewer(palette = "Dark2", guide = F) + 
  theme_minimal(base_size = 14)+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

summarised_all = aggregate(by = list(summarised$desc), x = summarised$count, mean)
summarised_all = cbind(summarised_all, col = classification$V3[match(summarised_all$Group.1, classification$V2)])
summarised_all$Group.1 = factor(summarised_all$Group.1, summarised_all$Group.1)
summarised_all$col = as.character(summarised_all$col)

total_core = round(sum(summarised_all$x[11:13]) / sum(summarised_all$x) * 100, digits = 0)
total_lineage_specific = sum(summarised_all$x[8:10])/ sum(summarised_all$x) * 100
total_varied = sum(summarised_all$x[4:7])/ sum(summarised_all$x) * 100
total_rare_inter  = sum(summarised_all$x[1:3])/ sum(summarised_all$x) * 100

B1 = ggplot(summarised_all,aes(fill = summarised_all$Group.1, x= "", y = summarised_all$x))+ geom_bar(stat = "identity", color = "black") + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = summarised_all$col, guide = F) +  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()) +
  theme(axis.text.x=element_blank()) 
B1
legend = as_ggplot(get_legend(A))
A = A + theme(legend.position = "None")
A = A + ggtitle("A")
B1 = B1 + ggtitle("B")
B = B + ggtitle("C")
grid.arrange(A,B1,B, layout_matrix = rbind(c(1,2),
                                               c(3,3),
                                               c(3,3)))

             