library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(ggpubr)
library(dplyr)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

complete_presence_absence = fread("../4_pairwise_roary/231019_corrected//complete_presence_absence.csv", sep = ",", header = T, stringsAsFactors = F)
classification = read.table("gene_classification.csv", sep = ",", header = F, stringsAsFactors = F, comment.char = "", quote = "")

strains = colnames(complete_presence_absence)[-1]
genes = complete_presence_absence[,1][-1]
clusters = unlist(complete_presence_absence[1,])[-1]

### Section to generate the random sampling
# for (rep in 1:20){
#   print("rep")
#   print(rep)
#   random_samples = c()
#   for (curr in unique(clusters)) {
#     samples = which(clusters == curr)
#     samples = sample(x = samples ,size = 20, replace = T)
#     random_samples = c(random_samples, samples)
#   }
#   res =  data.frame( cluster = numeric(0), desc = character(0), count = numeric(0), stringsAsFactors = F)
#   for (i in random_samples){
#     i = i+1
#     cluster = unlist(complete_presence_absence[1,..i])
#     vec = complete_presence_absence[,..i][-1]
#     vec = unlist(genes[which(vec>0)])
#     for (class in unique(classification$V2)) {
#       curr_count = length(intersect(classification$V1[which(classification$V2 == class)], vec))
#       res = rbind(res, data.frame( cluster = cluster, desc = class, count = curr_count, stringsAsFactors = F))
#     }
#   }
#   write.table(res, file = paste(rep, "_class_per_isolate.csv", sep = ""), sep = ",", row.names = F, col.names = T, quote = F) ## randomly using only 15 per cluster
# }

## add the phylogroup information from the graphics file
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "")

read_one_file <- function(rep_num){
  res = read.table(file = paste(rep_num,"_class_per_isolate.csv", sep = ""), sep = ",", header = T, stringsAsFactors = F, comment.char = "", quote = "")
  summarised = aggregate(.~cluster+desc, res, mean)
  return(summarised)
}

all = read_one_file(1)
for (num in 2:20){
  
  all = rbind(all, read_one_file(num))
}

summarised_all= all %>%
  group_by(cluster, desc) %>% 
  summarise_each(list(mean = mean, sd = sd))
summarised_all = data.frame(summarised_all, stringsAsFactors = F) 
summarised_all = cbind(summarised_all, phylogroup= graphics$Phylogroup[match(summarised_all$cluster, graphics$Cluster)])

is_outlier <- function(x) {
  ret_val = rep(0, length(x))
  ret_val[which(x < quantile(x, 0.25) - 1.5 * IQR(x) )] = 1.2
  ret_val[x > quantile(x, 0.75) + 1.5 * IQR(x)] = -0.5
  return(ret_val)
}

summarised_all = summarised_all %>%
  group_by(desc, phylogroup) %>%
  mutate(vjust = is_outlier(mean))


summarised_all$outlier = summarised_all$cluster 
summarised_all$outlier[which(summarised_all$vjust == 0)] = NA

summarised_all$cluster = factor(summarised_all$cluster, as.character(unique(unlist(summarised_all$cluster))))

all_order = c("Real core", "Ubiquitous (soft-core)","Missing in one","40-45 specific",
              "25-39 specific","10-24 specific","2-9 specific","Core and specific",
              "40-46 varied","25-39 varied", "10-24 varied","2-9 varied", "Intermediate and rare",
              "Multicluster intermediate", "Intermediate and specific",
              "Multicluster rare","Rare and specific",
              "Secondary (specific)","Secondary (varied)", "Secondary (intermediate)",
              "Secondary (rare)")

summarised_all$desc = factor(summarised_all$desc, rev(all_order))
cols = read.table("descs_template.csv", sep = ",", comment.char = "",
                  stringsAsFactors = F, header = T)

summarised_all = cbind(summarised_all, outlier = rep(NA, dim(summarised_all)[1]))
summarised_all$cluster = as.character(unlist(summarised_all$cluster))


# ggplot(summarised_all, aes(x = cluster, y = mean, fill = desc)) + geom_bar(stat = "identity", color= "black", size = 0.2) +
#   theme_classic(base_size = 12) + scale_fill_manual(values = rev(cols$Color), guide = F) 
# #  facet_grid(. ~ phylogroup, scales='free',switch = "x", space = "free_x")
# 

summarised_all$desc = factor(summarised_all$desc, all_order)
summarised_all$phylogroup = factor(summarised_all$phylogroup, c("B1","C","A","E","D","F","B2","U"))

# ggplot(summarised_all, aes(x = desc, y = mean, label = cluster, color = phylogroup)) + geom_text(size = 5, hjust = 0, nudge_x = 0.05) + geom_point()+
#   scale_color_brewer(palette = "Dark2") + theme_bw(base_size = 14)+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("Gene class") + ylab("Count") +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1)

## If I only want to include some of the gene classes
# summarised_filtered = summarised_all[which(summarised_all$desc %in% c("Real core","Missing in one","25-39 specific","10-24 specific","Core and specific",
#                                                               "Intermediate and specific", "Multicluster rare","Rare and specific","2-9 varied","Intermediate and rare",
#                                                               "Secondary (specific)")),]

summarised_all$Class = cols$Main.Label[match(summarised_all$desc, cols$Var1)]
plot_list = list()
for (curr_class in unique(summarised_all$Class)) {
  curr = summarised_all[summarised_all$Class == curr_class,]
  p = ggplot(curr, aes( x = phylogroup, y = mean, color = phylogroup, label = outlier, vjust = vjust)) + geom_boxplot() +   geom_text(na.rm = TRUE)+
    facet_grid(~desc, scales = "free",  space = "free",switch = "x",  labeller = label_wrap_gen(width=5)) + scale_color_brewer(palette = "Dark2", guide = F) + 
    theme_bw(base_size = 12)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
    xlab("") + ylab("Genes") +scale_y_continuous(expand = c(0.2,0,0.2,0)) +
    scale_x_discrete(expand = c(0.1,0.1))
  plot_list[[curr_class]] = p
}


summarised_typical = aggregate(by = list(summarised_all$desc), x = summarised_all$mean, mean)
summarised_typical = cbind(summarised_typical, col = classification$V3[match(summarised_typical$Group.1, classification$V2)])
summarised_typical$Group.1 = factor(summarised_typical$Group.1, rev(summarised_typical$Group.1))
summarised_typical$col = as.character(summarised_typical$col)

total_num_genes = sum(summarised_typical$x)
total_core = round(sum(summarised_typical$x[1:3]) / total_num_genes, digits = 2)
total_lineage_specific = round(sum(summarised_typical$x[4:8]) / total_num_genes, digits = 2)
total_inter = round(sum(summarised_typical$x[9:10]) / total_num_genes, digits = 2)
total_rare = round(sum(summarised_typical$x[11:12]) / total_num_genes, digits = 2)
total_varied =  round(sum(summarised_typical$x[13:17]) / total_num_genes, digits = 2)
total_secondary =  round(sum(summarised_typical$x[18:21]) / total_num_genes, digits = 2)


A = ggplot(summarised_typical,aes(fill = summarised_typical$Group.1, x= "", y = summarised_typical$x))+ geom_bar(stat = "identity", color = "black", lwd = 0.1) + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = rev(summarised_typical$col), guide = F) +  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()) +
  theme(axis.text.x=element_blank()) 

## arrange nicely on a grid
grid.arrange(A + ggtitle("A"), 
             plot_list[["Core"]] + ggtitle("B"),
             plot_list[["Specific"]] + ggtitle("C"),
             plot_list[["Varied"]] + ggtitle("D"),
             plot_list[["Intermediate"]] + ggtitle("E"),
             plot_list[["Rare"]] + ggtitle("F"),
             plot_list[["Secondary"]] + ggtitle("G"),
             layout_matrix = rbind(c(1,1,1,2,2,2,NA,NA),
                                   c(1,1,1,3,3,3,3,3),
                                   c(NA,NA,NA,4,4,4,4,4),
                                   c(5,5,6,6,7,7,7,7)))

