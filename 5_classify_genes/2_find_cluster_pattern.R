library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(data.table)
library(ggpubr)
library(dplyr)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

complete_presence_absence = fread("../4_pairwise_roary/231019_corrected//complete_presence_absence.csv", sep = ",", header = T, stringsAsFactors = F)
classification = read.table("classification_v2.csv", sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")

strains = colnames(complete_presence_absence)[-1]
genes = complete_presence_absence[,1][-1]
clusters = unlist(complete_presence_absence[1,])[-1]

# ### Section to generate the random sampling
# for (rep in 2:20){
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
#     vec = unlist(genes[which(vec>0)]) ## get all genes that this sample has
#     for (class in unique(classification$fill)) {
#       curr_count = length(intersect(classification$gene[which(classification$fill == class)], vec))
#       res = rbind(res, data.frame( cluster = cluster, desc = class, count = curr_count, stringsAsFactors = F))
#     }
#   }
#   write.table(res, file = paste(rep, "_class_per_isolate.csv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F) ## randomly using only 15 per cluster
# }

## add the phylogroup information from the graphics file
graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter4/figures/cluster_graphics.csv", sep = ",",
                      header = T, comment.char = "")

read_one_file <- function(rep_num){
  res = read.table(file = paste(rep_num,"_class_per_isolate.csv", sep = ""), sep = "\t", header = T, stringsAsFactors = F, comment.char = "", quote = "")
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

colours = read.table("colours_v2.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)

summarised_all$desc = factor(summarised_all$desc, colours$Class)
summarised_all = cbind(summarised_all, outlier = rep(NA, dim(summarised_all)[1]))
summarised_all$cluster = as.character(unlist(summarised_all$cluster))


# ggplot(summarised_all, aes(x = cluster, y = mean, fill = desc)) + geom_bar(stat = "identity", color= "black", size = 0.2) +
#   theme_classic(base_size = 12) + scale_fill_manual(values = colours$Colour, guide = F)+
#   facet_grid(. ~ phylogroup, scales='free',switch = "x", space = "free_x")


summarised_all$desc = factor(summarised_all$desc, colours$Class)
summarised_all$phylogroup = factor(summarised_all$phylogroup, c("B1","C","A","E","D","F","B2","U"))


## If I only want to include some of the gene classes
# summarised_filtered = summarised_all[which(summarised_all$desc %in% c("Real core","Missing in one","25-39 specific","10-24 specific","Core and specific",
#                                                               "Intermediate and specific", "Multicluster rare","Rare and specific","2-9 varied","Intermediate and rare",
#                                                               "Secondary (specific)")),]

summarised_typical = aggregate(by = list(summarised_all$desc), x = summarised_all$mean, mean)
summarised_typical$sd = aggregate(by = list(summarised_all$desc), x = summarised_all$mean, sd)$x
summarised_typical = cbind(summarised_typical, col = colours$Colour[match(summarised_typical$Group.1, colours$Class)])
summarised_typical$Group.1 = factor(summarised_typical$Group.1, rev(summarised_typical$Group.1))
summarised_typical$col = as.character(summarised_typical$col)
summarised_typical$Main = colours$Main.Class[match(summarised_typical$Group.1, colours$Class)]
sum(summarised_typical$x)
summarised_typical$min = (summarised_typical$x - summarised_typical$sd)
summarised_typical$max = summarised_typical$x + summarised_typical$sd
summarised_typical$min_percent =summarised_typical$min/sum(summarised_typical$x)*100
summarised_typical$max_percent =summarised_typical$max/sum(summarised_typical$x)*100


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



summarised_all$Class = colours$Main.Class[match(as.character(summarised_all$desc), colours$Class)]
plot_list = list()
plot_list[["A"]] = A + ggtitle("A")

write.table(summarised_all, file = "genes_per_ecoli.csv", sep = "\t",
            col.names = T, row.names = F,quote = F)
write.table(summarised_typical, file = "typical_ecoli.csv", sep = "\t",
            col.names = T, row.names = F,quote = F)


i = 2
for (curr_class in colours$Class) {
  curr = summarised_all[summarised_all$desc == curr_class,]
  curr_min = max(0,summarised_typical$min[summarised_typical$Group.1 == curr_class])
  curr_max = summarised_typical$max[summarised_typical$Group.1 == curr_class]
  curr_mean = summarised_typical$x[summarised_typical$Group.1 == curr_class]
  curr$name = curr_class
  text_colour = "white"
  if (curr_class %in% c("Intermediate and rare","Core, intermediate and rare", "Core and rare", "Cluster specific rare", "Cluster specific intermediate")) {
    text_colour = "black"
  }
  max_colour = colours$Colour[colours$Class == curr_class]
  title = chartr("123456789", "ABCDEFGHI",i)
  if (i == 10) {
    title = "J"
  } else if (i == 11) { 
    title = "K"
  }else if  (i == 12) {
    title ="L"
  }
  p = ggplot(curr, aes( x = phylogroup, y = mean, color = phylogroup, label = outlier, vjust = vjust)) + geom_boxplot() +   geom_text(na.rm = TRUE)+
    # facet_grid(~desc, scales = "free_y",  space = "free",switch = "x",  labeller = label_wrap_gen(width=5)) + 
    scale_color_brewer(palette = "Dark2", guide = F) + 
    theme_bw(base_size = 12)+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
    xlab("") + ylab("Genes") +scale_y_continuous(expand = c(0.2,0,0.2,0)) +
    scale_x_discrete(expand = c(0.1,0.1)) + ggtitle(curr_class) +
    annotate("rect", xmin = 0, xmax = 9, ymin = curr_min, ymax = curr_max, 
             alpha = .2) +geom_hline(yintercept = curr_mean) + facet_grid(. ~ name) +
    theme(strip.background = element_rect(fill=max_colour),
          strip.text = element_text(size=12, colour=text_colour)) + ggtitle(title)
  plot_list[[curr_class]] = p
  i = i + 1
}




## arrange nicely on a grid

do.call("grid.arrange", c(plot_list, ncol=4))



