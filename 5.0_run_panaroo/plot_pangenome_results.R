library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(vegan)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/5.0_run_panaroo/")

cluster_sizes = read.table("../2_dists_roary_analysis/cluster_sizes_updated.csv", sep = ",",
                           header = T, stringsAsFactors = F)
cluster_order = cluster_sizes$Cluster[order(as.numeric(cluster_sizes$Size), decreasing = T)]


variable_order = c("rare", "inter", "soft_core", "core")
cols = brewer.pal(n = 5, "Blues")[-1]
labs = c("Rare (<15%)", "Intermediate (15%-95%)", "Soft Core (95%-99%)", "Soft Core (>99%)")

### create the panaroo summary file from the presence absence files ###
presence_absence_files = list.files("presence_absence", pattern = "*.Rtab", full.names = T)
panaroo_summary = data.frame(Cluster = character(0),
                             Type = character(0),
                             Count = numeric(0))
for (f in presence_absence_files) {
  curr_cluster = strsplit(strsplit(f, split = "/", fixed = T)[[1]][2], split = "_", fixed = T)[[1]][1]
  curr = read.table(f, stringsAsFactors = F, comment.char = "", header = T, row.names = 1)
  curr = rowSums(curr) / dim(curr)[2]
  core = length(which(curr >= 0.99))
  soft_core = length(which(curr >= 0.95 & curr < 0.99))
  inter = length(which(curr >= 0.15 & curr < 0.95))
  rare = length(which(curr < 0.15))
  panaroo_summary = rbind(panaroo_summary, 
                          data.frame(Cluster = rep(curr_cluster, 4),
                                     Type = variable_order,
                                     Count = c(rare, inter, soft_core, core)))
}

panaroo_summary$Type = factor(panaroo_summary$Type, variable_order)
panaroo_summary$Cluster = factor(panaroo_summary$Cluster, cluster_order)
ggplot(panaroo_summary, aes(x = Cluster, y = Count, fill = Type)) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = cols, labels = labs, name = "") + theme_classic(base_size = 12) +
  scale_y_continuous(expand = c(0,0)) + xlab("Cluster") + 
  ylab("Genes") + ggtitle("C") + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=3, byrow = T))


## size point plot
sizes = cluster_sizes$Size[match(panaroo_summary$Cluster, cluster_sizes$Cluster)]
panaroo_summary = cbind(panaroo_summary, size = sizes)

ggplot(panaroo_summary, aes(y = Count, x = size, fill = Type)) + 
  geom_point(size = 3, pch=21, color = "black", alpha = 0.7) +
  theme_bw(base_size = 12) + xlab("Cluster size") + ylab("Genes") +
  scale_fill_manual(values = cols, guide = F) +  geom_smooth(method = "gam", formula = y ~ x, aes(color = Type)) +
  scale_color_manual(values = cols, guide = F) + ggtitle("D")

## compare Roary to Panaroo 
## count how many genes there are in panaroo and roary and draw a dot plot of one against the other
## need to run "analyse_dists.R"s
## here: run analyse_dists.R to get the roary outputs
merged = merge(panaroo_summary, roary_outputs, by.x = c("Cluster", "Type"), by.y = c("cluster", "variable"), 
      all.x = T, all.y = F)
colnames(merged) = c("cluster", "type","panaroo_count","rm1","rm2","roary_count","size")
merged = merged[,-c(4,5)]
ggplot(merged, aes(x = panaroo_count, y = roary_count, fill = type)) + geom_point(size = 3, color = "black", pch = 21, alpha = 0.7) +
  theme_bw(base_size = 16) + scale_fill_manual(values = cols)


# ### accumulation curves for each cluster of panaroo compared with roary
# ## gene accumilation curves using vegan
# df = data.frame(cluster = character(0),
#                 richness = numeric(0),
#                 genomes = numeric(0),
#                 sd = numeric(0))
# for (f in presence_absence_files){
#   cluster = strsplit(basename(f), split = "_", fixed = T)[[1]][1]
#   print(cluster)
#   mydata <- data.frame(t(read.table(f, header = T, row.names = 1, comment.char = "", quote = )))
#   sp <- specaccum(mydata, "random", permutations=100)
#   df = rbind(df, data.frame(cluster = rep(cluster, length(sp$sites)),
#                                  richness = sp$richness,
#                                  genomes = sp$sites,
#                                  sd = sp$sd))
# }
# # save the output
# write.table(df, "panaroo_accumilation_curves.csv", col.names = T, row.names = F, quote = F, sep = ',')

df = read.table("panaroo_accumilation_curves.csv", header = T, comment.char = "", sep = ",", stringsAsFactors = F)
df = cbind(df, min= df$richness-df$sd, max = df$richness+df$sd )
df$cluster = factor(df$cluster,cluster_sizes$Cluster)
df = df[which(df$genomes %% 12 == 0),] ## to visualise more clearly

graphics = read.table("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv", sep = ",",
                      comment.char = "", stringsAsFactors = F, header = T)
graphics = graphics[match( cluster_sizes$Cluster,graphics$Cluster),]
graphics = graphics[-1,]
ggplot(df, aes(x = genomes, y = richness, color = cluster, shape = cluster))+ geom_line(alpha = 0.8)+ geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 12) + xlab("Genomes") +
  ylab("Genes") +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                position=position_dodge(0.05)) +
  scale_color_manual(values = graphics$Color, guide = F) + scale_shape_manual(values = graphics$Shape, guide = F) + ggtitle("B") 


## plot the accumulation curves of the same cluster in a single plot
df2 = read.table("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis/accumilation_curves.csv", 
                 header = T, comment.char = "", sep = ",", stringsAsFactors = F)
df2 = cbind(df2, min= df2$richness-df2$sd, max = df2$richness+df2$sd )
df2$cluster = factor(df2$cluster,cluster_sizes$Cluster)
df2 = df2[which(df2$genomes %% 12 == 0),] ## to visualise more clearly

plots = list()
index = 1
for (curr_cluster in unique(df$cluster)){
  if (cluster_sizes$Size[cluster_sizes$Cluster == curr_cluster]<50){ next }
  first = df2[df2$cluster == curr_cluster, ]
  first = cbind(first, type = rep("roary", dim(first)[1]))
  second = df[df$cluster == curr_cluster, ]
  second = cbind(second, type = rep("panaroo", dim(first)[1]))
  curr = rbind(first, second)
  p = ggplot(curr, aes(x = genomes, y = richness, color = type))+ geom_line(alpha = 0.8)+ geom_point(size = 2, alpha = 0.8) +
    theme_bw(base_size = 12) + xlab("Genomes") +
    ylab("Genes") +
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                  position=position_dodge(0.05)) + ggtitle(curr_cluster) + scale_color_brewer(palette = "Set1", guide = F)
  plots[[index]] = p
  index = index + 1
}

legend = as_ggplot(get_legend(ggplot(curr, aes(x = genomes, y = richness, color = type))+ geom_line(alpha = 0.8)+ geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 12) + xlab("Genomes") +
  ylab("Genes") +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                position=position_dodge(0.05)) + ggtitle(curr_cluster) + scale_color_brewer(palette = "Set1")))
plots[[index]] = legend
n = length(plots)
nCol = floor(sqrt(n))
do.call("grid.arrange", c(plots, ncol=5))

