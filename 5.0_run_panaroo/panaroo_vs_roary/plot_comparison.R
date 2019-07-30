library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(vegan)

setwd("/Users/gh11/panaroo/")

t = read.table("summary.txt", sep = ",", header = T, stringsAsFactors = F)

order = rev(c("core", "soft_core", "inter", "rare"))
cols = brewer.pal(n=5, "Blues")[-1]

t$Type = factor(t$Type, order)
t$Name = factor(t$Name, rev(c("Roary", "strict", "moderate", "relaxed")))

ggplot(t, aes(x = Name, y = Count, fill = Type)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) + theme_bw(base_size = 16) +
  xlab("") + ylab("Num genes")

ggplot(t, aes(fill = Name, y = Count, x = Type)) + geom_bar(stat = "identity", position = "dodge", color = "black") + theme_bw(base_size = 16) +
  xlab("") + ylab("Num genes") + scale_fill_brewer(palette = "Set2")

df = data.frame(type = character(0),
                richness = numeric(0),
                genomes = numeric(0),
                sd = numeric(0))

dirs = c("relaxed","moderate", "strict", "roary")
for (d in dirs){
  f = paste(d, "gene_presence_absence.Rtab", sep = "/")
  mydata <- data.frame(t(read.table(f, header = T, row.names = 1, comment.char = "", quote = )))
  sp <- specaccum(mydata, "random", permutations=100)
  
  df = rbind(df, data.frame(type = rep(d, length(sp$sites)),
                            richness = sp$richness,
                            genomes = sp$sites,
                            sd = sp$sd))
}
df = cbind(df, min= df$richness-df$sd, max = df$richness+df$sd )
df$Type = factor(t$Type, c("Roary", "strict", "moderate", "relaxed"))
ggplot(df, aes(x = genomes, y = richness, color = type))+ geom_line(alpha = 0.8)+ geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 12) + xlab("Genomes") +
  ylab("Genes") + scale_color_brewer(palette = "Set2") +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                position=position_dodge(0.05))

### check how this compares across ALL poppunk clusters
df = read.table("roary_vs_panaroo.csv", sep = ",", header = T, stringsAsFactors = F, comment.char = "")
plots = list()
for (type in c("core", "soft_core", "inter", "rare")){
  df_core = df[which(df$gene_type == type),]
  df_core$cluster = factor(df_core$cluster, unique(df_core$cluster[order(df_core$size, decreasing = T)]))
  p = ggplot(df_core, aes(x = cluster, y = count, fill = tool)) + geom_bar(stat = "identity", position = "dodge", color = "black") +
    ggtitle(type) + theme_classic(base_size = 14) + scale_fill_brewer(palette = "Set3")
  plots[[type]] = p
}
do.call("grid.arrange", c(plots, ncol=1))

roary = df[which(df$tool == "roary"),]
panaroo = df[which(df$tool == "panaroo"),]
merged = merge(panaroo, roary, by = c("cluster","gene_type"))
merged = cbind(merged, diff = merged$count.x - merged$count.y)
merged$gene_type = factor(merged$gene_type, c("core","soft_core","inter","rare"))
ggplot(merged, aes(x = gene_type, y = diff, fill = gene_type)) + geom_violin() + geom_boxplot(width = 0.5) +
  theme_bw(base_size = 16) + geom_hline(yintercept = 0, color = "#ff6700", lty = 2, size = 1.4) +
  scale_fill_manual(values = rev(brewer.pal(n = 5, "Blues")[-1]), guide = F) +
  xlab("") + ylab("#panaroo - #roary")

resource_summary = read.table("resource_summary.csv", sep = ",", header = T, stringsAsFactors = F)
ggplot(resource_summary, aes(x = Size, y = Time.mins.)) + geom_point(size = 3, alpha = 0.6) +
  theme_bw(base_size = 16) + ylab("Minutes")
ggplot(resource_summary, aes(x = Size, y = Memory.MB.)) + geom_point(size = 3, alpha = 0.6) +
  theme_bw(base_size = 16) + xlab("Memory (MB)")
ggplot(resource_summary, aes(x = Num_genes, y = Memory.MB.)) + geom_point(size = 3, alpha = 0.6) +
  theme_bw(base_size = 16) + xlab("Memory (MB)")
ggplot(resource_summary, aes(x = Num_genes, y = Time.mins.)) + geom_point(size = 3, alpha = 0.6) +
  theme_bw(base_size = 16)+ ylab("Minutes")


lm(Memory.MB. ~ Size , data = resource_summary)

## filter network
g = read.table("edgelist", sep = " ", header = F, stringsAsFactors = F)
strings = paste(g[,3], g[,4],g[,5], g[,6], sep="")
summary = data.frame(table(strings))
ggplot(summary, aes(x = strings, y = Freq)) + geom_bar(stat = "identity") +
  ggtitle("Order: roary-strict-moderate-relaxed") + ylab("Number of interactions") + xlab("Presence of interaction")

### zoomed
ggplot(summary, aes(x = strings, y = Freq)) + geom_bar(stat = "identity") +
  ggtitle("Order: roary-strict-moderate-relaxed") + ylab("Number of interactions") + xlab("Presence of interaction") +
  scale_y_continuous(limits = c(0, 88096))


t2 = read.table("summary_diff_edges_updated.csv", sep = ",", header = T, colClasses = c("character", "numeric", "numeric"))
summary_refound = data.frame(table(t2$edge_type, t2$refound))
ggplot(summary_refound, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat = "identity") +
  xlab("edge type") + ggtitle("roary-strict-moderate-rare") + scale_fill_brewer(palette = "Set2", labels = c("annotated","refound"), name = "")

t2 = t2[which(t2$edge_type %in% c("0111","1000")),]
ggplot(t2, aes(x = edge_type, y = length_diff)) + geom_violin() + geom_boxplot(width = 0.5)

t3 = t2[which(t2$edge_type == "1000"),]
t3$alignment_score = as.numeric(t3$alignment_score)
ggplot(t3, aes(x = length_diff, y = alignment_score)) + geom_point() +
  ylab("local nucleotide identity")

t4 = t2[which(t2$edge_type == "0111"),]
t4$alignment_score = as.numeric(t4$alignment_score)
ggplot(t4, aes(x = length_diff, y = alignment_score)) + geom_point() + 
  ylab("local nucleotide identity")


## comparisons of different panaroo networks
graph_comps = read.table("graph_comparisons.csv", sep = ",", header = T, stringsAsFactors = F)
graphs_1 = {}
graphs_2 = {}
graphs_3 = {}
i = 1
for (t in unique(graph_comps$Threshold)) {
  percent_5 = graph_comps[which(graph_comps$Threshold == t),]
  percent_5 = cbind(percent_5, continous_regions = percent_5$Num_CCs - percent_5$Singletons)
  percent_5$Cluster =factor(percent_5$Cluster, percent_5$Cluster[order(percent_5$Size, decreasing = T)]) 
  p1 = ggplot(percent_5, aes(x = Cluster, y = continous_regions)) + geom_bar(stat = "identity") +
    theme_bw(base_size = 14) + ylab("Number of continous regions") + ggtitle(t)
  graphs_1[[i]] = p1
  p2 = ggplot(percent_5, aes(x = Cluster, y = Mean_degree))+ geom_bar(stat = "identity") +
    theme_bw(base_size = 14) + ylab("Mean degree") + ggtitle(t)
  graphs_2[[i]] = p2
  p3 = ggplot(percent_5, aes(x = Cluster, y = Singletons))+ geom_bar(stat = "identity") +
    theme_bw(base_size = 14) + ylab("Num singletons") + ggtitle(t)
  graphs_3[[i]] = p3
  
  ggplot(percent_5, aes(x = Size, y = Singletons))
  i = i + 1
}
do.call("grid.arrange", c(graphs_1, ncol=1))
do.call("grid.arrange", c(graphs_2, ncol=1)) ## the lower the degree the more broken up the graph is
do.call("grid.arrange", c(graphs_3, ncol=1)) ## singletons are like the rare variants - see more of them when the cluster is bigger


ggplot(graph_comps, aes(x = Num_CCs, y = Mean_degree)) + geom_point()





