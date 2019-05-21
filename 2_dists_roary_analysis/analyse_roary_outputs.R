library(ggplot2)
library(RColorBrewer)
library(vegan)

## gene accumilation curves using vegan

df = data.frame(cluster = character(0),
                richness = numeric(0),
                genomes = numeric(0),
                sd = numeric(0))

files = list.files(path = "poppunk_pangenome/2_dists_roary_analysis/gene_presence_absences/", full.names = T, pattern = "*.Rtab")
for (f in files){
  cluster = strsplit(basename(f), split = "_", fixed = T)[[1]][1]
  print(cluster)
  mydata <- data.frame(t(read.table(f, header = T, row.names = 1, comment.char = "", quote = )))
  sp <- specaccum(mydata, "random", permutations=100)
  
  df = rbind(df, data.frame(cluster = rep(cluster, length(sp$sites)),
                                 richness = sp$richness,
                                 genomes = sp$sites,
                                 sd = sp$sd))
}
## save the output
write.table(df, "poppunk_pangenome/2_dists_roary_analysis/accumilation_curves.csv", col.names = T, row.names = F, quote = F, sep = ',')

df = cbind(df, min= df$richness-df$sd, max = df$richness+df$sd )
df$cluster = factor(df$cluster, 1:39)
ggplot(df, aes(x = genomes, y = richness, color = cluster)) + geom_line(size = 2,alpha = 0.6) +
  theme_bw(base_size = 16) + xlab("Genomes") +
  ylab("Genes") +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.2,
                position=position_dodge(0.05)) +
  scale_color_manual( values = c(
    brewer.pal(n = 8, "Set1"),
    brewer.pal(n = 8, "Set2"),
    brewer.pal(n = 8, "Set3"),
    brewer.pal(n = 8, "Dark2"),
    brewer.pal(n = 7, "Paired"))
  )


write.table(tree, "../7_tree/tree.txt", col.names = F, row.names = F, quote = F, sep = ",")



