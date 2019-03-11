library(ggplot2)
library(RColorBrewer)
library(reshape2)

setwd("/Users/gh11/poppunk_pangenome/")

plasmids = read.table("plasmid_summary_all.csv", sep = ",", header = T, stringsAsFactors = F,
                      comment.char = "", quote = "")


plasmids = cbind(plasmids, total = rowSums(plasmids[,-c(1,2)]))

plasmids$cluster = factor(plasmids$cluster, 1:39)

### How many plasmids are there per strain? Does this change between clusters?
ggplot(plasmids, aes(x = total)) + geom_histogram(bins = 20) +
  theme_bw(base_size = 16) + xlab("Plasmids per strain") + ylab("Count")
median(plasmids$total)
mean(plasmids$total)

ggplot(plasmids, aes(x = cluster, y = total)) + geom_violin(fill = "#d3d3d3") + geom_boxplot(width = 0.1) +
  theme_classic(base_size = 16) + xlab("Cluster") + ylab("plasmids per strain")


### plasmid frequencies
## Get the frequency of each plasmid in each cluster
plasmid_freqs = data.frame(cluster = character(0), plasmid = character(0), freq = numeric(0))
for (i in 2:39){
  curr_cluster = plasmids[which(plasmids$cluster == i),]
  freqs = colSums(curr_cluster[,-c(1,2,dim(curr_cluster)[2])]) / dim(curr_cluster)[1]
  plasmid_freqs = rbind(plasmid_freqs, 
                        data.frame(cluster = rep(i, length(freqs)),
                                   plasmid = names(freqs),
                                   freq = freqs))
}


## count how many plasmids are core, soft_core, inter and rare within each cluster
## Do I find plasmids that are core to the collection or core in a specific group?
replicon_class = data.frame(cluster = character(0), type = character(0), cnt = numeric(0), stringsAsFactors = F)
for (i in 2:39) {
  curr_cluster = plasmid_freqs[which(plasmid_freqs$cluster == i),]

  cnt_rare = length(which(curr_cluster$freq <= 0.15 & curr_cluster$freq > 0))
  cnt_inter = length(which(curr_cluster$freq <= 0.90 & curr_cluster$freq > 0.15))
  cnt_soft_core = length(which(curr_cluster$freq <= 0.95 & curr_cluster$freq > 0.90))
  cnt_core = length(which(curr_cluster$freq > 0.95))
  
  replicon_class = rbind(replicon_class, 
                         data.frame(cluster = rep(i, 4), 
             type = c("rare", "inter", "soft_core", "core"),
             cnt = c(cnt_rare, cnt_inter, cnt_soft_core, cnt_core),
             stringsAsFactors = F))
}
replicon_class$cluster = factor(replicon_class$cluster, 1:39)
variable_order = c("rare", "inter", "soft_core", "core")
cols = brewer.pal(n = 5, "Blues")[-1]
labs = c("Rare", "Intermediate", "Soft core", "Core")

replicon_class$type = factor(replicon_class$type, variable_order)
ggplot(replicon_class, aes(x = cluster, y = cnt, fill = type)) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = cols, labels = labs) + theme_classic(base_size = 16) +
  scale_y_continuous(expand = c(0,0)) + xlab("Cluster") + 
  ylab("Plasmid replicons") 

plasmid_freqs$plasmid[which(plasmid_freqs$freq>0.90)]

ggplot(plasmid_freqs, aes(x = cluster, y = freq, fill = plasmid)) + geom_bar(stat = "identity")
sum(plasmid_freqs$freq[plasmid_freqs$cluster == 3])


for (p in unique(plasmid_freqs$plasmid)){
  curr = plasmid_freqs[plasmid_freqs$plasmid == p,]
  curr$cluster = factor(curr$cluster, 1:39)
  plt = ggplot(curr, aes(x = cluster, y = freq)) + geom_bar(stat = "identity") +
    theme_classic(base_size = 16) + scale_y_continuous(limits = c(0,1)) + ggtitle(p)
  print(plt)
}

### co occurence of plasmids / anti-occurent
## Are there plasmids that are mutually exclusive or always found together?
req_pval = 0.05/1378
for (i in 3:(dim(plasmids)[2]-2)){
  for (j in (i+1):(dim(plasmids)[2]-1)){
    test = fisher.test(plasmids[,i], plasmids[,j])
    if (test$p.value < req_pval) {
      df = data.frame(first = plasmids[,i], second =  plasmids[,j])
      first_zero = length(which(df$first == 0)) / dim(plasmids)[1]
      first_one = length(which(df$first == 1)) / dim(plasmids)[1]
      second_zero = length(which(df$second == 0)) / dim(plasmids)[1]
      second_one = length(which(df$second == 1)) / dim(plasmids)[1]
      
      if (first_zero >= 0.9 || second_zero >= 0.9){ next } ## I don't care about these (the majority...)
      if (first_one >= 0.9 || second_one >= 0.9){ next }
      
      df = data.frame(name = rep(c(colnames(plasmids)[i], colnames(plasmids)[j]),2),
                      val = c(0,0,1,1), count = c(first_zero, second_zero, first_one, second_one))
      df$val = factor(df$val, c(0,1))
      plt = ggplot(df, aes(x = name, y = val, fill = count)) + geom_tile(color = "black") +
        theme_bw(base_size = 16) + scale_fill_gradient(low = "white", high = "black", limits = c(0,1)) +
        ggtitle(paste(test$p.value)) 
      print(plt)
    }
  }
}



