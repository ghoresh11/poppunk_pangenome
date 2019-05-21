library(ggplot2)
library(RColorBrewer)
library(reshape2)

setwd("/Users/gh11/poppunk_pangenome/X_plasmids/")

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
for (i in 1:39){
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
for (i in 1:39) {
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

## the core and soft-core plasmids:
plasmid_freqs$plasmid[which(plasmid_freqs$freq>0.90)]


### Look at each replicon seperately and look at its frequency in each cluster
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
      ## count how many times they co-occur or are mutually exclusive
      both_present = length(which(plasmids[,i] == 1 & plasmids[,j] == 1))
      mutually_exclusive_A = length(which(plasmids[,i] == 1 & plasmids[,j] != 1))
      mutually_exclusive_B = length(which(plasmids[,j] == 1 & plasmids[,i] != 1))
      neither = length(which(plasmids[,i] == 0 & plasmids[,j] == 0))
      
     plt_df = data.frame(label = rep(c("both present", "mutually_exclusive", "neither present"), 2),
                         fill = c("both", "plasmidA_present", "both", "both", "plasmidB_present", "both"),
                         values = c(both_present, mutually_exclusive_A, neither, 0, mutually_exclusive_B, 0))
      p = ggplot(plt_df, aes(x = label, y = values, fill = fill)) + geom_bar(stat = "identity", color = "black") + 
        theme_classic(base_size = 16) + xlab("") + ylab("count") + 
        ggtitle(paste(colnames(plasmids)[i], colnames(plasmids)[j], sep = "\n"))  +
        scale_fill_brewer(palette = "Greys")
      print(p)
    }
  }
}



