library(ggplot2)
library(RColorBrewer)
library(gridExtra)

setwd("/Users/gh11/poppunk_pangenome/5_classify_genes/")

## classify genes according to their distribution pattern and see if there
## is a composition of genes which are "universal" to E. coli, or specific to 
## a PopPUNK cluster
gene_classes = read.table("../4_pairwise_roary/231019_corrected//gene_classes.csv", sep = ",", comment.char = "", header = T, 
                          stringsAsFactors = F, quote = "", row.names = 1)
gene_classes$gene_class = rep("", dim(gene_classes)[1])

remove_from_gene_classes <- function(vec){
  gene_classes = gene_classes[-which(row.names(gene_classes) %in% vec),]
  return (gene_classes)
}

update_gene_classes <- function(vec, name) {
  gene_classes$gene_class[which(row.names(gene_classes) %in% vec)] = name
  return (gene_classes)
}


## writing exactly according to the decision tree:
one = gene_classes[which(gene_classes$total_presence == 1),]
core_and_specific = row.names(one)[which(one$core == 1)]
gene_classes = update_gene_classes(core_and_specific, "Core and specific")
intermediate_and_specific =  row.names(one)[which(one$inter == 1)]
gene_classes = update_gene_classes(intermediate_and_specific, "Intermediate and specific")
rare_and_specific = row.names(one)[which(one$rare == 1)]
gene_classes = update_gene_classes(rare_and_specific, "Rare and specific")


fortyseven = gene_classes[which(gene_classes$total_presence == 47),]
complete_core = row.names(fortyseven)[which(fortyseven$core == 47)]
gene_classes = update_gene_classes(complete_core, "Real core")
ubiq_almost_core = row.names(fortyseven)[which(fortyseven$core != 47)]
gene_classes = update_gene_classes(ubiq_almost_core, "Ubiquitous (soft-core)")

other = gene_classes[which(gene_classes$total_presence>1 & gene_classes$total_presence < 47),]

## 0
no_core = other[which(other$core == 0),]
rare_but_multiclade = row.names(no_core)[which(no_core$inter == 0)]
gene_classes = update_gene_classes(rare_but_multiclade, "Multicluster rare")
inter_multiclade = row.names(no_core)[which(no_core$rare == 0)]
gene_classes = update_gene_classes(inter_multiclade, "Multicluster intermediate")
inter_and_rare = row.names(no_core)[which(no_core$rare != 0 & no_core$inter != 0)]
gene_classes = update_gene_classes(inter_and_rare, "Intermediate and rare")


### Other
get_values_in_range <- function(minimum, maximum) {
  range = other[which(other$core >= minimum & other$core<= maximum),]
  specific = row.names(range)[which(range$core == range$total_presence)]
  varied = row.names(range)[which(range$core != range$total_presence)]
  return(list(specific, varied))
}

rarely = get_values_in_range(1,9)
rarely_specific = rarely[[1]]
gene_classes = update_gene_classes(rarely_specific, "2-9 specific")
rarely_varied = rarely[[2]]
gene_classes = update_gene_classes(rarely_varied, "2-9 varied")

occasional = get_values_in_range(10,24)
occasional_specific = occasional[[1]]
gene_classes = update_gene_classes(occasional_specific, "10-24 specific")
occasional_varied = occasional[[2]]
gene_classes = update_gene_classes(occasional_varied, "10-24 varied")

sometimes = get_values_in_range(25,39)
sometimes_specific = sometimes[[1]]
gene_classes = update_gene_classes(sometimes_specific, "25-39 specific")
sometimes_varied = sometimes[[2]]
gene_classes = update_gene_classes(sometimes_varied, "25-39 varied")

mostly = get_values_in_range(40,45)
mostly_specific = mostly[[1]]
gene_classes = update_gene_classes(mostly_specific, "40-45 specific")
mostly_varied = mostly[[2]]
gene_classes = update_gene_classes(mostly_varied, "40-45 varied")


missing_in_one = rownames(other)[which(other$core == 46)]
gene_classes = update_gene_classes(missing_in_one, "Missing in one")

secondary = gene_classes[grepl(gene_classes$truncation_assignments, pattern = "secondary"),]
secondary$gene_class[grepl(secondary$gene_class, pattern = "varied")] = "Secondary (varied)"
secondary$gene_class[grepl(secondary$gene_class, pattern = "core")] = "Secondary (core)"
secondary$gene_class[grepl(tolower(secondary$gene_class), pattern = "rare")] = "Secondary (rare)"
secondary$gene_class[grepl(tolower(secondary$gene_class), pattern = "intermediate")] = "Secondary (intermediate)"
secondary$gene_class[grepl(secondary$gene_class, pattern = "specific")] = "Secondary (specific)"

gene_classes = gene_classes[!grepl(gene_classes$truncation_assignments, pattern = "secondary"),]
gene_classes = rbind(gene_classes, secondary)


all = data.frame(table(gene_classes$gene_class), stringsAsFactors = F)
all = all[order(all$Var1),]
all$Cat =  c(rep(c("Specific Core", "Varied"),5),"Intermediate","Core","Intermediate","Rare","Rare","Core",rep("Secondary variant",4),"Core" )
all = all[order(all$Cat),]
all_order = c("Real core", "Ubiquitous (soft-core)","Missing in one","40-45 specific",
              "25-39 specific","10-24 specific","2-9 specific","Core and specific",
              "Multicluster intermediate", "Intermediate and specific",
              "Multicluster rare","Rare and specific",
              "40-45 varied","25-39 varied", "10-24 varied","2-9 varied",
              "Intermediate and rare","Secondary (specific)","Secondary (intermediate)",
              "Secondary (rare)","Secondary (varied)")
all = all[match(all_order, all$Var1),]

all$Color = c(rev(brewer.pal(n=4, "Blues")[-1]), ## core
              rev(brewer.pal(n = 6, "Purples")[-1]), ## specific core
              rev(brewer.pal(n = 3, "Reds")[-1]), ## intermediate
              rev(brewer.pal(n = 3, "Oranges")[-1]), ## rare
              rev(brewer.pal(n = 5, "Greens")), ## varied
              c("#DECBE4", "#FBB4AE", "#FED9A6" ,"#B3E2CD")) ## secondary

#all$Color[which(!all$Cat %in% c("Core","Specific Core"))] = "#d3d3d3"
#all$Color[which(!all$Cat %in% c("Rare"))] = "#d3d3d3"
all$Color[which(!all$Cat %in% c("Intermediate","Varied"))] = "#d3d3d3"

#write.table(all, file = "descs_template.csv", sep = ",", quote = F, row.names = F, col.names = T)

## get percentages for text
sum(all$Freq[all$Cat == "Core"]) / sum(all$Freq)
sum(all$Freq[all$Cat == "Specific Core"]) / sum(all$Freq)
sum(all$Freq[all$Cat == "Intermediate"]) / sum(all$Freq)
sum(all$Freq[all$Cat == "Rare"]) / sum(all$Freq)
sum(all$Freq[all$Cat == "Varied"]) / sum(all$Freq)
sum(all$Freq[all$Cat == "Secondary variant"]) / sum(all$Freq)


all$Cat = factor(all$Cat, rev(c("Core","Specific Core","Intermediate","Rare","Varied","Secondary variant")))
all$Var1 = factor(all$Var1, all$Var1 ) ## TODO change the order
# df$Count = as.numeric(df$Count)
# df$Color = as.character(df$Color)
A = ggplot(all, aes(x = Cat, fill = Var1, y= Freq)) + geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = all$Color, name = "") +
  theme_classic(base_size = 12) + ylab("Genes") + coord_flip() + scale_y_continuous(expand = c(0.01,0.01))+
  theme(legend.position = "None") + xlab("")+guides(fill=guide_legend(ncol=3,byrow=F))

A
legend = as_ggplot(get_legend(A))
legend
A = A + theme(legend.position = "None")
A

write.table(x = data.frame(row.names(gene_classes), gene_classes$gene_class, all$Color[match(gene_classes$gene_class, all$Var1)]),file = "gene_classification.csv",
            col.names = F, row.names = F, sep = ",", quote = F)

