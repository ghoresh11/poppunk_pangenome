

setwd("/Users/gh11/poppunk_pangenome/10_gene_sharing/sharing_on_tree/")

gene_class = "Core, intermediate and rare"


df_real = read.table(paste(gene_class,"_itol_pie.txt", sep = ""), stringsAsFactors = F,
                     sep = ",", skip = 8, comment.char = "", quote = "", header = F)

randoms_specific = read.table(paste("randoms/random_results/",gene_class,"_specific.csv", sep = ""), sep = ",",
                              header = F, comment.char = "", quote = "", stringsAsFactors = F, row.names = 1)
randoms_missing = read.table(paste("randoms/random_results/",gene_class,"_missing.csv", sep = ""), sep = ",",
                              header = F, comment.char = "", quote = "", stringsAsFactors = F, row.names = 1)


df_real$signif_specific = rep(1,dim(df_real)[1])
df_real$signif_missing = rep(1,dim(df_real)[1])

for (i in 1:length(df_real$V1)) {
  id = df_real$V1[i]
  curr_row = which(rownames(randoms_specific) == id)
  df_real$signif_specific[i] = length(which(randoms_specific[curr_row,] > df_real$V4[i]))/dim(randoms_specific)[2]
  df_real$signif_missing[i] = length(which(randoms_missing[curr_row,] > df_real$V5[i]))/dim(randoms_missing)[2]
}


label_file = paste("randoms/signif/",gene_class, "_ALL.txt",sep = "")
file_open =  file(label_file)
writeLines(c("DATASET_TEXT","SEPARATOR COMMA","DATASET_LABEL,labels","COLOR,#d3d3d3",
             "DATA"), file_open)
close(file_open)
df_out = data.frame(id = df_real$V1,
                    label = paste(df_real$signif_specific, "/", df_real$signif_missing),
                    pos = rep(0.3, length(df_real$V1)),
                    col = rep("black", length(df_real$V1)),
                    type = rep("bold", length(df_real$V1)),
                    size = rep(2, length(df_real$V1)),
                    angle = rep(270, length(df_real$V1)[1]), stringsAsFactors = F)
write.table(x = df_out, file = label_file, sep = ",", quote = F, col.names = F, row.names = F, append = T)



df = read.table("randoms/observed_vs_random.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = T)
ggplot(df, aes(x = Random, y = Observed, label =Gene.class)) + geom_text() + geom_smooth(method = "lm",se = F)

