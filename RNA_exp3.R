library(DESeq2)
library(tidyverse)
# install.packages("Hmisc") PLEASE INSTAL THIS AND COMMENT IT OUT


id_frame<- read_csv("samplesNameFile.csv")
id_frame$sampleNames <- paste(id_frame$cell_type, id_frame$replicate, sep = "_")

directory <- "rna_proj"
# directory is where your htseq-count output files are located.

sampleFiles <- paste0(file.path(directory, c(list.files(directory))))
# samplesFiles is a variable which points to your htseq-count output files,
file.exists(sampleFiles)
#checking if files in sampleFiles exist

# creating dataframe for DESeq2
sampleNames <- id_frame$sampleNames
condition <- id_frame$cell_type
## One for one for your sample type
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design = ~ condition)

keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq <- ddsHTSeq[keep,]
# pre-filter out genes that have 0 read count across the cell type

dds <- DESeq(ddsHTSeq)
#collapse biological replicates
collapse <- collapseReplicates(dds, condition)
count_collapse <- collapse@assays@data@listData[["counts"]]

library(Hmisc)

matrix <- rcorr(as.matrix(count_collapse), type = "spearman")
print(matrix)
cormat <- matrix[["r"]]
col.order <- c("LTHSC_1","STHSC_1","CMP_1","GMP_1", "BMM_1", "GN_1", "MONO_1")
cormat <- cormat[col.order,rev(col.order)]
rownames(cormat) <- c("LT-HSC","ST-HSC","CMP","GMP", "Macrophage", "Granulocyte", "Monocyte")
colnames(cormat) <- rev(c("LT-HSC","ST-HSC","CMP","GMP", "Macrophage", "Granulocyte", "Monocyte"))
# 
# library(reshape2)
# melted_cormat <- melt(cormat)

# library(corrplot)
# corrplot_try <- cor(cormat)
# col <- colorRampPalette(c("blue", "white", "red")) (20)

library(ComplexHeatmap)
pdf("dendrogram.pdf")
# Plot the result
Heatmap(corrplot_try, 
        name = "Correlation", column_names_rot = 45) #title of legend
dev.off()

# heatmap(x = corrplot_try, col = col, symm = TRUE)

# library(ggplot2)
# cols <- rev(rainbow(7)[-7])
# ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(colour="black",size=0.25) +
#   labs(x="",y="") +
#   scale_fill_gradientn(colours = cols) + #, limits = c(0.4, 1)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# ggsave("heatmap.pdf")
# 
# 
# # Euclidean distance
# dist <- dist(cormat[ , c(4:7)] , diag=TRUE)
# 
# # Hierarchical Clustering with hclust
# hc <- hclust(dist)
# 
# pdf("dendrogram.pdf")
# # Plot the result
# plot(hc,
#      main="RNAseq Cluster Dendrogram",
#      ylab="",
#      xlab="",
#      axes=F,
#      sub="")
# dev.off()

