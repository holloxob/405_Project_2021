library(tidyverse)

jaccard_cor <- read_csv("jaccard.csv")

jaccard_cor <- column_to_rownames(jaccard_cor, "Name")

jaccard_cor <- as.matrix(jaccard_cor)

library(reshape2)
melted_jaccard_cor <- melt(jaccard_cor)

library(ggplot2)
cols <- rev(rainbow(7)[-7])

ggplot(data = melted_jaccard_cor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(colour="black",size=0.25) +
  labs(x="",y="") +
  scale_fill_gradientn(colours = cols, limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


ggsave("chip_heatmap.pdf")


# Euclidean distance
dist <- dist(jaccard_cor[ , c(4:8)] , diag=TRUE)

# Hierarchical Clustering with hclust
hc <- hclust(dist)

pdf("chip_dendrogram.pdf")
# Plot the result
plot(hc,
     main="H3K4me1 Cluster Dendrogram",
     ylab="",
     xlab="",
     axes=F,
     sub="")
dev.off()


