library(tidyverse)
library(readxl)
library(dplyr)
library(reshape2) # need to be installed first

# load RNA data from excel
RNAexp <- read_excel("_1256271tableS2.xlsx", skip = 1)

# combine unique id number with gene name
RNAexp$gene_id <- paste(RNAexp$UNIQUD, "-", RNAexp$NAME)

# relocate gene_id from the last to the first column so can make it into row id
RNAexp <- RNAexp %>% 
  relocate(gene_id) %>% 
  column_to_rownames("gene_id")

# remove columns that we dont care about
RNAexpNarrow <- select(RNAexp, "LT-HSC":MPP, CMP:Mono)

# take out the row if all the variable in the row is 0
RNAexpNarrow_filtered <- filter_all(RNAexpNarrow, any_vars(. != 0))

# Below is modded based on this website
# see website for explanation http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

cormat <- round(cor(RNAexpNarrow_filtered),2)

melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white")
  
  