library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

id_frame<- read_csv("samples.csv")
id_frame$sampleNames <- paste(id_frame$cell_type, id_frame$replicate, sep = "_")

directory <- "rna_proj"
# directory is where your htseq-count output files are located.

sampleFiles <- paste0(file.path(directory, c(list.files(directory))))
# samplesFiles is a variable which points to your htseq-count output files,
file.exists(sampleFiles)

sampleNames <- id_frame$sampleNames

condition <- id_frame$cell_type

# One for one for your sample type

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       design = ~ condition)

keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq <- ddsHTSeq[keep,]

# pre-filter out genes that have 0 read count across the cell type

dds <- DESeq(ddsHTSeq)
res <- results(dds)

collapse <- collapseReplicates(dds, condition)
test <- collapse@assays@data@listData[["counts"]]


