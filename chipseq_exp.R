library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
# BiocManager::install("plyranges")
library(plyranges)

# below uses the bed files before intersect
STHSC_1 <- read_narrowpeaks("./chipseq_postintersect/STHSC/SRR1521829.narrowPeak.genes.bed")
STHSC_2 <- read_narrowpeaks("./chipseq_postintersect/STHSC/SRR1521831.narrowPeak.genes.bed")
STHSC <- GRangesList(STHSC_1, STHSC_2)
peak_grangeslist <- GRangesList(STHSC)
peak_coverage <- coverage(peak_grangeslist)
covered_ranges <- slice(STHSC_peak_coverage, lower = 2, rangesOnly = TRUE)


# below uses the bed files after intersect
STHSC_peak_files <- list.files("./chipseq_preintersect/STHSC", full.names = TRUE)
STHSC_peak_granges <- lapply(STHSC_peak_files, import)
STHSC_peak_grangeslist <- GRangesList(STHSC_peak_granges)
STHSC_peak_coverage <- coverage(STHSC_peak_grangeslist)
STHSC_covered_ranges <- slice(STHSC_peak_coverage, lower = 2, rangesOnly = TRUE)

