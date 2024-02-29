
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

Tabs <- unlist(read.table(file = "Tabs", as.is = T))

names(Tabs) <- sapply(Tabs, function(x) unlist(strsplit(x, split = "/"))[10])

#SamTab <- read.delim(file = "RNAseq_samples.txt", as.is = T)
#
#rownames(SamTab) <- paste0("S", SamTab$ID)
#
#Tabs <- Tabs[names(Tabs) %in% rownames(SamTab)[SamTab$Time_point == "0 Month" & SamTab$Tissue == "Blood"]]


TabsLis <- lapply(Tabs, function(x) read.delim(file = x, as.is = T))

names(TabsLis) <- gsub("_aln.dedup.bam.RG.bam.ase.output.table", "", names(TabsLis))

#names(TabsLis) <- sapply(Tabs, function(x) unlist(strsplit(x, split = "/"))[10])

#names(TabsLis) <- gsub("_BS", "", names(TabsLis))

TabsLis <- lapply(TabsLis, function(x) {rownames(x) <- paste(paste0("chr", x$contig), x$position, x$refAllele, x$altAllele, sep = ":"); x})

saveRDS(TabsLis, file = "TabsLis.rds")

