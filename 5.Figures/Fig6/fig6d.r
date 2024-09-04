
Args <- commandArgs(trailingOnly = TRUE)

CountF <- Args[1]
OrdersS <- Args[2]
SamInfoF <- Args[3]

library(ggplot2)
library(RColorBrewer)
library(ggrepel)

Colors <- c(brewer.pal(8, "Set2")[c(-1, -6, -8)], brewer.pal(9, "Set1")[c(-1, -6, -9)], brewer.pal(12, "Set3")[c(-2)], brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))

Stat_Count <- read.delim(file = Args[1], row.names = 1, comment.char = "#", check.names = F)

Sum <- colSums(Stat_Count)

Stat_Perc <- data.frame(t(apply(Stat_Count, 1, function(x) x/Sum)), check.names = F)

Sum2 <- rowSums(Stat_Count)

Stat_Perc2 <- data.frame(apply(Stat_Count, 2, function(x) x/Sum2), check.names = F)

write.table(Stat_Count, file = "Stat_Count.xls", col.names = NA, quote = FALSE, sep = "\t")

write.table(Stat_Perc, file = "Stat_Perc.xls", col.names = NA, quote = FALSE, sep = "\t")

write.table(Stat_Perc2, file = "Stat_Perc2.xls", col.names = NA, quote = FALSE, sep = "\t")


Stat_Count <- data.frame(t(Stat_Count))

Stat_Perc <- data.frame(t(Stat_Perc))

Stat_Perc2 <- data.frame(Stat_Perc2)

Sample.info <- read.delim(file = SamInfoF, row.names = 1, as.is = T, comment.char = "#")

Tab <- list()

for(rr in rownames(Stat_Count)){

	Tab[[rr]] <- data.frame(Cell = colnames(Stat_Count), Sam = rr, Count = unlist(Stat_Count[rr, ]))


}

Tab <- Reduce(rbind, Tab)

Tab <- Tab[Tab$Sam %in% rownames(Sample.info), ]

Tab$Sam <- factor(as.character(Tab$Sam), levels = rownames(Sample.info))

Tab$Cell <- factor(as.character(Tab$Cell), levels = unlist(strsplit(OrdersS, split = ";")))

p<-ggplot(data=Tab, aes(x=Sam, y=Count, fill = Cell)) +

        geom_bar(stat="identity", width = 0.5) +

        theme_minimal() +

        scale_fill_manual(values=Colors) +

	xlab("") +

        ylab("Cell Count") +

        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.title = element_blank(), legend.text = element_text(size = 12))

ggsave("CellNumber_Distribution.pdf", p, width = 8)


p<-ggplot(data=Tab, aes(x=Cell, y=Count, fill = Sam)) +

        geom_bar(stat="identity", width = 0.5) +

        theme_minimal() +

        scale_fill_manual(values=Colors) +

        ylab("Cell Count") +

        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())

ggsave("CellNumber_Distribution2.pdf", p, width = 14)


Tab <- list()

for(rr in rownames(Stat_Perc)){

        Tab[[rr]] <- data.frame(Cell = colnames(Stat_Perc), Sam = rr, Count = unlist(Stat_Perc[rr, ]))


}

Tab <- Reduce(rbind, Tab)

Tab$Sam <- factor(as.character(Tab$Sam), levels = rownames(Sample.info))

Tab$Cell <- factor(as.character(Tab$Cell), levels = unlist(strsplit(OrdersS, split = ";")))

p<-ggplot(data=Tab, aes(x=Sam, y=Count, fill = Cell)) +

        geom_bar(stat="identity", width = 0.5) +

        theme_minimal() +

        scale_fill_manual(values=Colors) +

        ylab("Propotion of each cell type") +

	xlab("") +

        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), legend.title = element_blank(), legend.text = element_text(size = 12))

ggsave("CellPercentage_Distribution.pdf", p, width = 8)


Tab <- list()

Stat_PercCon <- apply(Stat_Perc, 2, function(x) tapply(x, Sample.info[rownames(Stat_Perc), ], mean))

for(rr in rownames(Stat_PercCon)){

        Tab[[rr]] <- data.frame(Cell = colnames(Stat_PercCon), Sam = rr, Count = unlist(Stat_PercCon[rr, ]))


}

Tab <- Reduce(rbind, Tab)

Tab$Sam <- factor(as.character(Tab$Sam), levels = Sample.info$Treatment[!duplicated(Sample.info$Treatment)])

Tab$Cell <- factor(as.character(Tab$Cell), levels = unlist(strsplit(OrdersS, split = ";")))

p<-ggplot(data=Tab, aes(x=Sam, y=Count, fill = Cell)) +

        geom_bar(stat="identity", width = 0.5) +

        theme_classic() +

        scale_fill_manual(values = Colors) +

        ylab("Propotion of each cell type") +

        xlab("") +

        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "plain"), axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 20), legend.title = element_blank(), legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"))

ggsave("CellPercentage_DistributionMean.pdf", p, width = 5, height = 6)

