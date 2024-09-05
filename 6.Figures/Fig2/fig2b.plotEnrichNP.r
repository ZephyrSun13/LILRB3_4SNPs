
Args <- commandArgs(trailingOnly = T)

EnrichTab <- Args[1]
PlotTerms <- Args[2]
ColorPlot <- Args[3]
DBs <- Args[4]

library(ggplot2)
library(RColorBrewer)

Colors <- read.delim(file = "ColorLibrary", row.names = 1, header = F, as.is = T)

Enrich <- read.delim(file = EnrichTab, as.is = T)

Enrich <- Enrich[Enrich$DB %in% unlist(strsplit(DBs, split = ";")), ]

ForPlot <- read.delim(file = PlotTerms, as.is = T, header = F)

ForPlot <- ForPlot[ForPlot[[1]] %in% Enrich$Term, ]

Enrich <- Enrich[Enrich$Term %in% ForPlot[[1]], ]

Tab <- rbind(data.frame(Term = ForPlot[match(Enrich$Term, ForPlot[[1]]), 2], P = -log10(Enrich$P.value), Rate = -log10(Enrich$P.value)/Enrich$Count*Enrich$Pos, Group = "Pos"), data.frame(Term = ForPlot[match(Enrich$Term, ForPlot[[1]]), 2], P = -log10(Enrich$P.value), Rate = -log10(Enrich$P.value)/Enrich$Count*Enrich$Neg, Group = "Neg"))

Tab$Term <- factor(as.character(Tab$Term), levels = rev(ForPlot[[2]]))

Tab$Group <- factor(Tab$Group, levels = c("Neg", "Pos"))

pp <- ggplot(Tab, aes(x = Term, y = Rate, fill = Group)) +

	geom_bar(stat="identity") +

	scale_fill_manual(values = Colors[unlist(strsplit(ColorPlot, split = ";")), ]) +

	theme_classic() +

	ylab("-log10(P value)") +

	xlab("") +

	coord_flip() +

	scale_y_continuous(expand = c(0, 0))+

	theme(legend.title = element_blank(), axis.text.y = element_text(size = 10))

ggsave(paste0(EnrichTab, PlotTerms, ".plot.pdf"), pp, width = 0.1 * max(sapply(Tab$Term, function(x) nchar(as.character(x)))) + 2, height = length(unique(Tab$Term)) * 0.2 + 0.7)

