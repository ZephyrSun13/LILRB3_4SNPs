
Args <- commandArgs(trailingOnly = TRUE)

library(fgsea)
library(ggplot2)
library(data.table)
library(RColorBrewer)

ListF <- Args[1] ## Gene list (ranked)
gmtF <- Args[2]  ## GMT file
KeyS <- Args[3]  ## Key word
ThresP <- as.numeric(Args[4]) ## P.adj value threshold
GSEAParaI <- as.numeric(Args[5]) ## 1

plotEnrichmentSun <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) 
{
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
        returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    Colors <- brewer.pal(12, "Set3")
    g <- ggplot(toPlot, aes(x = x, y = y)) + 
	geom_point(color = Colors[1], size = 0.1) + 
	geom_hline(yintercept = max(tops), colour = Colors[4], linetype = "dashed") + 
	geom_hline(yintercept = min(bottoms), colour = Colors[4], linetype = "dashed") + 
	geom_hline(yintercept = 0, colour = "black") + 
	geom_area(fill = Colors[1]) + 
	theme_bw() + 
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
        theme(panel.border = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 23)) + 
        labs(x = "rank", y = "enrichment score")
    g
}

ranks <- read.table(ListF, header=TRUE, colClasses = c("character", "numeric"))

names(ranks) <- c("ID", "t")

ranks <- setNames(ranks$t, ranks$ID)

str(ranks)

pathways <- gmtPathways(gmtF)

str(head(pathways))

fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize  = 15, maxSize  = 500, nperm = 10000, gseaParam = GSEAParaI)

saveRDS(fgseaRes, file = paste0(KeyS, "_fgseaRes.rds"))

head(fgseaRes[order(pval), ])

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]

topPathways <-c(topPathwaysDown, rev(topPathwaysUp))

pp <- plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)

ggsave(paste0(KeyS, "_GSEATab.pdf"), pp, width = 12)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < ThresP], pathways, ranks)

mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

pp <- plotGseaTable(pathways[mainPathways], ranks, fgseaRes, gseaParam = 0.5)

ggsave(paste0(KeyS, "GSEATabMain.pdf"), pp, width = 12)

fwrite(fgseaRes, file=paste0(KeyS, "_fgseaRes.txt"), sep="\t", sep2=c("", " ", ""))

