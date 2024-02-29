
library(ggplot2)

Colors <- read.delim(file = "ColorLibrary", header = F, as.is = T, row.names = 1)

Tab <- read.delim(file = "DCGL.over.xls.SigSNPTab.sele.xls.ann.xls.Gene", as.is = T, header = F)

Tab <- Tab[order(Tab[[2]], decreasing = T), , drop = F]

Tab <- Tab[1:20, ] 

Tab[[1]] <- factor(Tab[[1]], levels = Tab[[1]])

colnames(Tab) <- c("Gene", "Frequency")


        pp <- ggplot(data=Tab, aes_string(x="Gene", y="Frequency")) +

                geom_bar(stat="identity", width=0.5, fill = Colors["Blue", ]) +

		theme_classic() +

		theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), axis.title.x=element_blank())

	ggsave("GeneFreq.pdf", pp, height = 3)

