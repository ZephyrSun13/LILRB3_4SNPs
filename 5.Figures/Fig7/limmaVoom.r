
Args <- commandArgs(trailingOnly = T)

ExprFile <- Args[1]
FactorFile <- Args[2]
IfPaired <- Args[3]
Key <- Args[4]
Variables <- Args[5]
ThresFC <- as.numeric(Args[6])
ThresP <- as.numeric(Args[7]) ## 0.05
IDType <- Args[8]
Key2 <- Args[9]
AdjustBy <- Args[10]
Method <- Args[11]
ThresCPM <- as.numeric(Args[12]) ## 1
AnnF <- Args[13]

library(limma)
library(multtest)
library(pheatmap)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)

Colors <- read.delim(file = "ColorLibrary", row.names = 1, header = F, as.is = T)

Expr <- read.table(file = ExprFile, sep = "\t", header = T, row.names = 1, check.names = FALSE)

if(IDType == "Ensembl"){

	rownames(Expr) <- sapply(rownames(Expr), function(x) unlist(strsplit(x, split = "\\."))[1])

}

dim(Expr)

Factor <- read.delim(file = FactorFile, header = T, stringsAsFactors = T, row.names = 1)

type <- IfPaired

Focus <- Key 

Cases <- unlist(strsplit(Variables, split = ";"))

Factor <- Factor[!is.na(Factor[[Focus]]), , drop = F]

Factor <- Factor[Factor[[Focus]] %in% Cases, ]

Factor[[Focus]] <- factor(Factor[[Focus]], levels = c(Cases[2], Cases[1]))

Expr <- Expr[, intersect(rownames(Factor), names(Expr))]

Factor <- Factor[intersect(rownames(Factor), names(Expr)), ]

print("Raw")

dim(Expr)

Expr[1:4, 1:4]

dge <- DGEList(counts = Expr)

dge <- calcNormFactors(dge)

drop <- which(apply(cpm(dge), 1, max) >= ThresCPM)

Expr.norm <- cpm(dge)

dge <- dge[drop,] 

print("Filtered")

dim(dge) # number of genes left


mds <- plotMDS(dge)

Tab <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = Factor[names(mds$x), Focus])

pp <- ggplot(Tab, aes(x = Dim1, y = Dim2, color = Group)) +

	geom_point() +

	geom_text_repel(aes(label = rownames(Tab)), size = 2) +

	theme_minimal()

ggsave(paste0(Focus, "_MDSPlot.pdf"), height = 4, width = 4)

if(AdjustBy != "No"){

	ADJ <- unlist(strsplit(AdjustBy, split = ";"))

	for(nn in ADJ){

		Tab <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = Factor[names(mds$x), nn])
		
		pp <- ggplot(Tab, aes(x = Dim1, y = Dim2, color = Group)) +
		
		        geom_point() +
		
		        geom_text_repel(aes(label = rownames(Tab)), size = 2) +
		
		        theme_minimal()
		
		ggsave(paste0(nn, "_MDSPlot.pdf"))

	}

}

print("Starting limma")

if(type=="paired"){

	pair <- factor(Factor[[Focus]], levels = c(Cases[2], Cases[1]))

	stage <- Factor$Identity

	design<-model.matrix(~pair+stage)

	fit<-lmFit(Expr, design)

	fit <- eBayes(fit, trend = F)

        result <- topTable(fit, coef=grep("pair", colnames(design), value = TRUE), adjust="BH", number = 50000)

        result <- result[order(result$P.Val), ]

        names(result) <- c("Log2Rat", "AveExpr", "t", "limma.p", "limma.padj", "B")

	write.table(result, file = paste("DEG", "diff.xls", sep="_"), sep = "\t", quote = F, col.names = NA)

	result <- result[abs(result$Log2Rat) >= as.numeric(Args[6]) & result$limma.padj <= as.numeric(Args[7]), ]

	Ann <- read.delim(file = AnnF, stringsAsFactors = F)


                if(IDType == "Symbol"){

                        result.ann <- merge(result, Ann, by.x = 0, by.y = 1, all = T)

                        colnames(result.ann)[1] <- "Symbol"

                }else if(IDType == "Ensembl"){

                        result.ann <- merge(result, Ann, by.x = 0, by.y = 2, all = T)

                        colnames(result.ann)[1] <- "Ensembl"

                }

                result.ann <- merge(result.ann, Expr, by.x = 1, by.y = 0)

        write.table(result.ann, file = paste("DEG", "diff.ann.xls", sep="_"), sep = "\t", quote = F, row.names = F)

}else{

	if(AdjustBy != "No"){

		AdjustBy <- unlist(strsplit(AdjustBy, split = ";"))

		design <- model.matrix(as.formula(paste0(" ~ 0 + ", Focus, " + ", paste(AdjustBy, collapse = " + "))), Factor)

		colnames(design)[1:2] <- levels(Factor[[Focus]])


	}else{

		design <- model.matrix(~ 0 + Factor[[Focus]])

		colnames(design) <- levels(Factor[[Focus]])

	}

	if(Method == "voom"){

		pdf("Voom.pdf")

			dge2 <- voom(dge, design, plot = T)

		dev.off()

	}

	fit <- lmFit(dge2, design)

	if(length(Cases) == 2){

		Contrast <- paste0(Cases[1], "-", Cases[2])

	}else{

		Contrast <- c(paste0(Cases[1], "-", Cases[2]), paste0(Cases[1], "-", Cases[3]), paste0(Cases[2], "-", Cases[3]))

	}

	contrast.matrix <- makeContrasts(contrasts=Contrast, levels=design)

	fit <- contrasts.fit(fit, contrast.matrix)

	fit <- eBayes(fit, trend = F)

	for(i in 1:length(Contrast)){

		result <- topTable(fit, coef=i, adjust="BH", number = 50000)

		result <- result[order(result$P.Val), ]

		names(result) <- c("Log2Rat", "AveExpr", "t", "limma.p", "limma.padj", "B")		

        	write.table(result, file = paste("DEG", "diff.xls", sep="_"), sep = "\t", quote = F, col.names = NA)

	        Ann <- read.delim(file = AnnF, stringsAsFactors = F)

		if(IDType == "Symbol"){

			rownames(Ann) <- Ann$Symbol

			result.ann <- cbind(result, Ann[rownames(result), ])

		}else if(IDType == "Ensembl"){

                        rownames(Ann) <- Ann$Ensembl

                        result.ann <- cbind(result, Ann[rownames(result), ])

		}

	        result.ann <- cbind(result.ann, Expr[rownames(result.ann), ], Expr.norm[rownames(result.ann), ])

                write.table(result.ann, file = paste("DEG", "diff.ann.xls", sep="_"), sep = "\t", quote = F, col.names = NA)

	}

}

if(AdjustBy != "No"){

        Factor <- Factor[, c(Focus, AdjustBy), drop = FALSE]

}else{

        Factor <- Factor[, Focus, drop = FALSE]

}

Factor <- Factor[order(Factor[[Focus]]), , drop = FALSE]

if(Method == "voom"){

SumFac <- colSums(Expr)

Expr <- t(apply(Expr, 1, function(x) x/SumFac*1000000))

ExprTmp <- log2(Expr[rownames(result)[result$limma.p <= ThresP], rownames(Factor)] + 1)

if(IDType == "Ensembl"){

	RowNames <- result.ann[match(rownames(ExprTmp), result.ann$Ensembl), "Symbol"]

	ExprTmp <- ExprTmp[!is.na(RowNames), ]

	rownames(ExprTmp) <- RowNames[!is.na(RowNames)]

}

write.table(ExprTmp, file = "Expr.sig.tpm.log2.xls", sep = "\t", quote = F, col.names = NA)

ExprTmp <- t(apply(ExprTmp, 1, function(x) (x-mean(x))/sd(x)))

write.table(Factor, file = "FactorTmp.xls", sep = "\t", quote = F, col.names = NA)

}else{

        ExprTmp <- Expr[rownames(result)[result$limma.p <= ThresP], rownames(Factor)]

	ExprTmp <- t(apply(ExprTmp, 1, function(x) (x-mean(x))/sd(x)))

}

pdf(paste0("DEG", "_DEG.pdf"))

        pheatmap(ExprTmp, cluster_rows = TRUE, cluster_cols = F, annotation_col = Factor, show_rownames = FALSE, show_colnames = FALSE)

dev.off()

Tab <- result.ann[, 1:8]

Tab <- Tab[!duplicated(Tab$Symbol) & !is.na(Tab$Symbol), ]

rownames(Tab) <- Tab$Symbol

Tab$Group <- "No"

Tab$Group[Tab$Log2Rat > 0 & Tab$limma.p <= ThresP] <- "Up"

Tab$Group[Tab$Log2Rat < 0 & Tab$limma.p <= ThresP] <- "Down"


TabSig <- Tab[Tab$limma.p <= ThresP, ]

TabSig <- TabSig[c(head(order(TabSig$Log2Rat, decreasing = T), 20), head(order(TabSig$Log2Rat, decreasing = F), 20)), ]

pp <- ggplot(Tab, aes(y=-log10(limma.p), x=Log2Rat, color = Group)) +

        geom_point(alpha = 1) +

	scale_color_manual(values = Colors[c("Blue", "Yellow", "Red"), ]) +

        geom_text_repel(data = TabSig, aes(label = Symbol), size = 3) +

        xlab("log2(Fold change)") +

        ylab("-log10(P value)") +

        geom_hline(yintercept=-log10(ThresP), linetype="dashed", size = 0.5) +

        theme_minimal() +

        theme(legend.position = "none")

ggsave(paste0("DEG", "_Volcano.pdf"), width = 3, height = 3)


Plots <- list()

Factor[[Focus]] <- as.factor(Factor[[Focus]])

for(gg in result.ann$Symbol[1:9]){

        Factor$Expr <- ExprTmp[gg, rownames(Factor)]

        pp <- ggplot(Factor, aes_string(x = Focus, y = "Expr")) +

                geom_boxplot(outlier.shape = NA, fill = brewer.pal(8, "Set2")[3]) +

                geom_jitter(shape = 16, position = position_jitter(0.2)) +

                ylab(gg) +

                theme_minimal()

        Plots[[gg]] <- pp

}

ggsave(paste0("DEG", "_TopGene.pdf"), arrangeGrob(grobs = Plots, ncol = 3))

