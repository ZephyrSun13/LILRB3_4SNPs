
Args <- commandArgs(trailingOnly = T)

library(survival)
library(survminer)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(pROC)

SNP <- Args[1]

Colors2 <- read.delim(file = "ColorLibrary", header = F, row.names = 1)


AltRate <- readRDS(file = "AltRate.rds")

colnames(AltRate) <- gsub("_BS", "", colnames(AltRate))


Map <- read.delim(file = "AltSNPGene.map.uniq", row.names = 5, as.is = T, header = F)


Clinical <- read.delim(file = "GOCAR_Clinical.xls", row.names = 1, as.is = T)

rownames(Clinical) <- paste0("P", gsub("Rb", "", rownames(Clinical)))

Clinical$graft_loss_time_update <- Clinical$graft_loss_time_update/30


	Expr <- unlist(AltRate[SNP, ])

	quantile(Expr, na.rm = T)

	Clinical_tt <- Clinical[names(Expr), ]
	
	Clinical_tt$Expr <- Expr[rownames(Clinical_tt)]

	Clinical_tt <- Clinical_tt[Clinical_tt$Batch == "Pair", ]

	Tab <- Clinical_tt

	Tab$Genetic_Rec_Race <- factor(Tab$Genetic_Rec_Race)

	levels(Tab$Genetic_Rec_Race) <- c("AA", "European", "Hispanic")	

	pp <- ggplot(Tab[!is.na(Tab$Genetic_Rec_Race), ], aes(x = Genetic_Rec_Race, y = Expr, fill = Genetic_Rec_Race)) +

		geom_boxplot() +

                geom_jitter(shape = 16, position = position_jitter(0.2)) +

		scale_fill_manual(values = Colors2[c("Red", "Blue", "SteelDark"), ]) +

		xlab("") + ylab("AEF of rs549267286") +

                theme_classic() +

		theme(legend.position = "null", axis.text.x = element_text(angle = 45, hjust = 1))

	ggsave("RaceExprBox.pdf", pp, width = 4, height = 4)


	Tab <- Tab[order(Tab$Expr, decreasing = T), ]

	Tab$ID <- factor(rownames(Tab), levels = rownames(Tab))

       pp <- ggplot(Tab[!is.na(Tab$Genetic_Rec_Race), ], aes(x = ID, y = Expr, fill = Genetic_Rec_Race)) +

                geom_bar(stat = "identity") +

                scale_fill_manual(values = Colors2[c("Red", "Blue", "SteelDark"), ]) +

                xlab("") + ylab("AEF of rs549267286") +

                theme_classic() +

                theme(axis.text.x = element_blank(), legend.title = element_blank())

        ggsave("RaceExprBar.pdf", pp, width = 8, height = 2)


        Plots <- list()

        ColorsTmp <- Colors2[c("Red", "Blue", "SteelDark"), ]

        names(ColorsTmp) <- c("AA", "European", "Hispanic")

        for(rr in c("AA", "Hispanic", "European")){

                pp <- ggplot(Tab[which(Tab$Genetic_Rec_Race == rr), ], aes(x = ID, y = Expr)) +

                        geom_boxplot() +

                	geom_bar(stat = "identity", fill = Colors2["SteelDark", ]) +

                        xlab("") + ylab("AEF of rs549267286") +

			ylim(0, 50) +

                        ggtitle(rr) +

                        theme_classic() +

                        theme(legend.position = "null", axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5))

		Plots[[rr]] <- pp

        }

        ggsave("RaceExprBar.sep.pdf", arrangeGrob(grobs = Plots, ncol = 3), width = 8, height = 2)


	
	if(length(which(Clinical_tt$Expr == 0))>=nrow(Clinical_tt)/2){

	Clinical_tt$ExprG <- "Medium"

        Clinical_tt$ExprG[which(Clinical_tt$Expr >= mean(Clinical_tt$Expr[which(Clinical_tt$Expr > 0)]))] <- "High"

	Clinical_tt$ExprG[Clinical_tt$Expr == 0] <- "Low"

	Clinical_tt$ExprG <- factor(Clinical_tt$ExprG, levels = c("Low", "Medium", "High"))

	}else{

                vTert = quantile(Clinical_tt$Expr, c(0:3/3))

                Clinical_tt$ExprG= with(Clinical_tt, cut(Expr, vTert, include.lowest = T, labels = c("Low", "Medium", "High")))

		levels(Clinical_tt$ExprG) <- c("Low", "Medium", "High")

	}

                fit <- survfit(Surv(graft_loss_time_update, death_censored_graft_loss_event) ~ ExprG, data = Clinical_tt)

                Colors <- brewer.pal(9, "Set1")

		pp <- ggsurvplot(fit, palette = c(Colors[2], Colors[3], Colors[1]), ggtheme = theme_minimal(), xlab = "Months", conf.int = F, pval = TRUE, pval.size = 6, legend.labs=c("Low", "Medium", "High"), font.x = 18, font.y = 18, font.legend = 18, font.tickslab = 18, legend.title = "", risk.table = T, font.title = c(16, "bold", "darkblue"))

                ggsave(paste0(SNP, "_", "Rate", "_Survival_Univariate_AS.group.pair.pdf"), print(pp), width = 5, height = 5)

