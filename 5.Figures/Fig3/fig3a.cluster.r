
## R/4.0.3

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


AltRate <- readRDS(file = "AltRate.rds")

Colors2 <- read.delim(file = "ColorLibrary", header = F, row.names = 1)

DCGL <- read.delim(file = "DCGL.over.xls.SigSNPTab.sele.xls.ann.xls", row.names = 5)


Clinical <- read.delim(file = "GOCAR_Clinical.xls", row.names = 1, as.is = T)

rownames(Clinical) <- paste0("P", gsub("Rb", "", rownames(Clinical)))


AnnCol <- Clinical[colnames(AltRate), c("Gender", "Age", "death_censored_graft_loss_event", "Genetic_Rec_Race")]

AnnCol$Genetic_Rec_Race <- factor(AnnCol$Genetic_Rec_Race)

levels(AnnCol$Genetic_Rec_Race) <- c("AA", "European", "Hispanic")

AnnCol$Gender <- factor(AnnCol$Gender)

levels(AnnCol$Gender) <- c("F", "M")

AnnCol$death_censored_graft_loss_event <- factor(AnnCol$death_censored_graft_loss_event)

levels(AnnCol$death_censored_graft_loss_event) <- c("No", "Yes")

AnnRow <- DCGL[c(grep("LILR", DCGL$Gene)), c("Gene", "Func", "Rate_Pr...z..")]

colnames(AnnRow) <- c("Gene", "Syn", "Pval")

Tab <- AltRate[rownames(DCGL)[c(grep("LILR", DCGL$Gene))], ]

column_ha = HeatmapAnnotation(DCGL = AnnCol$death_censored_graft_loss_event, Race = AnnCol$Genetic_Rec_Race, Gender = AnnCol$Gender, Age = AnnCol$Age, col = list(Race = c("AA" = Colors2["Red", ], "Hispanic" = Colors2["SteelDark", ], "European" = Colors2["Blue", ]), DCGL = c("Yes" = Colors2["HeatPurple8", ], "No" = Colors2["HeatWt", ])))

row_ha = rowAnnotation(Gene = AnnRow$Gene, Syn = AnnRow$Syn, Pval = AnnRow$Pval)

#Cols = colorRamp2(c(min(Tab, na.rm = T), min(Tab[Tab>0], na.rm = T), max(Tab, na.rm = T)), brewer.pal(9, "YlOrRd")[c(3, 5, 9)])

Cols = colorRamp2(c(0, 50, 100), brewer.pal(9, "YlOrRd")[c(3, 5, 9)])

pdf("DCGLSNP.cluster.lilr.only.pdf", width = 10)

        Heatmap(Tab, col = Cols, na_col = brewer.pal(8, "Set2")[8], row_names_side = "left", row_dend_side = "left", column_names_side = "bottom", column_dend_side = "top", top_annotation = column_ha, cluster_rows = T, cluster_columns = T, left_annotation = row_ha, show_column_names = F, show_heatmap_legend = T, row_names_gp = gpar(fontsize = 11))

dev.off()

