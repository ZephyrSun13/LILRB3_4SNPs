
library(ggplot2)

Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)

AltMat <- readRDS(file = "AltMat.rds")

AltMat[AltMat<=5] <- 0

RefMat <- readRDS(file = "RefMat.rds")

RefMat[RefMat<=5] <- 0

SNPMat <- readRDS(file = "SNPMat.rds")

#colnames(SNPMat) <- paste0(colnames(SNPMat), "_BS")


RefMat.ori <- RefMat[rownames(RefMat) %in% rownames(AltMat), ]

RefMat <- matrix(0, nrow = nrow(AltMat), ncol = ncol(AltMat), dimnames = list(rownames(AltMat), colnames(AltMat)))

RefMat[rownames(RefMat.ori), colnames(RefMat.ori)] <- RefMat.ori


AltRate <- AltMat/(AltMat + RefMat) * 100

saveRDS(AltRate, file = "AltRate.rds")


#AltRate <- readRDS(file = "AltRate.rds")

Tab <- data.frame(Rates = AltRate[AltRate > 0 & !is.na(AltRate)])

pp <- ggplot(Tab, aes(x = Rates)) +

	geom_density(fill = Colors["SteelDark", ]) +

	theme_classic()

ggsave("Rate_Density.pdf", pp, height = 4, width = 6)


AltRate <- AltRate[rownames(AltRate) %in% rownames(SNPMat), colnames(SNPMat)]

SNPMat <- SNPMat[rownames(AltRate), colnames(AltRate)]


Tab <- data.frame(Rates = AltRate[SNPMat == 1])

pp <- ggplot(Tab, aes(x = Rates)) +

        geom_density(fill = Colors["SteelDark", ]) +

	xlab("Alt fraction") +

        theme_classic()

ggsave("Rate_Density_heta.pdf", pp, height = 2, width = 3)

