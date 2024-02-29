
TabsLis <- readRDS(file = "TabsLis.rds")

sort(unlist(lapply(TabsLis, nrow)))

sort(unlist(lapply(TabsLis, function(x) length(which(x$refAllele == x$altAllele)))))

TabsLis <- lapply(TabsLis, function(x) x[x$refAllele != x$altAllele, ])


SNPs <- unique(unlist(lapply(TabsLis, rownames)))

length(SNPs)

#RefMat <- Reduce(cbind, lapply(TabsLis, function(x) x[SNPs, "refCount"]))

RefMat <- matrix(NA, nrow = length(SNPs), ncol = length(TabsLis), dimnames = list(SNPs, names(TabsLis)))

for(nn in names(TabsLis)){

	#print(nn)

	RefMat[match(rownames(TabsLis[[nn]]), SNPs), nn] <- TabsLis[[nn]][["refCount"]]

}

RefMat <- RefMat[apply(RefMat, 1, function(x) length(which(x>=10))) >= 5, ]

saveRDS(RefMat, file = "RefMat.rds")



AltMat <- matrix(NA, nrow = length(SNPs), ncol = length(TabsLis), dimnames = list(SNPs, names(TabsLis)))

for(nn in names(TabsLis)){

        #print(nn)

        AltMat[match(rownames(TabsLis[[nn]]), SNPs), nn] <- TabsLis[[nn]][["altCount"]]

}

AltMat <- AltMat[apply(AltMat, 1, function(x) length(which(x>=10))) >= 5, ]

saveRDS(AltMat, file = "AltMat.rds")


AltMat[AltMat<=5] <- 0


