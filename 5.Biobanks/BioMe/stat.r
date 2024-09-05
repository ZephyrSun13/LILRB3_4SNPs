
WES <- read.delim(file = "LILR.WES.vcf.tab", row.names = 1, as.is = T, check.names = F)

rownames(WES) <- gsub(" ", "_", rownames(WES))

WESGQ <- read.delim(file = "LILR.WES.vcf.GQ.tab", row.names = 1, as.is = T, check.names = F)

rownames(WESGQ) <- gsub(" ", "_", rownames(WESGQ))


for(ss in rownames(WES)){

	WES[ss, which(as.numeric(unlist(WESGQ[ss, colnames(WES)])) < 13)] <- "./."

}

saveRDS(WES, file = "WES.qc.rds")


GSA <- read.delim(file = "LILR.GSA.vcf.tab", row.names = 1, as.is = T, check.names = F)

rownames(GSA) <- gsub(" ", "_", rownames(GSA))

#colnames(GSA) <- sapply(colnames(GSA), function(x) paste(unlist(strsplit(x, split = "_"))[1:3], collapse = "_"))

Map <- read.table(file = "RGNIDupdate.txt", as.is = T, sep = " ")

Map <- Map[!is.na(Map[[1]]), ]

Map <- Map[!duplicated(Map[[1]]), ]

#Map <- Map[Map[[3]] %in% rownames(CodeTab), ]

rownames(Map) <- paste(Map[[1]], Map[[2]], sep = "_")

GSA <- GSA[, colnames(GSA) %in% rownames(Map)]

colnames(GSA) <- Map[colnames(GSA), 3]

saveRDS(GSA, file = "GSA.rds")


table(unlist(WES["19_54217126_C_T", ]))

table(data.frame(t(WES[c("19_54217126_C_T", "19_54217180_A_G"), ])))

Freq <- apply(WES, 1, function(x) length(which(!(x %in% c("0/0", "./.")))))


table(unlist(GSA["19_54721006_G_A", ]))

table(unlist(GSA["19_54721007_T_G", ]))

table(unlist(GSA["19_54721008_T_G", ]))

table(data.frame(t(GSA[c("19_54721007_T_G", "19_54721049_A_G"), ])))

Inter <- intersect(colnames(GSA), colnames(WES))

table(data.frame(WES = WES["19_54217126_C_T", Inter]), GSA = GSA["19_54721007_T_G", Inter])

#table(data.frame(LILR995 = unlist(WES["19_54217126_C_T", ]), LILR143 = unlist(WES["19_54217143_C_T", ])))

#  ./.   0/0   0/1   1/1 
#  491 29407   877    38 

#WES["19_54217126_C_T", which(as.numeric(unlist(WESGQ["19_54217126_C_T", ])) < 13)] <- "./."

#table(unlist(WES["19_54217126_C_T", ]))

#  ./.   0/0   0/1   1/1 
#  620 29296   859    38 


table(unlist(WES["22_36265988_T_G", ]))

table(unlist(WES["22_36265995_AATAATT_A", ]))

#table(t(WES[c("19_54217126_C_T", "22_36265988_T_G"), ]))

table(data.frame(LILR = unlist(WES["19_54217126_C_T", ]), G1 = unlist(WES["22_36265988_T_G", ])))

table(data.frame(LILR = unlist(WES["19_54217126_C_T", ]), G2 = unlist(WES["22_36265995_AATAATT_A", ])))

table(data.frame(G1 = unlist(WES["22_36265988_T_G", ]), G2 = unlist(WES["22_36265995_AATAATT_A", ])))


GSA <- read.delim(file = "LILR.GSA.vcf.tab", row.names = 1, as.is = T, check.names = F)

rownames(GSA) <- gsub(" ", "_", rownames(GSA))

table(unlist(GSA["19_54721007_T_G", ]))

#  ./.   0/0   0/1 
#  736 30868    97 

table(unlist(GSA["22_36662034_T_G", ]))

table(unlist(GSA["22_36662041_AATAATT_A", ]))

table(data.frame(LILR = unlist(GSA["19_54721007_T_G", ]), G1 = unlist(GSA["22_36662034_T_G", ])))

table(data.frame(LILR = unlist(GSA["19_54721007_T_G", ]), G2 = unlist(GSA["22_36662041_AATAATT_A", ])))


#Demo <- read.table(file = "/sc/private/regen/data/BioMe/BRSSPD/Demographics.txt", sep = "\\|", as.is = T, head = T)

Demo <- readLines("/sc/private/regen/data/BioMe/BRSSPD/Demographics.txt")

Demo <- lapply(Demo, function(x) unlist(strsplit(x, split = "\\|")))

#Demo <- Reduce(cbind, Demo)

DemoTab <- matrix(NA, nrow = length(Demo) - 1, ncol = length(Demo[[1]]))

colnames(DemoTab) <- Demo[[1]]

for(nn in 2:length(Demo)){

	if(length(Demo[[nn]]) == 9){

		DemoTab[nn -1, ] <- c(Demo[[nn]][1:3], NA, Demo[[nn]][4:9])

	}else{

		DemoTab[nn - 1, ] <- Demo[[nn]]

	}

}

DemoTab <- data.frame(DemoTab)

rownames(DemoTab) <- DemoTab$subject_id

length(intersect(colnames(WES), rownames(DemoTab)))

#Diss <- setdiff(colnames(WES), rownames(DemoTab))

#Cols <- which(colnames(WES) %in% Diss)

#colnames(WES)[Cols] <- sapply(colnames(WES)[Cols], function(x) paste(unlist(strsplit(x, split = "_"))[1:2], collapse = "_"))

OverL <- intersect(colnames(WES), rownames(DemoTab))

DemoTab <- DemoTab[OverL, ]

DemoTab$LILR995 <- unlist(WES["19_54217126_C_T", rownames(DemoTab)])

DemoTab$G1 <- unlist(WES["22_36265988_T_G", rownames(DemoTab)])

DemoTab$G2 <- unlist((unlist(WES["22_36265995_AATAATT_A", rownames(DemoTab)])))

DemoTab$APOL1 <- NA

DemoTab$APOL1[DemoTab$G1 == "0/1" & DemoTab$G2 == "0/0"] <- "WT/G1"

DemoTab$APOL1[DemoTab$G1 == "0/0" & DemoTab$G2 == "0/1"] <- "WT/G2"

DemoTab$APOL1[DemoTab$G1 == "1/1"] <- "G1/G1"

DemoTab$APOL1[DemoTab$G2 == "1/1"] <- "G2/G2"

DemoTab$APOL1[DemoTab$G1 == "0/1" & DemoTab$G2 == "0/1"] <- "G1/G2"

DemoTab$APOL1[DemoTab$G1 == "0/0" & DemoTab$G2 == "0/0"] <- "WT/WT"


Pheno <- read.table(file = "/sc/private/regen/data/BioMe/BRSSPD/Phenotype_CDA.txt", sep = "|", as.is = T, head = T)

rownames(Pheno) <- Pheno$subject_id


TabMerge <- merge(DemoTab, Pheno, by = "subject_id")

#write.table(TabMerge, file = "TabMerge.xls", sep = "\t", quote = F, row.names = F)


table(TabMerge[, c("self_reported_race", "LILR995")])

table(TabMerge[, c("self_reported_race", "G1")])

table(TabMerge[, c("self_reported_race", "G2")])


#TabAA$H_CKD_CASE_CONTROL[TabAA$H_CKD_CASE_CONTROL == -9] <- NA

#TabAA$H_CKD_CASE_CONTROL[TabAA$H_CKD_CASE_CONTROL == 2] <- 0

#fit <- glm(as.formula(paste0("H_CKD_CASE_CONTROL ~ LILR995")), data = TabAA, family = "binomial")

TabMerge$APOL1Risk <- NA

TabMerge$APOL1Risk[TabMerge$APOL1 == "WT/WT"] <- "No"

TabMerge$APOL1Risk[TabMerge$APOL1 %in% c("WT/G1", "WT/G2", "G1/G1", "G2/G2", "G1/G2")] <- "Yes"

TabMerge$LILRRisk <- NA

TabMerge$LILRRisk[TabMerge$LILR995 == "0/0"] <- "No"

TabMerge$LILRRisk[TabMerge$LILR995 %in% c("0/1", "1/1")] <- "Yes"

write.table(TabMerge, file = "TabMerge.xls", sep = "\t", quote = F, row.names = F)


TabAA <- TabMerge[which(TabMerge$self_reported_race %in% c("African American", "AFRICAN AMERICAN (BLACK)")), ]

SNPs <- c("LILR995", "G1", "G2", "APOL1Risk", "LILRRisk")

Outcomes <- c(grep("CONTROL", colnames(TabAA), value = T), grep("FLAG", colnames(TabAA), value = T))

Result <- matrix(NA, nrow =1 , ncol = 7)

for(ss in SNPs){

	TabTmp <- TabAA[TabAA[[ss]] != "./.", ]

for(oo in Outcomes){

	#TabTmp[[oo]][TabTmp[[oo]]== -9] <- NA

	TabTmp[[oo]][TabTmp[[oo]]== -9] <- 0

	TabTmp[[oo]][TabTmp[[oo]]== 2] <- 0

	tryCatch({

	fit <- glm(as.formula(paste0(oo, " ~ ", ss)), data = TabTmp, family = "binomial")

	Result <- rbind(Result, cbind(confint(fit), summary(fit)$coefficients, oo))

	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

}

write.table(Result, file = "Result.xls", sep = "\t", quote = F, col.names = NA)

table(TabAA[, c("LILR995", "AKI_CASE_CONTROL")])

table(TabAA[, c("G1", "DH_CKD_CASE_CONTROL")])

table(TabAA[, c("G2", "DH_CKD_CASE_CONTROL")])

table(TabAA[, c("APOL1", "DH_CKD_CASE_CONTROL")])


Result2 <- matrix(NA, nrow =1 , ncol = 7)

for(ss in SNPs){

        #TabTmp <- TabAA[TabAA[[ss]] != "./.", ]

	TabTmp <- TabAA[, Outcomes]

	#TabTmp[TabTmp == -9] <- NA

	TabTmp[TabTmp == 2] <- 0

	Freq <- apply(TabTmp, 2, function(x) Reduce(cbind, tapply(factor(x), TabAA[[ss]], function(y) table(y))))

	write.table(Freq[["AKI_CASE_CONTROL"]], file = "AKI_CASE_CONTROL.xls", sep = "\t", quote = F, col.names = NA)

for(oo in Outcomes){


        tryCatch({

        fit <- glm(as.formula(paste0(oo, " ~ ", ss)), data = TabTmp, family = "binomial")

        Result2 <- rbind(Result2, cbind(confint(fit), summary(fit)$coefficients, oo))

        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

}

