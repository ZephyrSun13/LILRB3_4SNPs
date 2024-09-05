
Args <- commandArgs(trailingOnly = T)

WESF <- Args[1] ## WES.rds
Pop <- Args[2] ## AA
CodeTabF <- Args[3]	## ../CodeTab.all.rds.filt.rds
ModelS <- Args[4] ## Add
SNPS <- Args[5] ## 19_54217126_C_T 
AdjS <- Args[6] ## + gender + Race
DemoS <- Args[7] ## /sc/arion/projects/zhangw09a/Data/Zeguo_Sun/project/6.Transplant/8.LILR/14.BioBank/TabMerge.xls 

if(ModelS == "Add"){

	GenoMatch <- c(0, 1, 2)

	names(GenoMatch) <- c("0/0", "0/1", "1/1")

}else if(ModelS == "Dom"){

        GenoMatch <- c(0, 1, 1)

        names(GenoMatch) <- c("0/0", "0/1", "1/1")

}else if(ModelS == "Rec"){

        GenoMatch <- c(0, 0, 1)

        names(GenoMatch) <- c("0/0", "0/1", "1/1")

}else{

	stop("Wrong model")

}

WES <- readRDS(file = WESF)

CodeTab <- readRDS(file = CodeTabF)

TabMerge <- read.delim(file = DemoS, row.names = 1, as.is = T)

TabMerge <- TabMerge[which(rownames(TabMerge) %in% rownames(CodeTab)), ]

TabMerge$Race <- "Other"

TabMerge$Race[(TabMerge$self_reported_race %in% c("African American", "AFRICAN AMERICAN (BLACK)"))] <- "AA"

TabMerge$Race[(TabMerge$self_reported_race %in% c("Hispanic"))] <- "Hispanic"

TabMerge$Race[(TabMerge$self_reported_race %in% c("European American", "CAUCASIAN (WHITE)"))] <- "European"

TabMerge$Race[(TabMerge$self_reported_race %in% c("Asian", "ASIAN"))] <- "Asian"

TabMerge$Race <- factor(TabMerge$Race, levels = c("Other", "AA", "Hispanic", "European", "Asian"))


if(Pop == "AA"){

        Tab <- TabMerge[which(TabMerge$self_reported_race %in% c("African American", "AFRICAN AMERICAN (BLACK)")), ]

}else if(Pop == "His"){

	Tab <- TabMerge[which(TabMerge$self_reported_race %in% c("Hispanic")), ]

}else if(Pop == "AAHis" ){

        Tab <- TabMerge[which(TabMerge$self_reported_race %in% c("African American", "AFRICAN AMERICAN (BLACK)", "Hispanic")), ]

}else if(Pop == "All"){

        Tab <- TabMerge

}else{

	stop("No such population")

}


WES <- WES[, rownames(Tab)]

dim(WES)

#Freq <- apply(WES, 1, function(x) length(which(!(x %in% c("0/0", "./.")))))
#
#Freq <- data.frame(QCCount = Freq, QCFreq = Freq/ncol(WES))
#
#Freq$MAF <- apply(WES, 1, function(x) {(length(which(x == "1/1"))*2 + length(which(x == "0/1")))/(length(which(x != "./.")) * 2)})
#
#WES <- WES[which(Freq$MAF >= 0.01), ]

Outcomes <- colnames(CodeTab)

ResultM <- matrix(NA, nrow = length(Outcomes), ncol = 9, dimnames = list(Outcomes, c("NA", "N0", "N1", "N2", "CI1", "CI2", "OR", "Pval", "N")))

#ResultOR <- matrix(NA, nrow = nrow(WES) , ncol = length(Outcomes), dimnames = list(rownames(WES), Outcomes))

#for(ss in rownames(WES)){

	#if(grep(ss, rownames(WES)) %% 50 == 0){print(grep(ss, rownames(WES)))}

        TabTmp <- Tab

	if(SNPS %in% rownames(WES)){

        	TabTmp$Genotype <- unlist(WES[SNPS, rownames(TabTmp)])

	}else{

		TabTmp$Genotype <- TabTmp[[SNPS]]

	}

	TabTmp$Genotype2 <- factor(TabTmp$Genotype, levels = c("./.", "0/0", "0/1", "1/1"))

	AlleFreq <- table(factor(TabTmp$Genotype, levels = c("./.", "0/0", "0/1", "1/1")))

	#table(factor(TabTmp$Genotype, levels = c("./.", "0/0", "0/1", "1/1")))

	#TabTmp <- TabTmp[TabTmp$Genotype != "./.", ]

        TabTmp$Genotype <- GenoMatch[TabTmp$Genotype]

	ResultM <- matrix(NA, nrow = length(Outcomes), ncol = 9, dimnames = list(Outcomes, c(paste(c("NA", "N0", "N1", "N2"), AlleFreq, sep = "_"), "CI1", "CI2", "OR", "Pval", "N")))

	for(oo in Outcomes){

		TabTmp$Phenotype <- CodeTab[rownames(TabTmp), oo]

	        tryCatch({

			if(AdjS == "No"){
        
                		fit <- glm(as.formula(paste0("Phenotype", " ~ ", "Genotype")), data = TabTmp, family = "binomial")

			}else{

				fit <- glm(as.formula(paste0("Phenotype", " ~ ", "Genotype", AdjS)), data = TabTmp, family = "binomial")

			}
                	
                	#ResultP[ss, oo] <- summary(fit)$coefficients["Genotype", "Pr(>|z|)"]
                	
                	#ResultOR[ss, oo] <- summary(fit)$coefficients["Genotype", "Estimate"]
                	
                	#fit <- glm(as.formula(paste0(oo, " ~ ", ss, " + Race")), data = TabTmp, family = "binomial")

                	ResultM[oo, ] <- c(tapply(TabTmp$Phenotype, TabTmp$Genotype2, sum), exp(suppressMessages(confint(fit)["Genotype", ])), exp(summary(fit)$coefficients["Genotype", "Estimate"]), summary(fit)$coefficients["Genotype", "Pr(>|z|)"], summary(fit)$df.null + 1)
                
        	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


	}

write.table(ResultM, file = paste0("CodeTab_", ModelS, "_", SNPS, "_", Pop, ".xls"), sep = "\t", quote = F, col.names = NA)

CodeMap <- readRDS(file = "/sc/arion/projects/zhangw09a/Data/Zeguo_Sun/project/6.Transplant/8.LILR/14.BioBank/CodeMap.rds")

write.table(cbind(CodeMap[rownames(ResultM)], ResultM), file = paste0("CodeTab_", ModelS, "_", SNPS, "_", Pop, ".ann.xls"), sep = "\t", quote = F, col.names = NA)

