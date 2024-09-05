
Args <- commandArgs(trailingOnly = T)

library(survival)
library(survminer)
library(dplyr)
library(lubridate)

EgfrN <- as.numeric(Args[1]) ## 15 or 30
GenoS <- Args[2]
Diseases <- unlist(strsplit(Args[3], split = ";"))
RaceS <- Args[4]
LabelS <- unlist(strsplit(Args[5], split = ";"))
LowerN <- as.numeric(Args[6]) ## 90
UpperN <- as.numeric(Args[7]) ## 365
AgeLimitN <- as.numeric(Args[8])

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", row.names = 1, as.is = T, header = F)


DiaGTab <- readRDS(file = "DiaGTab.egfr.rds")

DiaGTab <- DiaGTab %>% group_by(Subject_id)

DiaGTab <- DiaGTab %>% arrange(mdy(proc_code), .by_group = T)


Tab <- read.delim(file = "TabMerge.AAHis.xls.geno.xls", row.names = 1, as.is = T)

if(RaceS == "AAHis"){

        Tab <- Tab

}else if(RaceS == "AA"){

        Tab <- Tab[Tab$self_reported_race %in% c("African American", "AFRICAN AMERICAN (BLACK)"), ]

}else if(RaceS == "His"){

        Tab <- Tab[Tab$self_reported_race == "Hispanic", ]

}


CodeTab <- readRDS(file = "CodeTab.all.rds")

CodeTab <- CodeTab[rownames(CodeTab) %in% rownames(Tab), ]

Tab <- Tab[rownames(Tab) %in% rownames(CodeTab), ]


CodeDis <- read.table(file = "Diseases.code", row.names = 1, as.is = T)

CodeDis.lst <- lapply(CodeDis[[1]], function(x) unlist(strsplit(x, split = ";")))

names(CodeDis.lst) <- rownames(CodeDis)

for(nn in names(CodeDis.lst)){

        Tab[[nn]] <- as.integer(rowSums(CodeTab[rownames(Tab), CodeDis.lst[[nn]], drop = F]) > 0)

}

DiaGTab <- DiaGTab[DiaGTab$Subject_id %in% rownames(Tab), ]


Ques <- read.delim(file = "QuesTab.xls", row.names = 1, as.is = T)

Ques <- Ques[Ques$platform %in% DiaGTab$Subject_id, ]

Ques <- Ques[!duplicated(Ques$platform), ]

rownames(Ques) <- Ques$platform


Result <- matrix(0, ncol = 10, nrow = 1, dimnames = list("Header", c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N", "Event", "Disease")))

for(dd in Diseases){

        TabTmp <- DiaGTab[DiaGTab$Subject_id %in% rownames(Tab)[Tab[[dd]] == 1], ]

	print(dd)

	quantile(table(TabTmp$Subject_id), probs = seq(0, 1, by = 0.1))

	PIDs <- unique(TabTmp$Subject_id)

	RecCount <- table(TabTmp$Subject_id)

	print(quantile(RecCount, probs = seq(0, 1, by = 0.1)))

	RecDur <- tapply(TabTmp$proc_code, TabTmp$Subject_id, function(x) as.Date(x[length(x)], format = "%m/%d/%Y") - as.Date(x[1], format = "%m/%d/%Y"))

	print(quantile(RecDur, probs = seq(0, 1, by = 0.1)))

        SuvTab <- matrix(NA, nrow = length(PIDs), ncol = 2, dimnames = list(PIDs, c("Date", "Event")))

	for(pp in PIDs){

		TabPP <- TabTmp[TabTmp$Subject_id == pp, ]

		if(nrow(TabPP) <= 4) next

		if(as.Date(TabPP$proc_code[nrow(TabPP)], format = "%m/%d/%Y") - as.Date(TabPP$proc_code[1], format = "%m/%d/%Y") <= 10) next

		Lows <- as.Date(TabPP$proc_code[which(as.numeric(TabPP$quantity) < EgfrN)], format = "%m/%d/%Y")

		if(length(Lows) > 0){

		for(ll in 1:length(Lows)){

			if(any((Lows - Lows[ll]) >= LowerN & (Lows - Lows[ll]) <= UpperN)){

				SuvTab[pp, "Date"] <- as.character(Lows[ll])

				SuvTab[pp, "Event"] <- 1

				break

			}

		}

		}

		if(is.na(SuvTab[pp, "Date"])){

			SuvTab[pp, "Date"] <- as.character(as.Date(TabPP$proc_code[nrow(TabPP)], format = "%m/%d/%Y"))

                        SuvTab[pp, "Event"] <- 0

		}

	}


        SuvTab <- data.frame(SuvTab)

        #SuvTab$Time <- as.numeric(gsub("A_", "", SuvTab$Time))

        SuvTab$Event <- as.numeric(SuvTab$Event)

        SuvTab$Geno <- Tab[rownames(SuvTab), GenoS]

	SuvTab$YOB <- paste(Ques[rownames(SuvTab), "YEAR_OF_BIRTH"], "01", "01", sep = "-")

	SuvTab$Time <- as.Date(SuvTab$Date, format = "%Y-%m-%d") - as.Date(SuvTab$YOB, format = "%Y-%m-%d")

	SuvTab$Time <- SuvTab$Time/365

        write.table(SuvTab, file = paste0(GenoS, "_", EgfrN, "_", dd, "_", RaceS, "_", UpperN, ".suv.xls"), sep = "\t", quote = F, col.names = NA)

        SuvTab$Geno[which(SuvTab$Geno == "./.")] <- NA

	SuvTab$Event[SuvTab$Time > AgeLimitN] <- 0

        SuvTab$Time[SuvTab$Time > AgeLimitN] <- AgeLimitN

	#SuvTab <- SuvTab[SuvTab$Event == 1, ]


        fit <- coxph(Surv(Time, Event) ~ Geno, data = SuvTab)

        Result <- rbind(Result, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n, summary(fit)$nevent, dd))

        fit <- survfit(Surv(Time, Event) ~ Geno, data = SuvTab)

        pp <- ggsurvplot(fit,
                  pval = TRUE, conf.int = F,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "hv", # Specify median survival
                  ggtheme = theme_minimal(), # Change ggplot2 theme
                  palette = Colors[c("Blue", "Red", "Yellow", "HeatPurple7"), ],
                  xlab = "Age",
                  legend.title = "",
		  legend = "none",
                  legend.labs = LabelS,
                  fontsize = 6,
                  font.x = 15,
                  font.y = 15,
                  font.legend = 15,
                  font.tickslab = 15,
        )

        ggsave(file = paste0(GenoS, "_", EgfrN, "_", dd, "_", RaceS, "_", UpperN, ".suv.pdf"), print(pp), width = 5.5, height = 6)

}

write.table(Result, file = paste0(GenoS, "_", EgfrN, "_", RaceS, "_", UpperN, ".result.xls"), sep = "\t", quote = F, col.names = NA)

