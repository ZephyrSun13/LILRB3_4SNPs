
Args <- commandArgs(trailingOnly = T)


library(survival)
library(survminer)
library(dplyr)
library(lubridate)


EgfrN <- as.numeric(Args[1])
GenoS <- Args[2]
DisS <- Args[3]
RankN <- as.numeric(Args[4])
LowerN <- as.numeric(Args[5])
UpperN <- as.numeric(Args[6])
LabelS <- unlist(strsplit(Args[7], split = ";"))
AgeLimitN <- as.numeric(Args[8])


Colors <- read.delim(file = "ColorLibrary", row.names = 1, as.is = T, header = F)


DemoTab <- read.delim(file = paste0("AA", DisS, ".demo.xls"), row.names = 1, as.is = T)

DemoTab$APOL1Num <- paste0("N", DemoTab$APOL1Num)

DemoTab$Combine2 <- factor(DemoTab$Combine2, levels = c("No", "APOL1", "LILR", "Combined"))

DemoTab$YOB <- sapply(DemoTab$date_of_birth, function(x) unlist(strsplit(x, split = " "))[1])


DiaGTab <- read.delim(file = paste0("AA", DisS, ".measurement.xls"))

DiaGTab$PID <- paste0("P", DiaGTab$person_id)

DiaGTab$proc_code <- sapply(DiaGTab$measurement_datetime, function(x) unlist(strsplit(x, split = " "))[1])

DiaGTab <- DiaGTab %>% group_by(PID)

DiaGTab <- DiaGTab %>% arrange(ymd(proc_code), .by_group = T)



length(intersect(rownames(DemoTab), DiaGTab$PID))

table(DiaGTab$unit_source_value)

quantile(DiaGTab$value_as_number, na.rm = T)

head(sort(DiaGTab$value_as_number, decreasing = T), 20)

head(sort(DiaGTab$value_as_number, decreasing = F), 20)

DiaGTab <- DiaGTab[which(!is.na(DiaGTab$value_as_number)), ]

DiaGTab <- DiaGTab[which(DiaGTab$value_as_number < 200), ]

head(sort(DiaGTab$value_as_number, decreasing = T), 20)

head(sort(DiaGTab$value_as_number, decreasing = F), 20)

options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

head(DiaGTab)

head(DemoTab)


Result <- matrix(0, ncol = 10, nrow = 1, dimnames = list("Header", c("lower.95", "upper.95", "coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)", "N", "Event", "Disease")))

       TabTmp <- DiaGTab

        head(TabTmp$value_as_number)

        head(TabTmp$proc_code)

        quantile(table(TabTmp$PID), probs = seq(0, 1, by = 0.1))

        PIDs <- unique(TabTmp$PID)

        length(PIDs)

        SuvTab <- matrix(NA, nrow = length(PIDs), ncol = 2, dimnames = list(PIDs, c("Date", "Event")))

       for(pp in PIDs){

                TabPP <- TabTmp[TabTmp$PID == pp, ]

                Lows <- as.Date(TabPP$proc_code[which(as.numeric(TabPP$value_as_number) < EgfrN)], format = "%Y-%m-%d")

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

                        SuvTab[pp, "Date"] <- as.character(as.Date(TabPP$proc_code[nrow(TabPP)], format = "%Y-%m-%d"))

                        SuvTab[pp, "Event"] <- 0

                }

        }

        SuvTab <- data.frame(SuvTab)

        #SuvTab$Time <- as.numeric(gsub("A_", "", SuvTab$Time))

        SuvTab$Event <- as.numeric(SuvTab$Event)

        SuvTab$Geno <- DemoTab[rownames(SuvTab), GenoS]

        SuvTab$YOB <- DemoTab[rownames(SuvTab), "YOB"]

        SuvTab$Time <- as.Date(SuvTab$Date, format = "%Y-%m-%d") - as.Date(SuvTab$YOB, format = "%Y-%m-%d")

        SuvTab$Time <- SuvTab$Time/365

        write.table(SuvTab[order(SuvTab$Time, decreasing = T), ], file = paste0(GenoS, "_", EgfrN, "_", DisS, ".suv.xls"), sep = "\t", quote = F, col.names = NA)


        SuvTab$Geno[which(SuvTab$Geno == "./.")] <- NA

        SuvTab$Event[SuvTab$Time > AgeLimitN] <- 0

        SuvTab$Time[SuvTab$Time > AgeLimitN] <- AgeLimitN

        fit <- coxph(Surv(Time, Event) ~ Geno, data = SuvTab)

        Result <- rbind(Result, cbind(summary(fit)$conf.int[, 3:4, drop = FALSE], summary(fit)$coefficients, summary(fit)$n, summary(fit)$nevent, DisS))

        fit <- survfit(Surv(Time, Event) ~ Geno, data = SuvTab)

        #fit <- survfit(Surv(Time, Event) ~ Geno, data = SuvTab[SuvTab$Time <= 60, ])

	pdf(paste0(GenoS, "_", EgfrN, "_", DisS, "_", UpperN, "_", AgeLimitN, ".suv.pdf"), width = 5.5, height = 6)

        	ggsurvplot(fit,
                  pval = TRUE, conf.int = F,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "hv", # Specify median survival
                  ggtheme = theme_minimal(), # Change ggplot2 theme
                  palette = Colors[c("Blue", "Red", "Yellow", "HeatPurple7"), ],
                  xlab = "Age",
                  legend.title = "",
                  legend.labs = LabelS,
		  legend = "none",
                  fontsize = 6,
                  font.x = 15,
                  font.y = 15,
                  font.legend = 15,
                  font.tickslab = 15,
        )

	dev.off()

        #ggsave(file = paste0(GenoS, "_", EgfrN, "_", DisS, ".suv.pdf"), print(pp), width = 5.5, height = 6)

write.table(Result, file = paste0(GenoS, "_", EgfrN, "_", DisS, "_", UpperN, "_", AgeLimitN, ".result.xls"), sep = "\t", quote = F, col.names = NA)

