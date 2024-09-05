
Args <- commandArgs(trailingOnly = T)

File <- Args[1]

Demo <- read.delim(file = File, row.names = 1, as.is = T, header = T)


Demo$LILRRisk <- factor(Demo$LILRRisk, levels = c("No", "Yes"))

levels(Demo$LILRRisk) <- c("0/0", "0/1")

Demo$APOL1Risk <- factor(Demo$APOL1Risk, levels = c("No", "Yes"))

levels(Demo$APOL1Risk) <- c("0/0", "0/1")


Demo$APOL1[is.na(Demo$APOL1)] <- "./."

Demo$APOL1Geno <- factor(Demo$APOL1)

levels(Demo$APOL1Geno) <- c("./.", "1/1", "1/1", "1/1", "0/1", "0/1", "0/0")

table(Demo[, c("LILR995", "APOL1Geno")])


Demo$APOL1Risk2 <- Demo$APOL1Geno

levels(Demo$APOL1Risk2) <- c("./.", "0/1", "0/0", "0/0")


Demo$Combine <- "./."

Demo$Combine[which((Demo$APOL1Geno %in% c("1/1", "0/1")) & (Demo$LILR995 %in% c("1/1", "0/1")))] <- "1/1"

Demo$Combine[which((Demo$APOL1Geno %in% c("0/0")) & (Demo$LILR995 %in% c("1/1", "0/1")))] <- "1/0"

Demo$Combine[which((Demo$APOL1Geno %in% c("1/1", "0/1")) & (Demo$LILR995 %in% c("0/0")))] <- "0/1"

Demo$Combine[which((Demo$APOL1Geno %in% c("0/0")) & (Demo$LILR995 %in% c("0/0")))] <- "0/0"


Demo$Combine2 <- "./."

Demo$Combine2[which((Demo$APOL1Geno %in% c("1/1")) & (Demo$LILR995 %in% c("1/1", "0/1")))] <- "1/1"

Demo$Combine2[which((Demo$APOL1Geno %in% c("0/0", "0/1")) & (Demo$LILR995 %in% c("1/1", "0/1")))] <- "1/0"

Demo$Combine2[which((Demo$APOL1Geno %in% c("1/1")) & (Demo$LILR995 %in% c("0/0")))] <- "0/1"

Demo$Combine2[which((Demo$APOL1Geno %in% c("0/0", "0/1")) & (Demo$LILR995 %in% c("0/0")))] <- "0/0"


write.table(Demo, file = paste0(File, ".geno.xls"), col.names = NA, quote = F, sep = "\t")

