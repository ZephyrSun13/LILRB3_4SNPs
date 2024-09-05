
WES <- readRDS(file = "WES.rds")

DemoTab <- read.delim(file = "All_demo.xls", as.is = T)

rownames(DemoTab) <- paste0("P", DemoTab$person_id)

colnames(WES) <- paste0("P", colnames(WES))

length(intersect(colnames(WES), rownames(DemoTab)))


OverL <- intersect(colnames(WES), rownames(DemoTab))

DemoTab <- DemoTab[OverL, ]

DemoTab$LILR995 <- unlist(WES["chr19_54217126_C_T", rownames(DemoTab)])

DemoTab$G1 <- unlist(WES["chr22_36265988_T_G", rownames(DemoTab)])

DemoTab$G2 <- unlist((unlist(WES["chr22_36265995_AATAATT_A", rownames(DemoTab)])))

DemoTab$APOL1 <- NA

DemoTab$APOL1[DemoTab$G1 == "0/1" & DemoTab$G2 == "0/0"] <- "WT/G1"

DemoTab$APOL1[DemoTab$G1 == "0/0" & DemoTab$G2 == "0/1"] <- "WT/G2"

DemoTab$APOL1[DemoTab$G1 == "1/1"] <- "G1/G1"

DemoTab$APOL1[DemoTab$G2 == "1/1"] <- "G2/G2"

DemoTab$APOL1[DemoTab$G1 == "0/1" & DemoTab$G2 == "0/1"] <- "G1/G2"

DemoTab$APOL1[DemoTab$G1 == "0/0" & DemoTab$G2 == "0/0"] <- "WT/WT"

table(DemoTab[, c("race", "LILR995")])

table(DemoTab[, c("race", "APOL1")])


