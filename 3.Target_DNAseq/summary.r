
Tab.lst <- unlist(read.table(file = "Tabs", as.is = T))

Tabs <- lapply(Tab.lst, function(x) read.delim(file = x, row.names = 1))

names(Tabs) <- sapply(Tab.lst, function(x) paste(unlist(strsplit(unlist(strsplit(x, split = "/"))[8], split = "_"))[1:3], collapse = "."))

Tabs <- lapply(Tabs, function(x) {rownames(x) <- sapply(rownames(x), function(y){paste(unlist(strsplit(y, split = " "))[1:2], collapse = "_")}); x})


Mat <- matrix(nrow = 400, ncol = length(Tabs), dimnames = list(paste0("Target_", 1:400), names(Tabs)))

for(nn in names(Tabs)){

	Mat[rownames(Tabs[[nn]]), nn] <- as.character(Tabs[[nn]][[1]])

}

write.table(Mat, file = "Mat.DP4.xls", quote = F, sep = "\t", col.names = NA)


Target160Rate <- sapply(Mat["Target_160", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

Target160Count <- sapply(Mat["Target_160", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})

Target188Rate <- sapply(Mat["Target_188", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

Target188Count <- sapply(Mat["Target_188", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})

Target199Rate <- sapply(Mat["Target_199", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

Target199Count <- sapply(Mat["Target_199", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})

rs549267286Rate <- sapply(Mat["Target_200", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

rs549267286Count <- sapply(Mat["Target_200", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})

Target201Rate <- sapply(Mat["Target_201", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

Target201Count <- sapply(Mat["Target_201", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})

Target202Rate <- sapply(Mat["Target_202", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

Target202Count <- sapply(Mat["Target_202", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})

Target242Rate <- sapply(Mat["Target_242", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); (Tmp[3] + Tmp[4])/(Tmp[1] + Tmp[2])})

Target242Count <- sapply(Mat["Target_242", ], function(x) {Tmp <- as.numeric(unlist(strsplit(x, split = ","))); paste0("(", Tmp[3], "+", Tmp[4], ")", "/", "(", Tmp[1], "+", Tmp[2], ")")})


SamMap <- read.delim(file = "SamMap", sep = "\t", as.is = T, row.names = 2)

rs549267286Tab <- data.frame(Sam = paste0(SamMap[names(rs549267286Rate), ], "Rb"), SeqID = names(rs549267286Rate), Rate = rs549267286Rate, Count = rs549267286Count[names(rs549267286Rate)], Rate_160 = Target160Rate[names(rs549267286Rate)], Count_160 = Target160Count[names(rs549267286Rate)], Rate_188 = Target188Rate[names(rs549267286Rate)], Count_188 = Target188Count[names(rs549267286Rate)], Rate_199 = Target199Rate[names(rs549267286Rate)], Count_199 = Target199Count[names(rs549267286Rate)], Rate_201 = Target201Rate[names(rs549267286Rate)], Count_201 = Target201Count[names(rs549267286Rate)], Rate_202 = Target202Rate[names(rs549267286Rate)], Count_202 = Target202Count[names(rs549267286Rate)], Rate_242 = Target242Rate[names(rs549267286Rate)], Count_242 = Target242Count[names(rs549267286Rate)])

write.table(rs549267286Tab, file = "rs549267286Tab.xls", sep = "\t", quote = F, col.names = NA)

