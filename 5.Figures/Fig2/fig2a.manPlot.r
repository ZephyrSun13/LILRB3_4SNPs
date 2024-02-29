
Args <- commandArgs(trailingOnly = T)

#library(qqman)
library(ggplot2)
library(RColorBrewer)
library(dplyr)


Tab <- read.delim(file = Args[1], as.is = T, row.names = 1)

#Tab <- Tab[Tab$Pr...z.. >= 1e-20, ]

gwasResults <- data.frame(SNP = rownames(Tab), CHR = sapply(rownames(Tab), function(x) gsub("chr", "", unlist(strsplit(x, split = ":"))[1])), Pos = sapply(rownames(Tab), function(x) as.numeric(unlist(strsplit(x, split = ":"))[2])), P = Tab$Pr...z..)

gwasResults$CHR <- as.character(gwasResults$CHR)

gwasResults <- gwasResults[with(gwasResults, order(CHR, Pos)), ]

Counts <- table(gwasResults$CHR)

gwasResults$BP <- unlist(lapply(gwasResults$CHR[!duplicated(gwasResults$CHR)], function(x) 1:Counts[x]))


gwasResults <- gwasResults[!is.na(gwasResults$P), ]

gwasResults$CHR <- factor(as.character(gwasResults$CHR), levels = c(as.character(1:22), "X", "Y"))

don <- gwasResults %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

#don$CHR <- factor(as.character(don$CHR), levels = c(as.character(1:22), "X", "Y"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", header = F, row.names = 1, as.is = T)

pp <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=CHR), alpha=0.8, size=0.6) +
    scale_color_manual(values = rep(Colors[c("Blue2", "Orange"), ], 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
    # Custom the theme:

    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = Colors["Red", ], size = 1) +

    xlab("Chromosome") +

    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 15)
    )

ggsave(paste0(Args[1], ".man.pdf"), width = 9, height = 4)

