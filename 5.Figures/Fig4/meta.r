
library(survival)
library(survminer)
library(ggplot2)
library(RColorBrewer)
library(metafor)

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", row.names = 1, as.is = T, header = F)


#GOCAR <- read.delim(file = "../2.Align/AllAfrican_American_SurvPlot.cox.xls", as.is = T)

GOCAR <- read.delim(file = "../2.Align/Tab_APOL1.adj.xls", as.is = T)

#CTOT <- read.delim(file = "../2.Align_CTOT/Target_SurvPlot.cox.xls", as.is = T)

#Pitt <- read.delim(file = "../2.Align_Pitt/Target_SurvPlot.cox.xls", as.is = T)

Pitt <- read.delim(file = "../2.Align_Pitt/Target_SurvPlot.cox.xls")

dat <- data.frame(study = c("GOCAR", "Pitt"),
                  logHR = c(GOCAR[9, "coef"], Pitt[2, "coef"]),
                  selogHR = c(GOCAR[9, "se.coef."], Pitt[2, "se.coef."]),
                  n = c(83, 54),
                  stringsAsFactors = FALSE)



## fixed effect    

res1 <- rma(yi = logHR, sei = selogHR, data = dat, method = "FE")


## random effect

res2 <- rma(yi = logHR, sei = selogHR, data = dat, method = "REML")


## forest plot
res.plot <- res1

# pdf(file = "./res/DCGL/meta_GOCAR_CTOT_DCGL_fixed.pdf", width = 6, height = 3)
# png(filename = "./res/DCGL/meta_GOCAR_CTOT_DCGL_fixed.png",
    # width = 6, height = 3, units = "in", res = 300)
# tiff(filename = "./res/DCGL/meta_GOCAR_CTOT_DCGL_fixed.tiff",
    # width = 6, height = 3, units = "in", res = 300)

pdf("Meta.pdf", height = 4)

op <- par(mar = c(4, 0.5, 0.1, 0.5))
forest(res.plot, atransf = exp,
       slab = c("GoCAR", "SIRPA"),, #c("GOCAR (n = 343)", "CTOT (n = 117)"),
       header = c("Study", "Hazard Ratio [95% CI]"),
       xlab = "Hazard Ratio",
       mlab = "Meta-analysis (p = 0.0008)",
       xlim = c(-3, 6),
       at = log(c(0.5, 1:9))
       )
par(op)

dev.off()



Result <- matrix(NA, ncol = 4, nrow = 5, dimnames = list(c("Title", "GoCAR", "SIRPA", "Space", "Meta"), c("lower.95", "upper.95", "OR", "P")))

Result["GoCAR", ] <- unlist(GOCAR[9, c(2, 3, 5, 8)])

Result["SIRPA", ] <- unlist(Pitt[2, c(2, 3, 5, 8)])

Result["Meta", ] <- c(exp(c(res1$ci.lb, res1$ci.ub, res1$beta)), res1$pval)


Tab <- data.frame(Result)

Tab$Label <- factor(rownames(Tab), levels = rev(rownames(Tab)))

#Tab$colour <- rep(c("white", "gray95"), 3)[1:nrow(Tab)]

Tab$colour <- rep(c("white", Colors["HeatYellow1", ]), 3)[1:nrow(Tab)]


Tab.man <- Tab

#Tab.point <- Tab[which(Tab$upper.95 > 9), ]

#Tab.point$upper.95 <- 9

#Tab.man$upper.95[which(Tab.man$upper.95 > 9)] <- 9


fp <- ggplot(data=Tab.man, aes(y=Label, x=OR, xmin=lower.95, xmax=upper.95)) +

        geom_hline(aes(yintercept = Label, colour = colour), size = 10) +

        geom_pointrange(color = Colors["Red", ], size = 0.7) +

        geom_vline(xintercept=1, lty=2) +  # add a dotted line at x=1 after flip

        #geom_point(data = Tab.point, aes(y = Label, x = upper.95), shape = 62, color = Colors["Red", ], size = 4) +

        xlim(NA, 9) +

        #coord_trans(x = "log2") +

        scale_x_continuous(trans='log2') +

        scale_y_discrete(breaks = rownames(Tab.man), labels = c("", "GoCAR", "SIRPA", "", "Meta analysis")) +

        #coord_flip() +  # flip coordinates (puts labels on y axis)

        ylab("") + xlab("Hazard ratio (95% CI)") +

        theme_classic() + scale_colour_identity() + # use a white background

        theme(axis.ticks = element_blank(), axis.text.y = element_text(size = 10))


Tab$Range <- paste0(round(Tab$OR, 2), " (", round(Tab$lower.95, 2), ", ", round(Tab$upper.95, 2), ")")

Tab$Range[Tab$Range == "NA (NA, NA)"] <- NA

Tab$P <- as.character(round(Tab$P, 3))

Tab$P[Tab$P == "0"] <- "<0.001"

Tab$P[1] <- "P value"

Tab$Range[1] <- "Hazard ratio (95% CI)"

fp2 <- ggplot(data = Tab, aes(y = Label)) +

        geom_hline(aes(yintercept = Label, colour = colour), size = 10) +

        #geom_text(aes(x = 0, label = Label), hjust = 0) +

        geom_text(aes(x = 0, label = P), hjust = 0, size = 3) +

        geom_text(aes(x = 3, label = Range), hjust = 1, size = 3) +

        scale_colour_identity() +

        theme_void() +

        theme(plot.margin = margin(5, 0, 32, 0))


ggsave("Meta_DCGL.Cox.adj.pdf", ggarrange(fp, fp2, ncol = 2, widths = c(2.8, 1.2)), width = 7, height = 3)

