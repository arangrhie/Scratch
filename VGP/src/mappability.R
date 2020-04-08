library(ggplot2)
library(scales)
library(cowplot)

getwd()
setwd("../VGP")

gray = "#EEEEEE"
black = "#808080"
my_col = c(gray, black)

dat = read.table("input/mappability.tab", header=T)
head(dat)

dat$Assembly = factor(dat$Assembly, levels = unique(dat$Assembly), labels = c("Taeniopygia_guttata-3.2.4", "bTaeGut1_v1.p"))

rna_dat = dat[dat$Type == "RNA" & (dat$Stat == "Total" | dat$Stat == "Unique"),]
atac_dat = dat[dat$Type == "ATAC" & (dat$Stat == "Total" | dat$Stat == "Unique"),]
head(rna_dat)

p <- ggplot(data=rna_dat, aes(x=Stat, y=Percent, fill=Assembly)) +
  geom_blank() +
  theme_classic() +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  scale_color_manual(values = my_col) +
  scale_fill_manual(values = my_col) +
  scale_y_continuous(expand = c(0, 0), limits=c(0, 100)) +
  theme(axis.text.y=element_text(face="bold", family = "Arial"),
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = "left")
p
legend <- get_legend(p)

plot_mappability <- function(dat, x_label) {
  ggplot(dat, aes(x=Stat, y=Percent, fill=Assembly)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Percent-sd, ymax=Percent+sd), width=.2,
                  position=position_dodge(.9)) +
    ylab("Avg. Reads Mapped (%)") +
    xlab(x_label) +
    theme_classic() +
    theme(#axis.title.x = element_blank(),
          axis.text = element_text(face="bold", family = "Arial"),
          legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_color_manual(values = my_col) +
    scale_fill_manual(values = my_col) +
    scale_y_continuous(expand = c(0, 0), limits=c(0, 100))
} 
a <- plot_mappability(rna_dat, "RNA-Seq")
b <- plot_mappability(atac_dat, "ATAC-Seq")
g <- arrangeGrob(legend, a, b, nrow = 1)
ggsave(file = "output/mappability.png", width = 8, height = 2.5, g)
