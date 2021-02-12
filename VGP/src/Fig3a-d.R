library(ggplot2)
library(scales)
library(cowplot)
library(gridExtra)

getwd()
setwd("../VGP")

gray = "#EEEEEE"
black = "#808080"
gray_light = "#ffffff"
black_light = "#a6a6a6"
my_col = c(gray, black)
my_col_3 = c(gray, black, black)
my_col_4 = c(gray_light, black_light, black_light, gray, black, black)

plot_mappability <- function(dat, label) {
  ggplot(dat, aes(x=Stat, y=Percent, fill=Assembly)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Percent-sem, ymax=Percent+sem), width=.2,
                  position=position_dodge(.9)) +
    ylab("Avg. Reads Mapped (%)") +
    #xlab(label) +
    theme_classic() +
    theme(
      axis.title = element_text(face="bold", size=6),
      axis.text = element_text(face="bold", size=5),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_blank()) +
    scale_color_manual(values = my_col) +
    scale_fill_manual(values = my_col) +
    scale_x_discrete(name=label,
                     labels=c("Total\n", "Unique\n")) +
    scale_y_continuous(expand = c(0, 0), limits=c(50, 100), oob = rescale_none)
}

plot_annotation_2 <- function(dat, label, y_limit) {
  ggplot(dat, aes(x=Order, y=Count, fill=Label)) + 
    geom_bar(position=position_stack(), stat="identity", color="black") +
    theme_classic() +
    theme(
      axis.title = element_text(face="bold", size=6),
      axis.text = element_text(face="bold", size=5),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_blank()) +
    scale_fill_manual(values = my_col_4) +
    scale_x_continuous(name=label,
                       breaks=c(1.5, 5, 8.5),
                       labels=c("Anna's\nHummingbird", "Zebra\nfinch", "Platypus")) +
    scale_y_continuous(expand = c(0, 0), labels = comma, limits=c(0, y_limit))
}

plot_annotation <- function(dat, label, y_limit) {
  ggplot(dat, aes(x=Order, y=Count, fill=Label)) + 
    geom_bar(position="dodge", stat="identity", color="black") +
    #xlab(label) +
    theme_classic() +
    theme(
      axis.title = element_text(face="bold", size=6),
      axis.text = element_text(face="bold", size=5),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_blank()) +
    scale_fill_manual(values = my_col_3) +
    scale_x_continuous(name=label,
                       breaks=c(1.5, 5, 8.5),
                       labels=c("Anna's\nHummingbird", "Zebra\nfinch", "Platypus")) +
    scale_y_continuous(expand = c(0, 0), labels = comma, limits=c(0, y_limit))
}

dat = read.table("input/Fig3/mappability.tab", header=T)
head(dat)

# Labels
# dat$Assembly = factor(dat$Assembly, levels = unique(dat$Assembly), labels = c("Taeniopygia_guttata-3.2.4", "bTaeGut1_v1.p"))
dat$Assembly = factor(dat$Assembly, levels = unique(dat$Assembly), labels = c("Previous", "VGP"))

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
  theme(
    axis.title = element_text(face="bold", size=5),
    axis.text = element_text(face="bold", size=5),
    axis.line.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    legend.position = "left",
    legend.text = element_text(size = 5),
    legend.title = element_text(face="bold", size = 5))
p
legend <- get_legend(p)

a <- plot_mappability(rna_dat, "RNA-Seq")
b <- plot_mappability(atac_dat, "ATAC-Seq")

## Annotation results ##
dat = read.table("input/Fig3/annotation.tab", header=T)
CDS_dat=dat[dat$Category !="PartialCodingGenes",]
CDS_dat$Label = factor(CDS_dat$Label,
                       levels = c("Previous_Total", "VGP_Total", "VGPTrio_Total", "Previous_Fully", "VGP_Fully", "VGPTrio_Fully"))

partial_dat=dat[dat$Category=="PartialCodingGenes",]

head(CDS_dat)
tail(CDS_dat)

c <- plot_annotation_2(dat = CDS_dat, label = "Total vs. Fully supp. CDS", y_limit = 50000)
d <- plot_annotation(partial_dat, "Partial Coding Genes", 6000)

# g <- arrangeGrob(legend, a, b, c, d, nrow = 1)
plot_grid(legend, a, b, c, d, nrow = 1, rel_widths = c(0.08, 0.18, 0.18, 0.28, 0.28))
ggsave(file = "output/pub/Fig3ad.pdf", width = 120, height = 30, units = "mm")
#ggsave(file = "output/Fig5ad.png", width = 6.5, height = 1.5)

