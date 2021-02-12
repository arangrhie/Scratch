library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

getwd()
setwd("VGP")

## Fig. 1 Scaffolding strategy stats summary

# Color
black = "black"
gray = "#D3D3D3"
dark_gray = "#555555"
red = "#E41A1C"
blue = "#377EB8" # light blue = "#56B4E9"
green = "#4DAF4A"
purple = "#984EA3"  # purple = "#CC79A7"
orange = "#FF7F00"  # orange = "#E69F00"
yellow = "#FFFF33"

dat=read.table("input/Fig1_scaffolding.txt", header = TRUE)
head(dat)
#tail(dat)

# Rename and lock in order
dat$renamed=gsub("_", " + ", dat$Technologies)
dat$Technologies=gsub("=", " ", dat$renamed)
dat$Technologies
dat$Technologies = factor(dat$Technologies, levels = rev(dat$Technologies[order(dat$Sort)]))

# For plotting labels at the end of the bars
plot_text <- function(label_comma, x_aes, ymax) {
    # Adjust nudge_y for smaller values
    if (ymax > 1000) {
        nudge_y_val = 100
    } else if (ymax > 20) {
        nudge_y_val = 1
    } else {
        nudge_y_val = 0.2
    }
    
    geom.text.size = 6 * 0.35 # 1 pt ~ 0.35 mm
    if (label_comma) {
        geom_text(
            aes(label = comma(x_aes, accuracy = 1)),
            color = black,
            nudge_y = nudge_y_val,
            hjust = 0,
            size = geom.text.size
        )
    } else {
        geom_text(
            aes(label = x_aes),
            color = black,
            nudge_y = nudge_y_val,
            hjust = 0,
            size = geom.text.size
        )
    }
}

# Plot the bars flipped, with the Y axis labels
plot_bar_flipped_wi_ylabel <- function(x_aes=NULL, x_label=NULL, color=NULL, label_comma = TRUE) {
    # Get max X axis. Note the coord gets flipped
    ymax=max(x_aes)*1.42
    
    ggplot(data=dat, aes(x=Technologies, y=x_aes)) +
        theme_classic() +
        geom_bar(stat="identity", fill=color) +
        ylab(x_label) +
        plot_text(label_comma, x_aes, ymax) +
        theme(axis.title=element_text(size=7, face="bold"),
              axis.text=element_text(size=5, face="bold"),
              axis.title.y=element_blank(),
              axis.text.y=element_text(size=6, face="bold")) +
        coord_flip() +
        scale_y_continuous(labels = comma, limits = c(0, ymax))
}

# Plot the bars flipped
plot_bar_flipped <- function(x_aes=NULL, x_label=NULL, color=NULL, label_comma = TRUE) {
  
  # Get max X axis. Note the coord gets flipped
  ymax=max(x_aes)*1.42
  
  ggplot(data=dat, aes(x=Technologies, y=x_aes)) +
    theme_classic() +
    geom_bar(stat="identity", fill=color) +
    ylab(x_label) +
    plot_text(label_comma, x_aes, ymax) +
    theme(axis.title=element_text(size=7, face="bold"),
          axis.text=element_text(size=5, face="bold"),
          axis.title.y=element_blank(),       # Remove Y labels
          axis.text.y=element_blank()) +
    coord_flip() +
    scale_y_continuous(labels = comma, limits = c(0, ymax))

}

# Fig. 1a-d
a = plot_bar_flipped_wi_ylabel(dat$ContigNG50, "Contig NG50 (Mb)", dark_gray, label_comma = FALSE)
b = plot_bar_flipped(dat$ScaffoldNG50, "Scaffold NG50 (Mb)", orange, label_comma = FALSE)
c = plot_bar_flipped(dat$NumGaps, "Num. Gaps", gray)
d = plot_bar_flipped(dat$NumMisjoins, "Num. Mis-joins", red)
grid.arrange(a, b, c, d, nrow = 1, widths = c(1.8,1,1,1))

g <- arrangeGrob(a, b, c, d, nrow = 1, widths = c(1.8,1,1,1))
ggsave(file = "output/Fig1a_d.pdf", width = 170, height = 50, g, units = "mm")


###########################################################################
# Fig. 1h Chromosome Size Comparison between Assembled vs. Karyotype
dat=read.table("input/Fig1_karyotype.txt", header = TRUE)
head(dat)

# Highlight Z and W
Z=dat[dat$Chromosome=="Z",]
W=dat[dat$Chromosome=="W",]

# Log-log scale
ggplot(dat, aes(x = AssembledChrSize, y = KaryotypeSize)) +
  geom_smooth(method='lm', se=F, color = "grey", size = 0.5) +
  geom_point(size=0.5) + theme_classic() + 
  geom_point(data = Z, size=0.6, color = "red") +
  geom_point(data = W, size=0.6, color = "red") +
  scale_x_log10(expand = c(0, 0.05), limits=c(1.5, 210),
                breaks = c(3, 10, 30, 100)) +
  scale_y_log10(expand = c(0, 0.05), limits=c(1.5, 210),
                breaks = c(3, 10, 30, 100)) +
  xlab("Assembled Chr. Size (Mb)") + ylab("Karyotype Size (Mb)") +
  theme(axis.title = element_text(size=6, face = "bold"),
        axis.text  = element_text(size=6))
ggsave("output/Fig1_karyotype_loglog.png", width = 1.5, height = 1.5)
ggsave("output/Fig1_karyotype_loglog.pdf", width = 40, height = 40, units = "mm", device=cairo_pdf)

# Absolute scale
ggplot(dat, aes(x = AssembledChrSize, y = KaryotypeSize)) +
    geom_smooth(method='lm', se=F, color = "grey") +
    geom_point() + theme_classic() + 
    scale_x_continuous(expand = c(0, 0), limits=c(0, 210)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(0, 210)) +
    xlab("Assembled Chr Size (Mb)") + ylab("Karyotype Size (Mb)")
fit=lm(KaryotypeSize ~ AssembledChrSize, data = dat)
summary(fit)
ggsave("output/Fig1_karyotype.png", width = 3, height = 3)