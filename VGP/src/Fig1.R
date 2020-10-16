library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

getwd()
setwd("../VGP")
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

# Lock in order
dat$Technologies = factor(dat$Technologies, levels = rev(dat$Technologies[order(dat$Sort)]))

# Plot only the y label, seems like the order is still alphabetical though. ggplot2 bug?
plot_ylab_only <- function(x_aes=NULL) {
  ggplot(data=dat, aes(x=gsub("_", " + ", Technologies), y=x_aes)) +
    geom_blank() +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_text(size=8, face="bold", family = "Arial"),
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    coord_flip()
}

ylabels = plot_ylab_only(dat$ContigNG50)
ylabels # Why is the order alphabetical :(
#ggsave(file = "output/Fig1a_d_y.png", width = 2, height = 4, ylabels)

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
  
  if (label_comma) {
    geom_text(aes(label=comma(x_aes, accuracy = 1)), color=black, nudge_y = nudge_y_val, hjust = 0, size = 3)
  } else {
    geom_text(aes(label=x_aes), color=black, nudge_y = nudge_y_val, hjust = 0, size = 3)
  }
}

# Plot the bars flipped
plot_bar_flipped <- function(x_aes=NULL, x_label=NULL, color=NULL, label_comma = TRUE) {
  
  # Get max X axis. Note the coord gets flipped
  ymax=max(x_aes)*1.4
  
  ggplot(data=dat, aes(x=Technologies, y=x_aes)) +
    theme_classic() +
    geom_bar(stat="identity", fill=color) +
    ylab(x_label) +
    plot_text(label_comma, x_aes, ymax) +
    theme(axis.title=element_text(size=10, face="bold"),
          axis.text=element_text(size=8, face="bold"),
          axis.title.y=element_blank(),       # Remove Y labels
          axis.text.y=element_blank()) +
    coord_flip() +
    scale_y_continuous(labels = comma, limits = c(0, ymax))

}

# Fig. 1a-d
a = plot_bar_flipped(dat$ContigNG50, "Contig NG50 (Mb)", dark_gray, label_comma = FALSE)
b = plot_bar_flipped(dat$ScaffoldNG50, "Scaffold NG50 (Mb)", orange, label_comma = FALSE)
c = plot_bar_flipped(dat$NumGaps, "Num. Gaps", gray)
d = plot_bar_flipped(dat$NumMisjoins, "Num. Mis-joins", red)
grid.arrange(a, b, c, d, nrow = 1)

g <- arrangeGrob(a, b, c, d, nrow = 1)
ggsave(file = "output/Fig1a_d.png", width = 8, height = 3, g)

# Fig. 1h Chromosome Size Comparison between Assembled vs. Karyotype
dat=read.table("input/Fig1_karyotype.txt", header = TRUE)
head(dat)

ggplot(dat, aes(x = AssembledChrSize, y = KaryotypeSize)) +
  geom_smooth(method='lm', se=F, color = "grey") +
  geom_point() + theme_classic() + 
  scale_x_continuous(expand = c(0, 0), limits=c(0, 210)) +  # for both axes to remove unneeded padding
  scale_y_continuous(expand = c(0, 0), limits=c(0, 210)) +
  xlab("Assembled Chr Size (Mb)") + ylab("Karyotype Size (Mb)")
fit=lm(KaryotypeSize ~ AssembledChrSize, data = dat)
summary(fit)
ggsave("output/Fig1_karyotype.png", width = 3, height = 3)

fancy_scientific <- function(d) {
  # turn in to character string in scientific notation
  d <- format(d, scientific = TRUE)
  # quote the part before the exponent to keep all the digits and turn the 'e+' into 10^ format
  d <- gsub("^(.*)e\\+", "'\\1'%*%10^", d)
  # convert 0x10^00 to 0
  d <- gsub("\\'0[\\.0]*\\'(.*)", "'0'", d)
  # return this as an expression
  parse(text=d)
}

W=dat[dat$Chromosome=="W",]
W

Z=dat[dat$Chromosome=="Z",]
Z


ggplot(dat, aes(x = AssembledChrSize, y = KaryotypeSize)) +
  geom_smooth(method='lm', se=F, color = "grey", size = 0.5) +
  geom_point(size=0.5) + theme_classic() + 
  geom_point(data = Z, size=0.6, color = "red") +
  geom_point(data = W, size=0.6, color = "orange") +
  scale_x_log10(expand = c(0, 0.05), limits=c(1.5, 210),
                breaks = c(3, 10, 30, 100)) +
  scale_y_log10(expand = c(0, 0.05), limits=c(1.5, 210),
                breaks = c(3, 10, 30, 100)) +
  xlab("Assembled Chr. Size (Mb)") + ylab("Karyotype Size (Mb)") +
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=6, family = "Arial"))
ggsave("output/Fig1_karyotype_loglog.png", width = 1.5, height = 1.5)

ggplot(dat, aes(x = AssembledChrSize, y = KaryotypeSize)) +
  geom_smooth(method='lm', se=F, color = "grey") +
  geom_point() + theme_classic() + 
  scale_x_log10(expand = c(0, 0.05), limits=c(1.5, 210),
                breaks = c(3, 10, 30, 100)) +
  scale_y_log10(expand = c(0, 0.05), limits=c(1.5, 210),
                breaks = c(3, 10, 30, 100)) +
  xlab("Assembled Chr. Size (Mb)") + ylab("Karyotype Size (Mb)") +
  theme(axis.title = element_blank(),
        axis.text  = element_text(size=12, family = "Arial"))
ggsave("output/Fig1_karyotype_loglog_small.png", width = 1.5, height = 1.5)

breaks = trans_breaks("log10", function(x) 10^x)
labels = trans_format("log10", math_format(10^.x))
