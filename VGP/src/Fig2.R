library(ggplot2)
library(scales)

# Color: Set1
red = "#E41A1C"
green = "#4DAF4A"
blue = "#377EB8" # light blue = "#56B4E9"
purple = "#984EA3"  # purple = "#CC79A7"
gray = "#808080"
gray_dark = "#404040"
black = "#000000"
my_col3 = c(gray, gray_dark, black)
my_col = c(gray, black)

fit_log <- function(dat, x, y) {
  log.model <-lm(log(y) ~ x, dat)
  log.model.dat <- data.frame(x = x, y = exp(fitted(log.model)))
  return(log.model.dat)
}

LINE=0.3
POINT=0.2

plot_dots_exp <- function (dat, log.model, xlabel, ylabel, xmin, xmax, ymin, ymax) {
  p <-  ggplot() +
    theme_classic() +
    geom_segment(data = dat, aes(x = x, y = Before, xend = x, yend = After), color = "gray", size = LINE) +
    geom_line(data = log.model.before, aes(x, y), colour = "grey", size=LINE) +
    geom_line(data = log.model.after, aes(x, y), color = black, size=LINE) +
    geom_point(dat, mapping = aes(x = x, y = Before), colour = gray, size=POINT) +
    geom_point(dat, mapping = aes(x = x, y = After),  colour = black, size=POINT) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=5),
          axis.ticks = element_line(size=LINE),
          axis.line  = element_line(size=LINE)) +
    scale_x_continuous(expand = c(0, 0), limits=c(xmin, xmax)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(ymin, ymax)) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}

plot_dots <- function (dat, xlabel, ylabel, xmin, xmax, ymin, ymax) {
  p <-  ggplot() +
    theme_classic() +
    geom_segment(data = dat, aes(x = x, y = Before, xend = x, yend = After), color = "gray", size = LINE) +
    geom_smooth(data = dat, aes(x = x, y = Before), method='lm', se=F, color = "gray", size = LINE) +
    geom_smooth(data = dat, aes(x = x, y = After), method='lm', se=F, color = "black", size = LINE) +
    geom_point(dat, mapping = aes(x = x, y = Before), colour = gray, size = POINT) +
    geom_point(dat, mapping = aes(x = x, y = After),  colour = black, size = POINT) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=5),
          axis.ticks = element_line(size=LINE),
          axis.line  = element_line(size=LINE)) +
    scale_x_continuous(expand = c(0, 0), limits=c(xmin, xmax)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(ymin, ymax)) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}

### Plot only "After"
plot_dot <- function (dat, xlabel, ylabel, xmin, xmax, ymin, ymax) {
  p <-  ggplot() +
    theme_classic() +
    geom_smooth(data = dat, aes(x = x, y = After), method='lm', se=F, color = "black", size = LINE) +
    geom_point(dat, mapping = aes(x = x, y = After),  colour = black, size = POINT) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=5),
          axis.ticks = element_line(size=LINE),
          axis.line  = element_line(size=LINE)) +
    scale_x_continuous(expand = c(0, 0), limits=c(xmin, xmax)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(ymin, ymax)) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}

getwd()
setwd("../VGP")

dat=read.table("input/Fig2/gsize_rep.tab", header=T)
head(dat)
tail(dat)

#### Genome size and scaffold NG50
p1 <- ggplot(dat, aes(x = GenomeSize, y = ScaffoldNG50)) +
  geom_smooth(method='lm', se=F, color = black, alpha = 0.2, size=LINE) +
  geom_point(size=POINT) +
  theme_classic() +
  theme(axis.title = element_text(size=6, face = "bold"),
        axis.text  = element_text(size=5),
        axis.ticks = element_line(size=LINE),
        axis.line  = element_line(size=LINE)) +
  scale_x_continuous(expand = c(0, 0), limits=c(0.0, 6.0)) +  # for both axes to remove unneeded padding
  scale_y_continuous(expand = c(0, 0), limits=c(0.0, 600.0)) +
  xlab("Genome Size (Gbp)") + ylab("Scaffold NG50 (Mbp)")
#ggsave("output/Fig2c_genomesize_scaffNG50.png", width = 1.5, height = 1.5)
fit=lm(ScaffoldNG50 ~ GenomeSize, data = dat)
summary(fit)
# Take adjusted R-square and p-value from F-statistic

#dat=read.table("input/initial_submit/Fig3_purging.tab", header=T)
dat=read.table("input/Fig2/purging.tab", header=T)
head(dat)
tail(dat)

#### Repeat and Contig NG50, before / after curation
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       Before = dat[dat$Purging == "Before",]$ContigNG50,
                       After = dat[dat$Purging == "Curation",]$ContigNG50)
log.model.before <- fit_log(dat_diff, dat_diff$x, dat_diff$Before)
log.model.after <- fit_log(dat_diff, dat_diff$x, dat_diff$After)

p2 <- plot_dots_exp(dat_diff, log.model, "Repeat (%)", "Contig NG50 (Mbp)", 0, 70, 0, 30)
#ggsave("output/Fig3c_rep_contigNG50.png", width = 1.5, height = 1.5)
summary(lm(log(Before) ~ x, data = dat_diff))
summary(lm(log(After) ~ x, data = dat_diff))

### Repeat and gaps
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       After = dat[dat$Purging == "Curation",]$NumGap)
p3 <- plot_dot(dat_diff, "Repeat (%)", "Num. Gaps / Gb", 0, 60, 0, 1500)
summary(lm(After ~ x, data = dat_diff))
#ggsave("output/Fig2d_rep_gap.png", width = 1.5, height = 1.5)

#### heterozygosity and primary assembly size
log.model.dat.before <- fit_log(dat[dat$Purging == "Before",], dat[dat$Purging == "Before",]$Heterozygosity, dat[dat$Purging == "Before",]$Primary)
log.model.dat.after <- fit_log(dat[dat$Purging == "After",], dat[dat$Purging == "After",]$Heterozygosity, dat[dat$Purging == "After",]$Primary)

dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$Primary,
                       After = dat[dat$Purging == "After",]$Primary)
p4 <- plot_dots(dat_diff, "Heterozygosity (%)", "Primary Size (%)", 0.0, 1.8, 70, 110) +
  geom_hline(yintercept=100, linetype="dashed", color = "gray")
#ggsave("output/Fig2e_het_prim.png", width = 1.5, height = 1.5)
fit=lm(Primary ~ Heterozygosity, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(Primary ~ Heterozygosity, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and k-mer duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$DupKmer,
                       After = dat[dat$Purging == "After",]$DupKmer)
p5 <- plot_dots(dat_diff, "Heterozygosity (%)", "k-mer dup. (%)", 0.0, 1.8, 0, 30)
#ggsave("output/Fig2f_het_dupKmer.png", width = 1.5, height = 1.5)
fit=lm(DupKmer ~ Heterozygosity, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupKmer ~ Heterozygosity, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and BUSCO duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$DupBUSCO,
                       After = dat[dat$Purging == "After",]$DupBUSCO)
p6 <- plot_dots(dat_diff, "Heterozygosity (%)", "BUSCO dup. (%)", 0.0, 1.8, 0, 20)
#ggsave("output/Fig2g_het_dupBUSCO.png", width = 1.5, height = 1.5)
fit=lm(DupBUSCO ~ Heterozygosity, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupBUSCO ~ Heterozygosity, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and k-mer duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       Before = dat[dat$Purging == "Before",]$DupKmer,
                       After = dat[dat$Purging == "After",]$DupKmer)
dat_diff
p7 <- plot_dots(dat_diff, "Repeat (%)", "k-mer dup. (%)", 0, 60, 0, 30)
#ggsave("output/Fig2h_rep_dupKmer.png", width = 1.5, height = 1.5)
fit=lm(DupKmer ~ Repeat, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupKmer ~ Repeat, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and BUSCO duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       Before = dat[dat$Purging == "Before",]$DupBUSCO,
                       After = dat[dat$Purging == "After",]$DupBUSCO)
p8 <- plot_dots(dat_diff, "Repeat (%)", "BUSCO dup. (%)", 0, 60, 0, 20)
#ggsave("output/Fig2i_rep_dupBUSCO.png", width = 1.5, height = 1.5)
fit=lm(DupBUSCO ~ Repeat, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupBUSCO ~ Repeat, data = dat[dat$Purging == "After",])
summary(fit)

# Aggregate all to one figure
g <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8,  nrow = 2)
ggsave(file = "output/pub/Fig2.pdf", width = 120, height = 60, g, units = "mm", device=cairo_pdf)

################### Figures not used  ####################

#### heterozygosity and Contig NG50, before / after curation
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$ContigNG50,
                       After = dat[dat$Purging == "Curation",]$ContigNG50)
log.model.before <- fit_log(dat_diff, dat_diff$x, dat_diff$Before)
log.model.after <- fit_log(dat_diff, dat_diff$x, dat_diff$After)
summary(lm(log(Before) ~ x, data = dat_diff))
summary(lm(log(After) ~ x, data = dat_diff))
plot_dots_exp(dat_diff, log.model, "heterozygosity (%)", "Contig NG50 (Mbp)", 0, 1.8, 0, 30)
ggsave("output/Fig2X_het_contigNG50.png", width = 2, height = 2)

#### Heterozygosity and alternate assembly size
dat=read.table("input/deprecated/het_asmsize.tab", header=T)
ggplot(dat, aes(x = Heterozygosity, y = Alts)) +
  geom_hline(yintercept=100, linetype="dashed", color = "grey") +
  geom_smooth(method='lm', formula= (y ~ log(x)), se=F, color = "grey") +
  geom_point() + theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits=c(0,120)) +
  scale_x_continuous(expand = c(0, 0), limits=c(0.00,1.80)) +  # for both axes to remove unneeded padding
  xlab("Heterozygosity (%)") + ylab("Alt. assembly size (%)")
fit=lm(Alts ~ log(Heterozygosity), data = dat)
summary(fit)
ggsave("output/Fig2X_het_alt.png", width = 2, height = 2)

