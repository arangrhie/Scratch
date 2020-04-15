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

plot_dots_exp <- function (dat, log.model, xlabel, ylabel, xmin, xmax, ymin, ymax) {
  p <-  ggplot() +
    theme_classic() +
    geom_segment(data = dat, aes(x = x, y = Before, xend = x, yend = After), color = "gray") +
    geom_line(data = log.model.before, aes(x, y), size = 1, colour = "grey") +
    geom_line(data = log.model.after, aes(x, y), size = 1, color = black) +
    geom_point(dat, mapping = aes(x = x, y = Before), colour = gray) +
    geom_point(dat, mapping = aes(x = x, y = After),  colour = black) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_x_continuous(expand = c(0, 0), limits=c(xmin, xmax)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(ymin, ymax)) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}

plot_dots <- function (dat, xlabel, ylabel, xmin, xmax, ymin, ymax) {
  p <-  ggplot() +
    theme_classic() +
    geom_segment(data = dat, aes(x = x, y = Before, xend = x, yend = After), color = "gray") +
    geom_smooth(data = dat, aes(x = x, y = Before), method='lm', se=F, color = "gray") +
    geom_smooth(data = dat, aes(x = x, y = After), method='lm', se=F, color = "black") +
    geom_point(dat, mapping = aes(x = x, y = Before), colour = gray) +
    geom_point(dat, mapping = aes(x = x, y = After),  colour = black) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_x_continuous(expand = c(0, 0), limits=c(xmin, xmax)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(ymin, ymax)) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}

### Plot only "After"
plot_dot <- function (dat, xlabel, ylabel, xmin, xmax, ymin, ymax) {
  p <-  ggplot() +
    theme_classic() +
    geom_smooth(data = dat, aes(x = x, y = After), method='lm', se=F, color = "black") +
    geom_point(dat, mapping = aes(x = x, y = After),  colour = black) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.key = element_blank()) +
    scale_x_continuous(expand = c(0, 0), limits=c(xmin, xmax)) +  # for both axes to remove unneeded padding
    scale_y_continuous(expand = c(0, 0), limits=c(ymin, ymax)) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}

getwd()
setwd("../VGP")

dat=read.table("input/gsize_rep.tab", header=T)
head(dat)
tail(dat)

#### Genome size and scaffold NG50
ggplot(dat, aes(x = GenomeSize, y = ScaffoldNG50)) +
  geom_smooth(method='lm', se=F, size = 1, color = black, alpha = 0.2) +
  geom_point() +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), limits=c(0.0, 6.0)) +  # for both axes to remove unneeded padding
  scale_y_continuous(expand = c(0, 0), limits=c(0.0, 600.0)) +
  xlab("Genome Size (Gb)") + ylab("Scaffold NG50 (Mb)")
ggsave("output/genomesize_scaffNG50.png", width = 2, height = 2)
fit=lm(ScaffoldNG50 ~ GenomeSize, data = dat)
summary(fit)


dat=read.table("input/purging.tab", header=T)
head(dat)
tail(dat)

#### Repeat and Contig NG50, before / after curation
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       Before = dat[dat$Purging == "Before",]$ContigNG50,
                       After = dat[dat$Purging == "Curation",]$ContigNG50)
log.model.before <- fit_log(dat_diff, dat_diff$x, dat_diff$Before)
log.model.after <- fit_log(dat_diff, dat_diff$x, dat_diff$After)

p <- plot_dots_exp(dat_diff, log.model, "Repeat (%)", "Contig NG50 (Mb)", 0, 70, 0, 30)
p
ggsave("output/rep_contigNG50.png", width = 2, height = 2)
summary(lm(log(Before) ~ x, data = dat_diff))
summary(lm(log(After) ~ x, data = dat_diff))

#### heterozygosity and Contig NG50, before / after curation
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$ContigNG50,
                       After = dat[dat$Purging == "Curation",]$ContigNG50)
log.model.before <- fit_log(dat_diff, dat_diff$x, dat_diff$Before)
log.model.after <- fit_log(dat_diff, dat_diff$x, dat_diff$After)

p <- plot_dots_exp(dat_diff, log.model, "heterozygosity (%)", "Contig NG50 (Mb)", 0, 1.8, 0, 30)
p
summary(lm(log(Before) ~ x, data = dat_diff))
summary(lm(log(After) ~ x, data = dat_diff))
ggsave("output/het_contigNG50.png", width = 2, height = 2)

### Repeat and gaps
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       After = dat[dat$Purging == "Curation",]$NumGap)
p <- plot_dot(dat_diff, "Repeat (%)", "Num. Gaps / Gb", 0, 60, 0, 1500)
p
summary(lm(After ~ x, data = dat_diff))
ggsave("output/rep_gap.png", width = 2, height = 2)

#### heterozygosity and primary assembly size
log.model.dat.before <- fit_log(dat[dat$Purging == "Before",], dat[dat$Purging == "Before",]$Heterozygosity, dat[dat$Purging == "Before",]$Primary)
log.model.dat.after <- fit_log(dat[dat$Purging == "After",], dat[dat$Purging == "After",]$Heterozygosity, dat[dat$Purging == "After",]$Primary)

dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$Primary,
                       After = dat[dat$Purging == "After",]$Primary)
p <- plot_dots(dat_diff, "Heterozygosity (%)", "Primary Size (%)", 0.0, 1.8, 70, 110)
p + geom_hline(yintercept=100, linetype="dashed", color = "gray")
ggsave("output/het_prim.png", width = 2, height = 2)
fit=lm(Primary ~ Heterozygosity, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(Primary ~ Heterozygosity, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and k-mer duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$DupKmer,
                       After = dat[dat$Purging == "After",]$DupKmer)
p <- plot_dots(dat_diff, "Heterozygosity (%)", "k-mer dup. (%)", 0.0, 1.8, 0, 30)
p
ggsave("output/het_dupKmer.png", width = 2, height = 2)
fit=lm(DupKmer ~ Heterozygosity, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupKmer ~ Heterozygosity, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and BUSCO duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Heterozygosity,
                       Before = dat[dat$Purging == "Before",]$DupBUSCO,
                       After = dat[dat$Purging == "After",]$DupBUSCO)
p <- plot_dots(dat_diff, "Heterozygosity (%)", "BUSCO dup. (%)", 0.0, 1.8, 0, 20)
p
ggsave("output/het_dupBUSCO.png", width = 2, height = 2)
fit=lm(DupBUSCO ~ Heterozygosity, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupBUSCO ~ Heterozygosity, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and k-mer duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       Before = dat[dat$Purging == "Before",]$DupKmer,
                       After = dat[dat$Purging == "After",]$DupKmer)
p <- plot_dots(dat_diff, "Repeat (%)", "k-mer dup. (%)", 0, 60, 0, 30)
p
ggsave("output/rep_dupKmer.png", width = 2, height = 2)
fit=lm(DupKmer ~ Repeat, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupKmer ~ Repeat, data = dat[dat$Purging == "After",])
summary(fit)

#### heterozygosity and BUSCO duplication
dat_diff <- data.frame(x = dat[dat$Purging == "Before",]$Repeat,
                       Before = dat[dat$Purging == "Before",]$DupBUSCO,
                       After = dat[dat$Purging == "After",]$DupBUSCO)
p <- plot_dots(dat_diff, "Repeat (%)", "BUSCO dup. (%)", 0, 60, 0, 20)
p
ggsave("output/rep_dupBUSCO.png", width = 2, height = 2)
fit=lm(DupBUSCO ~ Repeat, data = dat[dat$Purging == "Before",])
summary(fit)
fit=lm(DupBUSCO ~ Repeat, data = dat[dat$Purging == "After",])
summary(fit)






#### Heterozygosity and alternate assembly size
ggplot(dat, aes(x = Heterozygosity, y = Alts)) +
  geom_hline(yintercept=100, linetype="dashed", color = "grey") +
  geom_smooth(method='lm', formula= (y ~ log(x)), se=F, color = "grey") +
  geom_point() + theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits=c(0,120)) +
  scale_x_continuous(expand = c(0, 0), limits=c(0.00,1.80)) +  # for both axes to remove unneeded padding
  xlab("Heterozygosity (%)") + ylab("Alt. assembly size (%)")
ggsave("het_alt.png", width = 2, height = 2)

fit=lm(Alts ~ log(Heterozygosity), data = dat)
summary(fit)


