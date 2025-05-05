library(ggplot2)
library(scales)

getwd()
setwd("T2T-CenSat")

bsat=rgb(250, 153, 255, 1, maxColorValue=256)
censat=rgb(0, 204, 204, 1, maxColorValue=256)
ct=rgb(224, 224, 224, 1, maxColorValue=256)
dhor=rgb(153,0,0,1, maxColorValue=256)
gsat=rgb(0,153,76,1, maxColorValue=256)
hor=rgb(250,0,0,1, maxColorValue=256)
hsat1=rgb(0,0,153, 1, maxColorValue=256)
hsat2=rgb(0,128,250, 1, maxColorValue=256)
hsat3=rgb(0,0,250, 1, maxColorValue=256)
hsat4=rgb(0,51,102, 1, maxColorValue=256)
hsat5=rgb(102,128,255, 1, maxColorValue=256)
mon=rgb(255,204,153, 1, maxColorValue=256)
no_cen=rgb(102,102,102,1, maxColorValue = 256)
rDNA=rgb(250,0,250, 1, maxColorValue=256)
tar=rgb(0,0,0, 1, maxColorValue=256)

cenPalette=c(bsat, censat, ct, dhor, gsat, hor, hsat1, hsat2, hsat3, hsat4, hsat5, mon, no_cen, rDNA, tar)

stats <- function(dat) {
    max(dat$Spacing)
    centypes = unique(dat$Category)
    centypes
    for (cenType in centypes) {
        print(cenType)
        print(summary(dat[dat$Category == cenType,]$Spacing))
    }
}

plot_density_log <- function(dat) {
    # Fix order
    dat$Category=factor(dat$Category, levels=levels(dat$Category))
    ggplot(data = dat, aes(x=Spacing, color = Category, fill=Category)) + 
        scale_color_manual(values = cenPalette) +
        scale_fill_manual(values = cenPalette) +
        geom_density(alpha=0.4, color="black", size = 0.1) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)),
                      limits = c(1, 1000000)) +
        theme_bw() +
        annotation_logticks() +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              legend.key.size = unit(4, "mm"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0))) +
        xlab("Spacing (bp)") + ylab("Density")

}

plot_violin <- function(dat) {
    ggplot(data = dat, aes(x=Category, y=Spacing, fill=Category)) +
        geom_violin(alpha=0.6, size=0.1) +
        geom_boxplot(alpha=0.4, size=0.1, width=0.1, outlier.size = 0.1, outlier.alpha = 0.4) +
        #scale_fill_manual(values = cenPalette) +
        theme_bw() +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)),
                      limits = c(1, 1000000)) +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              legend.key.size = unit(4, "mm"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0))) +
        ylab("Spacing (bp)")
}

plot_violin_by_k <- function(dat) {
    dodge <- position_dodge(width = 0.7)
    ggplot(data = dat, aes(x=Category, y=Spacing, fill=k)) +
        geom_violin(alpha=0.6, size=0.1, position = dodge) +
        geom_boxplot(alpha=0.4, size=0.1, width=0.1, outlier.size = 0.1, outlier.alpha = 0.4, position = dodge) +
        scale_fill_brewer(palette="Dark2") +
        theme_bw() +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)),
                      limits = c(1, 1000000)) +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              legend.key.size = unit(4, "mm"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0))) +
        ylab("Spacing (bp)")
}

## marker distance for all cenHap annotations
dat=read.table("input/k100.spacing.all", header=F)
names(dat) <- c("Spacing", "Category")
stats(dat)
plot_violin(dat)
ggsave("output/chm13.v1.cen_v2.marker_density.k100.violin.pdf", width = 4.5, height = 3, device=cairo_pdf)

dat=read.table("input/k51.spacing.all", header=F)
names(dat) <- c("Spacing", "Category")
stats(dat)
#plot_density_log(dat)
#ggsave("output/chm13.v1.cen_v2.marker_density.k51.pdf", width = 4, height = 3, device=cairo_pdf)
plot_violin(dat)
ggsave("output/chm13.v1.cen_v2.marker_density.k51.violin.pdf", width = 4.5, height = 3, device=cairo_pdf)

dat=read.table(gzfile("input/k21.spacing.all.gz"), header=F)
names(dat) <- c("Spacing", "Category")
stats(dat)
#plot_density_log(dat)
#ggsave("output/chm13.v1.cen_v2.marker_density.k21.pdf", width = 4, height = 3, device=cairo_pdf)
plot_violin(dat)
ggsave("output/chm13.v1.cen_v2.marker_density.k21.violin.pdf", width = 4.5, height = 3, device=cairo_pdf)

## marker distance, k=21, ct_*(*_arm) moved to no_cen
dat=read.table(gzfile("input/k21.spacing.gz"), header=F)
names(dat) <- c("Spacing", "Category")
stats(dat)
plot_violin(dat)
ggsave("output/chm13.v1.cen_v2.marker_density.k21.violin.no_ct_arm.pdf", width = 4.5, height = 3, device=cairo_pdf)

## HOR and HSat1-3 Spacing
dat=read.table(gzfile("input/hor_and_hsat.spacing.gz"), header=F)
names(dat) <- c("Spacing", "Category", "k")
stats(dat)
dat$k = factor(dat$k)
plot_violin_by_k(dat)
ggsave("output/chm13.v1.cen_v2.marker_density.hor_and_hsat.pdf", width = 4, height = 2, device=cairo_pdf)

## HOR spacing per chromosomes
dat=read.table("input/k21.hor.spacing")
names(dat) <- c("Spacing", "Category")
stats(dat)
dat$Category = factor(dat$Category, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                                             "chr21", "chr22", "chrX"))
plot_violin(dat)
ggsave("output/chm13.v1.cen_v2.marker_density.hor_per_chr.pdf", width = 5.4, height = 3, device=cairo_pdf)

dat=read.table("input/hor.stat", header=F)
names(dat) <- c("Chr", "k", "Num", "Total", "Min", "Avg", "N50", "Max")
head(dat)
dat$k = factor(dat$k, levels = c("k21", "k51", "k100"))
ggplot(dat, aes(x = Avg/1000, y = Max/1000, color = k)) +
    geom_point(data = dat, mapping = aes(group=interaction(Chr,k), shape = k), size=0.6) +
    scale_shape_manual(values = c(2,3,4), name = "k") +
    theme_bw() +
    theme(axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=5),
          legend.key.size = unit(4, "mm"),
          legend.text = element_text(size=5),
          legend.title = element_text(size=5),
          legend.margin=margin(c(0,0,0,0))) +
    xlab("Avg. (kbp)") + ylab("Max. (kbp)")
ggsave("output/chm13.v1.cen_v2.marker_density.hor_stat.pdf", width = 3, height = 2, device=cairo_pdf)

