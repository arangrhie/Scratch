library(ggplot2)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)

getwd()
setwd("T2T-Polishing")

plot_points_by_qual <- function(dat, title="", QUAL, category, xmax=0, ymax=0) {
    if (xmax==0) {
        xmax = max(dat$REF_CNT)
    }
    if (ymax==0) {
        ymax = max(dat$ALT_CNT)
    }
    p <- ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, color=QUAL)) +
        theme_classic() +
        geom_point(alpha=0.3, size = 0.2) +
        scale_colour_gradient(low = "white", high = "black") +
        ggtitle(title) +
        scale_x_continuous(labels = comma, limits = c(0, xmax)) +
        scale_y_continuous(labels = comma, limits = c(0, ymax)) +
        labs(color = category) +
        xlab("Ref. Count") + ylab("Alt. Count") +
        theme(
            plot.title = element_text(face="bold", size=7),
            axis.title = element_text(face="bold", size=6),
            axis.text = element_text(face="bold", size=5),
            legend.position = c(0.95,0.80),
            legend.key = element_blank(),
            legend.text = element_text(size = 5),
            legend.title = element_text(face="bold", size = 5),
            legend.key.size = unit(2, "mm"),
            legend.margin=margin(c(0,0,0,0)))
    return(p)
}

# v0.9
dat=read.table("input/ED_Fig4/ED_Fig4a_v0.9_combined.cnt", header=F)
names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GL", "GT", "REF_CNT", "ALT_CNT")
dat_filt=dat[dat$REF_CNT < 200.0 & dat$ALT_CNT < 200.0 & dat$GT == "HOM", ]
summary(dat_filt)

p1 = plot_points_by_qual(dat_filt, "v0.9", dat$GQ, "GQ", xmax = 200, ymax = 200)
p1

# v1.0
dat=read.table("input/ED_Fig4/ED_Fig4a_v1.0_combined.cnt", header=F)
names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GL", "GT", "REF_CNT", "ALT_CNT")
dat_filt=dat[dat$REF_CNT < 200.0 & dat$ALT_CNT < 200.0 & dat$GT == "HOM", ]
summary(dat_filt)

p2 = plot_points_by_qual(dat_filt, "v1.0", dat$GQ, "GQ", xmax = 200, ymax = 200)
p2

plot_grid(p1, p2, nrow = 1)
ggsave(file = "output/ED_Fig4a_hom_gq.200x.jpg", width = 2.8, height = 1.5)




plot_points <- function(dat, title="", xmax=0, ymax=0) {
    if (xmax==0) {
        xmax = max(dat$REF_CNT)
    }
    if (ymax==0) {
        ymax = max(dat$ALT_CNT)
    }
    p <- ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, fill=GT, color=GT)) +
        theme_classic() +
        geom_point(shape=21, alpha=0.3) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        scale_x_continuous(labels = comma, limits = c(0, xmax)) +
        scale_y_continuous(labels = comma, limits = c(0, ymax)) +
        ggtitle(title) +
        theme(
            plot.title = element_text(face="bold", size=7),
            axis.title = element_text(face="bold", size=7),
            axis.text  = element_text(face="bold", size=6),
            legend.position = "none",
            legend.key = element_blank())
    return(p)
}

# Get legends
legend_gt <- function(dat) {
    p<-ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, fill=GT, color=GT)) +
        theme_classic() +
        geom_point(shape=21, alpha=0.3) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        theme(
            legend.position = "left",
            legend.key = element_blank(),
            legend.text = element_text(size = 6),
            legend.title = element_text(face="bold", size = 7))
    return(get_legend(p))
}

legend_qual <- function(dat, QUAL, category) {
    p<-ggplot(dat, mapping = aes(x=REF_CNT, y=ALT_CNT, color=QUAL)) +
        theme_classic() +
        geom_point(alpha=0.3) +
        scale_colour_gradient(low = "white", high = "black") +
        #scale_color_gradient2(midpoint=mid, low="blue", mid="yellow",
        #                      high="red", space ="Lab" ) +
        labs(color = category) +
        theme(
            legend.position = "left",
            legend.key = element_blank(),
            legend.text = element_text(size = 6),
            legend.title = element_text(face="bold", size = 7))
    return(get_legend(p))
}

# Legends
legend1=legend_gt(chr_only.noMulti)
legend2=legend_qual(data_sum, mid_qual, data_sum$QUAL, "QUAL")
legend3=legend_qual(data_sum, mid_gq, data_sum$GQ, "GQ")

chr_only.noMulti=dat[dat$GT == "HOM",]
head(chr_only.noMulti)
chr_only.noMulti$ALT_CNT=as.numeric(chr_only.noMulti$ALT_CNT)

# Non-Missing
dat=read.table("input/0918/0918.hybrid.PASS.VAF0.5_GQ30.cnt", header=F)
dat=read.table("input/0918/0918.combined.cnt", header=F)

names(dat) <- c("Chr", "Pos", "Ref", "Alt", "QUAL", "FILTER", "GQ", "GL", "GT", "REF_CNT", "ALT_CNT")
summary(dat$GT)
non_chr_only.noMulti=dat[dat$GT == "HOM",]
head(non_chr_only.noMulti)
non_chr_only.noMulti$ALT_CNT=as.numeric(non_chr_only.noMulti$ALT_CNT)

# Legends
data_sum <- rbind(chr_only.noMulti, non_chr_only.noMulti)
mid_qual<-mean(data_sum$QUAL)
mid_gq<-mean(data_sum$GQ)
mid_qual
mid_gq
legend1=legend_gt(chr_only.noMulti)
legend2=legend_qual(data_sum, mid_qual, data_sum$QUAL, "QUAL")
legend3=legend_qual(data_sum, mid_gq, data_sum$GQ, "GQ")

# Subset to max 200x
SUBSET=200
chr_only.noMulti.sub=chr_only.noMulti[chr_only.noMulti$REF_CNT < SUBSET & chr_only.noMulti$ALT_CNT < SUBSET,]
non_chr_only.noMulti.sub=non_chr_only.noMulti[non_chr_only.noMulti$REF_CNT < SUBSET & non_chr_only.noMulti$ALT_CNT < SUBSET,]

p1 = plot_points_by_qual(chr_only.noMulti.sub, "0904", mid_gq, chr_only.noMulti.sub$GQ, xmax = 200, ymax = 200)
p2 = plot_points_by_qual(non_chr_only.noMulti.sub, "0918", mid_gq, non_chr_only.noMulti.sub$GQ, xmax = 200, ymax = 200)
plot_grid(legend3, p1, p2, nrow = 1, rel_widths = c(0.1, 0.45, 0.45))
ggsave(file = "output/hybrid_0904_0918.hom.200x.png", width = 4, height = 2)
ggsave(file = "output/combined_0904_0918.hom.200x.png", width = 6, height = 3)


p1 = plot_points(chr_only.noMulti, "0904")
p2 = plot_points(non_chr_only.noMulti, "0918")
p3 = plot_points_by_qual(chr_only.noMulti, "0904", mid_qual, chr_only.noMulti$QUAL)
p4 = plot_points_by_qual(non_chr_only.noMulti, "0918", mid_qual, non_chr_only.noMulti$QUAL)
p5 = plot_points_by_qual(chr_only.noMulti, "0904", mid_gq, chr_only.noMulti$GQ)
p6 = plot_points_by_qual(non_chr_only.noMulti, "0918", mid_gq, non_chr_only.noMulti$GQ)
plot_grid(legend1, p1, p2, legend2, p3, p4, legend3, p5, p6, nrow = 3, rel_widths = c(0.1, 0.45, 0.45))
ggsave(file = "output/hybrid_0904_0918.png", width = 6, height = 8)
ggsave(file = "output/combined_0904_0918.png", width = 6, height = 8)

# Subset to max 200x
SUBSET=200
chr_only.noMulti.sub=chr_only.noMulti[chr_only.noMulti$REF_CNT < SUBSET & chr_only.noMulti$ALT_CNT < SUBSET,]
non_chr_only.noMulti.sub=non_chr_only.noMulti[non_chr_only.noMulti$REF_CNT < SUBSET & non_chr_only.noMulti$ALT_CNT < SUBSET,]

p1 = plot_points(chr_only.noMulti.sub, "0904", xmax = 200, ymax = 200)
p2 = plot_points(non_chr_only.noMulti.sub, "0918", xmax = 200, ymax = 200)
p3 = plot_points_by_qual(chr_only.noMulti.sub, "0904", mid_gq, chr_only.noMulti$GQ, xmax = 200, ymax = 200)
p4 = plot_points_by_qual(non_chr_only.noMulti.sub, "0918", mid_gq, non_chr_only.noMulti$GQ, xmax = 200, ymax = 200)
plot_grid(legend1, p1, p2, legend3, p3, p4, nrow = 2, rel_widths = c(0.1, 0.45, 0.45))
ggsave(file = "output/hybrid_0904_0918.hom.200x.png", width = 6, height = 6)


p3 = plot_points_by_qual(chr_only.noMulti.sub, "0904", mid_qual, chr_only.noMulti$QUAL, xmax = 200, ymax = 200)
p4 = plot_points_by_qual(non_chr_only.noMulti.sub, "0918", mid_qual, non_chr_only.noMulti$QUAL, xmax = 200, ymax = 200)
p5 = plot_points_by_qual(chr_only.noMulti.sub, "0904", mid_gq, chr_only.noMulti$GQ, xmax = 200, ymax = 200)
p6 = plot_points_by_qual(non_chr_only.noMulti.sub, "0918", mid_gq, non_chr_only.noMulti$GQ, xmax = 200, ymax = 200)
plot_grid(legend1, p1, p2, legend2, p3, p4, legend3, p5, p6, nrow = 3, rel_widths = c(0.1, 0.45, 0.45))
ggsave(file = "output/hybrid_0904_0918.200x.png", width = 6, height = 8)







p1 = plot_points(chr_only.noMulti, "Merqury \"Missing\"")
p2 = plot_points(non_chr_only.noMulti, "\"Non-Missing\"")
p3 = plot_points(non_chr_only.noMulti.sub, "\"Non-Missing, <200x\"", xmax = SUBSET, ymax = SUBSET)

p4 = plot_points_by_qual(chr_only.noMulti, "Merqury \"Missing\"", mid_qual, chr_only.noMulti$QUAL)
p5 = plot_points_by_qual(non_chr_only.noMulti, "\"Non-Missing\"", mid_qual, non_chr_only.noMulti$QUAL)
p6 = plot_points_by_qual(non_chr_only.noMulti.sub, "\"Non-Missing, <200x\"", mid_qual, non_chr_only.noMulti.sub$QUAL, xmax = SUBSET, ymax = SUBSET)

p7 = plot_points_by_qual(chr_only.noMulti, "Merqury \"Missing\"", mid_gq, chr_only.noMulti$GQ)
p8 = plot_points_by_qual(non_chr_only.noMulti, "\"Non-Missing\"", mid_gq, non_chr_only.noMulti$GQ)
p9 = plot_points_by_qual(non_chr_only.noMulti.sub, "\"Non-Missing, <200x\"", mid_gq, non_chr_only.noMulti.sub$GQ, xmax = SUBSET, ymax = SUBSET)

plot_grid(legend1, p1, p2, p3, legend2, p4, p5, p6, legend3, p7, p8, p9, nrow = 3, rel_widths = c(0.1, 0.3, 0.3, 0.3))
#ggsave(file = "output/10x_deepvariant_pass_allele_count.png", width = 8, height = 8)
#ggsave(file = "output/10x_deepvariant_pass_realign.png", width = 8, height = 8)
ggsave(file = "output/PCR-Free.png", width = 8, height = 8)

# QUAL > mid_qual
mid_qual=30
chr_only.noMulti.qual=chr_only.noMulti[chr_only.noMulti$QUAL>mid_qual,]
non_chr_only.noMulti.qual=non_chr_only.noMulti[non_chr_only.noMulti$QUAL>mid_qual,]
non_chr_only.noMulti.sub.qual=non_chr_only.noMulti.sub[non_chr_only.noMulti.sub$QUAL>mid_qual,]

# GQ > mid_gq
mid_gq=30
chr_only.noMulti.gq=chr_only.noMulti[chr_only.noMulti$GQ>mid_gq,]
non_chr_only.noMulti.gq=non_chr_only.noMulti[non_chr_only.noMulti$GQ>mid_gq,]
non_chr_only.noMulti.sub.gq=non_chr_only.noMulti.sub[non_chr_only.noMulti.sub$GQ>mid_gq,]

# QUAL > mid_qual & GQ > mid_gq
summary(chr_only.noMulti)
summary(chr_only.noMulti.qual)
summary(chr_only.noMulti.gq)

summary(non_chr_only.noMulti)
summary(non_chr_only.noMulti.qual)
summary(non_chr_only.noMulti.gq)

p1 = plot_points_by_qual(chr_only.noMulti.qual, "Merqury \"Missing\"", mid_qual, chr_only.noMulti.qual$QUAL)
p2 = plot_points_by_qual(non_chr_only.noMulti.qual, "\"Non-Missing\"", mid_qual, non_chr_only.noMulti.qual$QUAL)
p3 = plot_points_by_qual(non_chr_only.noMulti.sub.qual, "\"Non-Missing, <200x\"", mid_qual, non_chr_only.noMulti.sub.qual$QUAL, xmax = SUBSET, ymax = SUBSET)

p4 = plot_points_by_qual(chr_only.noMulti.gq, "Merqury \"Missing\"", mid_gq, chr_only.noMulti.gq$GQ)
p5 = plot_points_by_qual(non_chr_only.noMulti.gq, "\"Non-Missing\"", mid_gq, non_chr_only.noMulti.gq$GQ)
p6 = plot_points_by_qual(non_chr_only.noMulti.sub.gq, "\"Non-Missing, <200x\"", mid_gq, non_chr_only.noMulti.sub.gq$GQ, xmax = SUBSET, ymax = SUBSET)

plot_grid(legend2, p1, p2, p3, legend3, p4, p5, p6, nrow = 2, rel_widths = c(0.1, 0.3, 0.3, 0.3))
#ggsave(file = "output/10x_deepvariant_pass_allele_count.filt.30.png", width = 8, height = 5)
ggsave(file = "output/10x_deepvariant_pass_realign.filt.30.png", width = 8, height = 5)
ggsave(file = "output/PCR-Free.30.png", width = 8, height = 5)
