library(ggplot2)
library(reshape2)
library(dplyr)

setwd("T2T-Ref/")

plot_dist <- function(dat = null, out) {
    
    dat_DJ_PHR = dat %>% select(Sample, DJ, Chr13_PHR_arm1) # , Chr13_PHR_arm1, Chr13_PHR_arm2)
    
    # Melt the data
    data_melt <- melt(dat_DJ_PHR, id.vars = "Sample", variable.name = "Category", value.name = "Copies")
    head(data_melt)
    max(data_melt$Copies)
    
    data_melt$Category <- factor(data_melt$Category, levels = c("DJ", "Chr13_PHR_arm1")) #, "Chr13_PHR_arm1", "Chr13_PHR_arm2"))
    
    ggplot(data = data_melt, aes(x = Category, y = Copies, color = Category)) +
        #geom_boxplot(width=0.3, outlier.shape = 1) +
        #geom_boxplot(width=0.05, outlier.shape = NA) +
        geom_jitter(shape=1, position=position_jitter(0.4), alpha = 0.3, size = 0.3) +
        geom_violin(trim = F, scale = "width", alpha = 0.5, linewidth = 0.2) +
        scale_color_brewer(palette="Set1") +
        scale_y_continuous(breaks = seq(0, 14, 1), limits = c(0, 14)) +
        theme_bw() +
        theme(legend.position = "top",
              legend.justification = "right",
              legend.box = "horizontal")
    
    ggsave(paste("output/", out, ".png", sep = "" ), height = 3.5, width = 3)
}

dat=read.table("input/inova.copynum_dip.4187.txt", header = T)
head(dat)
plot_dist(dat, "count_dist_DJ_PHR")

# awk 'NR==1 || $NF>3.2 && $NF<4.2' inova.copynum_dip.txt > inova.copynum_dip.PHR_4.txt
dat=read.table("input/inova.copynum_dip.4187.PHR_4.txt", header = T)
plot_dist(dat, "count_dist_DJ_PHR_4")

dat=read.table("input/inova.copynum_dip.4187.DJ_9.txt", header = T)
plot_dist(dat, "count_dist_DJ9_PHR")

dat=read.table("input/inova.copynum_dip.4187.rob.txt", header = T)
plot_dist(dat, "count_dist_DJ_PHR_ROB")



# Plot rDNA CN count
max(dat$rDNA)

ggplot(data = dat, aes(x = "rDNA", y = rDNA, color = "rDNA")) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.05, outlier.shape = NA) +
    geom_jitter(shape=1, position=position_jitter(0.3)) +
    scale_color_brewer(palette="Set1") +
    ylim(0, 700)






##### Backup ######


dat=read.table("input/count_dist.tsv", header = T)
head(dat)

dat$Category <- factor(dat$Category, levels = c("DJ", "PHR_keep", "chr13_PHR_arm1", "chr13_PHR_arm2"))
dat$Copies = dat$Count * 2
ggplot(data = dat, aes(x = Category, y = Copies, color = Category)) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.05) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    scale_color_brewer(palette="Set1") +
    scale_y_continuous(breaks = seq(0, 12, 1), limits = c(0.5, 12))

ggsave("output/count_dist.png", height = 4, width = 7)

dat=read.table("input/count_dist_rdna.tsv", header = T)
head(dat)
dat$Copies = dat$Count * 2
ggplot(data = dat, aes(x = rDNA, y = Copies, color = rDNA)) +
    geom_violin(trim = F) +
    geom_boxplot(width=0.05) +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    scale_color_brewer(palette="Set1") +
    ylim(0, 700)

ggsave("output/count_dist_rdna.png", height = 4, width = 4)