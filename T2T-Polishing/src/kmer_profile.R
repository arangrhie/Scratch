library(ggplot2)
library(scales)
library(data.table)
library(cowplot)
# install.packages('bit64')
# library("bit64")

getwd()
setwd("T2T-Polishing")

plot_heatmap <- function(dat = NULL, cnt, title = "", x_min = 0, x_max, y_max, pattern, counts) {
    p <- ggplot(dat, aes(x = V2, y = V1)) +
        geom_raster(aes(fill = cnt)) +
        theme_bw() +
        ggtitle(title) +
        labs(fill = counts) +
        scale_x_continuous("kmer multiplicity", expand = c(0, 0), limits=c(x_min, x_max)) +
        scale_y_continuous(paste(pattern, "in 21-mers", sep=" "), expand = c(0, 0), limits=c(0, y_max))
    return(p)
}

plot_heatmap_logx <- function(dat = NULL, cnt, title = "", x_min = 0, x_max, pattern) {
    p <- ggplot(dat, aes(x = V2, y = V1)) +
        geom_tile(aes(fill = cnt)) +
        theme_bw() +
        ggtitle(title) +
        labs(fill = "log(Count)") +
        scale_x_log10("kmer multiplicity",
                      breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_y_continuous(paste(pattern, "in 21-mers", sep=" "), expand = c(0, 0), limits=c(0, 11))
    return(p)
}

plot_by_multiplicity <- function(dat1, dat2, y_max, peak1, peak2, cp, pattern, outdir) {
    x_max1=peak1 * cp
    x_max2=peak2 * cp
    
    m1 <- plot_heatmap(dat1, dat1$LogCount,
                       paste("PacBio HiFi20k (", cp, "x)", sep=""),
                       x_min = 0, x_max = x_max1, y_max = y_max, pattern = pattern, counts = "log(Count)")
    m2 <- plot_heatmap(dat2, dat2$LogCount,
                       paste("Illumina PCR-Free (", cp, "x)", sep=""),
                       x_min = 0, x_max = x_max2, y_max = y_max, pattern = pattern, counts = "log(Count)")
    plot_grid(m1, m2, nrow = 1, rel_widths = c(0.5, 0.5))
    ggsave(paste("output/", outdir, "/", pattern, "_log_count_", cp, "x.png", sep=""),
           width = 8, height = 3)
    
    m1 <- plot_heatmap(dat1, dat1$V3,
                       paste("PacBio HiFi20k (", cp, "x)", sep=""),
                       x_min = 0, x_max = x_max1, y_max = y_max,
                       pattern = pattern, counts = "Count")
    
    m2 <- plot_heatmap(dat2, dat2$V3,
                       paste("Illumina PCR-Free (", cp, "x)", sep=""),
                       x_min = 0, x_max = x_max2, y_max = y_max,
                       pattern = pattern, counts = "Count")
    plot_grid(m1, m2, nrow = 1, rel_widths = c(0.5, 0.5))
    ggsave(paste("output/", outdir, "/", pattern, "_count_", cp, "x.png", sep=""),
           width = 8, height = 3)
}
# contour lines
# m + geom_density_2d()

# peak for normalizing copy numbers
peak1=31
peak2=105
y_max=21

########### 0915 ver. ############
dat1=fread("input/exclusive_0915/illumina.0.no_hifi_0.hifim.meryl.GA_TC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/exclusive_0915/hifi.0.no_illumina_0.illm.meryl.GA_TC.hist", header=F)
head(dat2)
# To handle "Don't know how to automatically pick scale for object of type integer64. Defaulting to continuous." error
# dat2$V3=dat2$V3/2
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GA", "exclusive_0915")
}

# Plot G/C of the missings kmer
dat1=fread("input/exclusive_0915/illumina.0.no_hifi_0.hifim.meryl.GC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/exclusive_0915/hifi.0.no_illumina_0.illm.meryl.GC.hist", header=F)
head(dat2)
# To handle "Don't know how to automatically pick scale for object of type integer64. Defaulting to continuous." error
# dat2$V3=dat2$V3/4
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)


for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GC", "exclusive_0915")
}



# Plot GA/TC dimer of the missings kmer
dat1=fread("input/missings/hifi.0.illm.meryl.GA_TC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/missings/illumina.0.hifim.meryl.GA_TC.hist", header=F)
head(dat2)
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

# y axis max
y_max=max(dat1$V1)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GA")
}

# Plot G/C of the missings kmer
dat2=fread("input/missings/hifi.0.illm.meryl.GC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/missings/illumina.0.hifim.meryl.GC.hist", header=F)
head(dat2)
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

y_max=max(dat1$V1)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GC")
}




########### 0904 ver. ############
# Plot GA/TC dimer of the missings kmer
dat1=fread("input/exclusive/hifi-only.meryl.GA_TC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/exclusive/illumina-only.meryl.GA_TC.hist", header=F)
head(dat2)
# To handle "Don't know how to automatically pick scale for object of type integer64. Defaulting to continuous." error
dat2$V3=dat2$V3/2
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

# y axis max
y_max=max(dat1$V1)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GA", "exclusive")
}

# Plot G/C of the missings kmer
dat2=fread("input/exclusive/hifi-only.meryl.GC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/exclusive/illumina-only.meryl.GC.hist", header=F)
head(dat2)
# To handle "Don't know how to automatically pick scale for object of type integer64. Defaulting to continuous." error
dat2$V3=dat2$V3/4
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

y_max=max(dat1$V1)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GC", "exclusive")
}



# Plot GA/TC dimer of the missings kmer
dat1=fread("input/missings/hifi.0.illm.meryl.GA_TC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/missings/illumina.0.hifim.meryl.GA_TC.hist", header=F)
head(dat2)
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

# y axis max
y_max=max(dat1$V1)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GA")
}

# Plot G/C of the missings kmer
dat2=fread("input/missings/hifi.0.illm.meryl.GC.hist", header=F)
dat1$LogCount=log(dat1$V3, base = 10)
summary(dat1)

dat2=fread("input/missings/illumina.0.hifim.meryl.GC.hist", header=F)
head(dat2)
dat2$LogCount=log(dat2$V3, base = 10)
summary(dat2)

y_max=max(dat1$V1)

for ( cp in c(3, 10, 50, 100, 500, 1000, 5000)) {
    plot_by_multiplicity(dat1, dat2, y_max, peak1, peak2, cp, "GC")
}


###################### Experimental ###########################
dat1=fread("input/kmer_profile/hifi20k.meryl.ga", header=F)
dat1=fread("input/kmer_profile/hifi20k.meryl.AG_CT.hist", header=F)
dat1=fread("input/kmer_profile/hifi20k.meryl.GC.hist", header=F)
head(dat1)
summary(dat1)
dat1$LogCount=log(dat1$V3, base = 10)
peak1=31
dat2=fread("input/kmer_profile/illumina.meryl.ga", header=F)
dat2=fread("input/kmer_profile/illumina.meryl.AG_CT.hist", header=F)
dat2=fread("input/kmer_profile/illumina.meryl.GC.hist", header=F)
head(dat2)
summary(dat2)
peak2=105
# Log count
dat2$LogCount=log(dat2$V3, base = 10)
# Absolute count, smooth if counts are too high - here the absolute counts aren't meaningful anymore.
dat2$Smooth=dat2$V3 / 2
y_max=max(dat1$V1)


cp = 2000
x_max=peak1 * cp
m1 <- plot_heatmap_logx(dat1, dat1$LogCount,
                   paste("PacBio HiFi20k (", cp, "x)", sep=""),
                   x_min = 0, x_max = x_max)
m1
x_max=peak2 * cp
m2 <- plot_heatmap_logx(dat2, dat2$LogCount,
                   paste("Illumina PCR-Free (", cp, "x)", sep=""),
                   x_min = 0, x_max = x_max)
plot_grid(m1, m2, nrow = 1, rel_widths = c(0.5, 0.5))
ggsave(paste("output/kmer_profile/GA_log_count_", cp, "x_logx.png", sep=""),
       width = 8, height = 3)
