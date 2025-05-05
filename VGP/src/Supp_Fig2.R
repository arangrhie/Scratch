require("ggplot2")
require("scales")
library(cowplot)

gray = "black"
red = "#E41A1C"
blue = "#377EB8" # light blue = "#56B4E9"
green = "#4DAF4A"
purple = "#984EA3"  # purple = "#CC79A7"
orange = "#FF7F00"  # orange = "#E69F00"
yellow = "#FFFF33"

ALPHA=0.4

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

merqury_col = c(gray, red, blue, green, purple, orange)

merqury_brw <- function(dat, direction=1) {
    merqury_colors=merqury_col[1:length(unique(dat))]
    if (direction == -1) {
        merqury_colors=rev(merqury_colors)
    }
    if (length(unique(dat)) == 1) {
        merqury_colors = c(red)
    }
    merqury_colors
}

format_theme <- function() {
    theme(#legend.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5, size = 6),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text=element_text(size=6))
}

plot_zero_line <- function(zero) {
    if (!is.null(zero)) {
        if (length(zero[,1]) == 1) {
            scale_fill_manual(values = c(red), name="k-mer")
        } else if (length(zero[,1]) == 2) {
            scale_fill_manual(values = c(blue, red), name="k-mer")
        } else if (length(zero[,1]) == 3) {
            scale_fill_manual(values = c(purple, blue, red), name="k-mer")
        } else {
            scale_fill_manual(values = merqury_brw(zero[,1]), name="k-mer")
        }
    }
}

plot_cutoff <- function(cutoff) {
    if (!is.null(cutoff)) {
        geom_vline(data = cutoff, aes(xintercept = cutoff[,2], colour = cutoff[,1]), show.legend = FALSE, linetype="dashed", size=LINE_SIZE)
    }
}


plot_zero_fill <- function(zero) {
    if (!is.null(zero)) {
        geom_bar(data=zero, aes(x=zero[,2], y=zero[,3], fill=zero[,1], colour=zero[,1], group=zero[,1]),
                 position="stack", stat="identity", show.legend = FALSE, width = 2, alpha=ALPHA)
    }
}

plot_zero_stack <- function(zero) {
    if (!is.null(zero)) {
        geom_bar(data=zero, aes(x=zero[,2], y=zero[,3], fill=zero[,1], colour=zero[,1]),
                 position="stack", stat="identity", show.legend = FALSE, width = 2, alpha=ALPHA)
    }
}

plot_line <- function(dat, name, x_max, y_max, zero, cutoff) {
    ggplot(data=dat, aes(x=kmer_multiplicity, y=Count, group=dat[,1], colour=dat[,1])) +
        geom_line() +
        ggtitle(name) +
        scale_color_manual(values = merqury_brw(dat[,1]), name="k-mer") +
        plot_zero_line(zero=zero) +
        plot_zero_fill(zero=zero) +
        plot_cutoff(cutoff) +
        theme_bw() +
        format_theme() +
        scale_y_continuous(labels=fancy_scientific) +
        coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max))
}

spectra_cn_plot  <-  function(hist, name="out", zero="", cutoff="", w=6, h=4.5, x_max=0, y_max=0, type="line", pdf=FALSE) {
    # Read hist
    dat=read.table(hist, header=TRUE)
    dat[,1]=factor(dat[,1], levels=unique(dat[,1]), ordered=TRUE) # Lock in the order
    
    # Read asm-only
    dat_0 = NULL
    if (zero != "") {
        dat_0=read.table(zero, header=FALSE)
        dat_0[,1]=as.factor(dat_0[,1])
        dat_0[,1]=factor(dat_0[,1], levels=rev(unique(dat_0[,1])), ordered=TRUE)
    }
    
    # Read cutoffs
    dat_cut = NULL
    if (cutoff != "") {
        dat_cut=read.table(cutoff, header = FALSE)
        dat_cut[,1]=as.factor(dat_cut[,1])
        dat_cut[,1]=factor(dat_cut[,1], levels=unique(dat_cut[,1]), ordered=TRUE)
    }
    
    # x and y max
    y_max_given=TRUE;
    if (y_max == 0) {
        y_max=max(dat[dat[,1]!="read-total" & dat[,1]!="read-only" & dat[,2] > 3,]$Count)
        y_max_given=FALSE;
    }
    if (x_max == 0) {
        x_max=dat[dat[,3]==y_max,]$kmer_multiplicity
        x_max=x_max*2.5
    }
    if (! y_max_given) {
        y_max=y_max*1.1
    }
    print(paste("x_max:", x_max, sep=" "))
    if (zero != "") {
        y_max=max(y_max, sum(dat_0[,3]*1.1))	# Check once more when dat_0 is available
    }
    print(paste("y_max:", y_max, sep=" "))
    
    outformat="png"
    if (pdf) {
        outformat="pdf"
    }
    
    if (type == "all" || type == "line") {
        print("## Line graph")
        p <- plot_line(dat, name, x_max, y_max, zero = dat_0, cutoff = dat_cut)
        #print(p)
        #save_plot(name=name, type="ln", outformat, h=h, w=w)
    }
    
    if (type == "all" || type == "fill") {
        print("## Area under the curve filled")
        p <- plot_fill(dat, name, x_max, y_max, zero = dat_0, cutoff = dat_cut)
        #print(p)
        #save_plot(name=name, type="fl", outformat, h=h, w=w)
    }
    
    if (type == "all" || type == "stack") {
        print("## Stacked")
        p <- plot_stack(dat, name, x_max, y_max, zero = dat_0, cutoff = dat_cut)
        #print(p)
        #save_plot(name=name, type="st", outformat, h=h, w=w)
    }
    p
}

getwd()
setwd("VGP")
genomes = read.table("input/Supp_Fig2/genomes.list.srt", header=F)
colnames(genomes) = c("Genome")
genomes
i=1
myplots = list()
for (g in genomes$Genome) {
    print(g)
    cn=paste("input/Supp_Fig2/", g, ".spectra-asm.hist", sep = "")
    z=paste("input/Supp_Fig2/", g, ".dist_only.hist", sep = "")
    p1 <- spectra_cn_plot(hist=cn, zero = z, type = "line", name = paste(g, " spectra-asm", sep = ""))
    if (g == "bTaeGut2_trio") {
        cn=paste("input/Supp_Fig2/", g, ".mat.spectra-cn.hist", sep = "")
        z=paste("input/Supp_Fig2/", g, ".mat.only.hist", sep = "")
        p2 <- spectra_cn_plot(hist = cn, zero = z, type = "line", name = paste(g, " mat spectra-cn", sep = ""))
    } else {
        cn=paste("input/Supp_Fig2/", g, ".pri.spectra-cn.hist", sep = "")
        z=paste("input/Supp_Fig2/", g, ".pri.only.hist", sep = "")
        p2 <- spectra_cn_plot(hist = cn, zero = z, type = "line", name = paste(g, " pri spectra-cn", sep = ""))
    }
    myplots[[i]] <- p1  # add each plot into plot list
    i = i+1
    myplots[[i]] <- p2
    i = i+1
}

p1 <- plot_grid(plotlist = myplots, ncol = 4)
#ggsave(file = "output/pub/Supp_Fig2_merqury.pdf", width = 6.5, height = 7, p1)
ggsave(file = "output/pub/Supp_Fig2_merqury.png", width = 6.5, height = 7, p1)

