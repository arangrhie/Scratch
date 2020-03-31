#!/usr/bin/env Rscript

require("argparse")
require("ggplot2")
require("scales")

parser <- ArgumentParser(description = "Make blob plots.")
parser$add_argument("-f", "--file", type="character", help=".count tdf; with headers; ie. <category> <seqId> <hap1Count> <hap2Count> <seqSize> (required)", default=NULL)
parser$add_argument("-o", "--output", type="character", default="hapmers.blob", help="output prefix [default %(default)s]")
parser$add_argument("-t", "--type", type="character", default="exact", help="available types: exact or scaled. [default %(default)s]" )
parser$add_argument("-m", "--ymax", type="double", default=0, help="maximum Y scale, auto detect max if 0 [default %(default)s]")
parser$add_argument("-x", "--xdim", type="double", default=6.5, help="width of output plot [default %(default)s]")
parser$add_argument("-y", "--ydim", type="double", default=6, help="height of output plot [default %(default)s]")
parser$add_argument("-p", "--pdf", dest='pdf', default=FALSE, action='store_true', help="set to get output in .pdf. [default .png]")
args <- parser$parse_args()

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

plot_blob_exact <- function(dat=NULL, ymax) {
  
  if (ymax == 0) {
    max_total = max(max(dat[,3]), max(dat[,4])) * 1.1
  } else {
    max_total = ymax
  }
  col_lab=names(dat)[1]
  x_lab=names(dat)[3]
  y_lab=names(dat)[4]
  print(x_lab)
  print(y_lab)
  
  p <- ggplot(dat, aes(x=dat[,3], y=dat[,4], fill=dat[,1], size=dat[,5])) +
    geom_point(shape=21, alpha=0.3) + theme_bw() +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +     # Set1 -> Red / Blue. Set2 -> Green / Orange.
    scale_x_continuous(labels = comma, limits = c(0, max_total)) +
    scale_y_continuous(labels = comma, limits = c(0, max_total)) +
    scale_size_continuous(labels = comma, range = c(1, 10), name = "Total k-mers (size)") +
    theme(legend.text = element_text(size=12),
          #legend.position = c(0.95,0.95),  # Modify this if the legend is covering your favorite circle
          legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
          #legend.box.just = "right",
          #legend.justification = c("right", "top"),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=11)) +
    guides( size = guide_legend(order = 1),
            fill = guide_legend(override.aes = list(size=5, alpha=1), order = 2, title = col_lab)) +
    xlab(x_lab) + ylab(y_lab)
  
  p
  return(p)
}

plot_blob_scaled <- function(dat=NULL, ymax = 0, plot_log=FALSE) {
  head(dat)
  dat$Scld=ifelse(dat[,3] == dat[,4], 0.5, 100*dat[,4]/(dat[,3]+dat[,4]))
  
  if (ymax == 0) {
    max_total = max(dat[,5]) * 1.1
  } else {
    max_total = ymax
  }
  col_lab=names(dat)[1]
  x_lab=paste(names(dat)[4] , "(%)")
  y_lab=names(dat)[5]
  print(max_total)
  print(x_lab)
  print(y_lab)
  
  p <- ggplot(dat, aes(x=Scld, y=Size, fill=dat[,1], size=dat[,5])) +
    geom_point(shape=21, alpha=0.3) + theme_bw() +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")

  if (!plot_log) {
    # Plot in standard continuous Y
    p <- p + 
      scale_y_continuous(labels=comma, limits = c(0, max_total))
  } else {
    # Plot in log scaled Y
    p <- p + 
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits = c(1, max_total)) +
      annotation_logticks(sides = "lr")
  }
  p <- p +
    scale_size_continuous(labels=comma, range = c(1, 10), name = "Total k-mers (size)") +
    theme(legend.text = element_text(size=11),
          legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12)) +
    guides( size = guide_legend(order = 1),
            fill = guide_legend(override.aes = list(size=5, alpha=1), order = 2, title = col_lab)) +
    xlab(x_lab) + ylab(y_lab)
  
  return(p)
}



blob_plot <- function(dat=NULL, out, ymax = 0, w=6.5, h=6, pdf=FALSE, type="exact", plot_log=FALSE) {

  dat=read.table(dat, header=TRUE)
  
  if (type == "exact") {
    print("Plot exact k-mer counts")
    p <- plot_blob_exact(dat, ymax)
  } else {
    print("plot scaled k-mer counts")
    w=8
    p <- plot_blob_scaled(dat, ymax, plot_log)
  }
  p
  
  outformat="png"
  if (pdf) {
    outformat="pdf"
  }  
  ggsave(file = paste(out, outformat, sep="."), height = h, width = w)

}

blob_plot(dat = args$file, out = args$output, ymax = args$ymax, w = args$xdim, h = args$ydim, pdf = args$pdf, type=args$type, plot_log = args$plot_log)

