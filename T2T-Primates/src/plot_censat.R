library(ggplot2)
library(scales)

getwd()
setwd("Acro/src")

# Load data
# Input
# Sample <tab> chrom <tab> Hap <tab> CenSat <tab> Size (in Mbp or bp)
dat1 = read.table("../input/all.chrom_breakdown.clean.hap1.tsv", sep = "\t", header = TRUE)
summary(dat1)

dat2 = read.table("../input/all.chrom_breakdown.clean.hap2.tsv", sep = "\t", header = TRUE)
summary(dat2)

dat3 = read.table("../input/censat_summary.clean.tsv", sep = "\t", header = TRUE)
summary(dat3)


# Color palette
aSat=rgb(153,0,0, maxColorValue=255)
bSat=rgb(250,153,255, maxColorValue=255)
ct=rgb(224,224,224, maxColorValue=255)
HSat1A=rgb(0,222,96, maxColorValue=255)
HSat1B=rgb(27,153,139, maxColorValue=255)
HSat2=rgb(0,128,250, maxColorValue=255)
HSat3=rgb(51,81,137, maxColorValue=255)
Other=rgb(0,204,204, maxColorValue=255)
rDNA=rgb(102,47,144, maxColorValue=255)
SD=rgb(255,146,0, maxColorValue=255)

# Has to match the order of dat$CenSat
censat_col = c(ct, aSat, bSat, HSat1A, HSat1B, HSat2, HSat3, Other, rDNA, SD)

# Plot function
plot_censat <- function(dat, y_limit) {
  dat$CenSat = factor(dat$CenSat, levels = c("ct", "aSat", "bSat", "HSat1A", "HSat1B", "HSat2", "HSat3", "Other", "rDNA", "SD"))

  ggplot(dat, aes(x=Sample_Chromosome, y=Size, fill=CenSat)) + 
    geom_bar(position=position_stack(), stat="identity") +
    theme_classic() +
	theme(
	  axis.title = element_text(face="bold", size=6),
	  axis.text = element_text(face="bold", size=5),
	  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	  legend.position = "bottom",
	  legend.text = element_text(size=5),
	  legend.key.size = unit(3, "mm"),
	  legend.title = element_blank()) +
    scale_fill_manual(values = censat_col) +
    #scale_x_continuous(labels=Chromosome) +
    scale_y_continuous(expand = c(0, 0), label=comma, limits=c(0, y_limit)) +
	xlab("Chromosome") +
	guides(fill = guide_legend(nrow = 1))
}

# Plot hap1
plot_censat(dat = dat1, y_limit = 50000000)
ggsave("../output/acro_censat_hap1.png", width=5, height=3)
ggsave("../output/acro_censat_hap1.pdf", width=5, height=3)


# Plot hap2
plot_censat(dat = dat2, y_limit = 50000000)
ggsave("../output/acro_censat_hap2.png", width=5, height=3)
ggsave("../output/acro_censat_hap2.pdf", width=5, height=3)

# Plot summary
plot_censat(dat = dat3, y_limit = 370000000)
ggsave("../output/acro_censat_summary.png", width=3, height=3)
ggsave("../output/acro_censat_summary.pdf", width=3, height=3)
