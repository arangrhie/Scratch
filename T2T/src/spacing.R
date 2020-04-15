library(ggplot2)
library(scales)

getwd()
setwd("../T2T")
greenPalette=c("#CCFFCC", "#009E50")
dat=read.table("input/markers.dxz1.dat", header=T)
head(dat)
max(dat$Spacing)
summary(dat[dat$Category == "DXZ1",]$Spacing)

#levels=levels(dat$Category)[c(2,1)]
dat$Category=factor(dat$Category, levels=levels(dat$Category)[c(1,2)])
ggplot(data = dat, aes(x=Spacing, fill=Category)) + 
  scale_fill_manual(values = greenPalette) +
  geom_density(alpha=0.7) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  annotation_logticks() +
  xlab("Spacing (bp)") + ylab("Density")
ggsave("output/relv0.7_marker_density.pdf", width = 4.5, height = 3.8)


### Local paths
# setwd("/Users/rhiea/CHM13/chrX_markers/relv0.7")  # release v0.7
# setwd("/Users/rhiea/CHM13/chrX_markers")  # release v0.6
# setwd("/Users/rhiea/perfect-polish/kmer-spacing") # polishing benchmarks
