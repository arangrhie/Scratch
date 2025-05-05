library(ggplot2)
library(directlabels)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)

#read in BED file
setwd("T2T-HG002Y")   # set to folder of your files
df = fread("input/AS-HOR-vs-chrY.track.bed", select = c(1:9), stringsAsFactors = TRUE, fill=TRUE, quote="", header=FALSE)
cols =  c("chr", "start", "stop", "region", "value", "strand", "start2", "stop2", "rgb")
colnames(df) <- cols
head(df)

#change all live array names
df$region <- gsub("S.*L.*", "Live", df$region)

#change all divergent array names
df$region <- gsub("\\bS.*d.*", "Dead", df$region)

#change all remaining SF names
df2=df
df$region <- gsub("\\bS.*", "Dead", df$region)

#reorder rows so that live arrays are plotted on top
df$region <- factor(df$region, levels = c("ct", "ALR/Alpha", "Live", "Dead", "BSR/Beta","GSAT", "GSATII", "GSATX", "HSat1A", "HSat1B", "HSat2", "HSat3"), ordered = T)

#centromere colors (HG00733 and NHPs)
myColors <- c("#FFCC99", "#D1ADCA",
              "#FFFFFF", "#BD6C3F",
              "#8175A6", "#8175A6", "#8175A6",
              "#469287", "#469287", 
              "#313564", "#FFC107", "#B02026")
names(myColors) <- levels(as.factor(c("ALR/Alpha", "BSR/Beta", "ct", "Dead", "GSAT", "GSATII", "GSATX", "HSat1A", "HSat1B", "HSat2", "HSat3", "Live")))

#filter to region of interest
#df2 <- df
#df2 <- filter(df2, start > 2000000)
#df2 <- filter(df2, stop < 2001000)

height=5
ggplot(data = df[order(df$region), ]) + 
    geom_segment(aes(x=start, y=chr, xend=stop+1000, yend=chr, color=region), alpha=1, size=height*2 ) +
    scale_color_manual(values=myColors) +
    #scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
    theme_classic() + theme(axis.line.y=element_line(colour="white")) +
    theme(legend.position="bottom") + theme(legend.key.width = unit(1,"cm")) + theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size=6))) + 
    xlab("Position (kbp)")

ggsave(paste0("HG002_ChrY_structure.png"), width=10, height=3)
#ggsave(paste0("browser_nhp_cens.png"), width=20, height=15)
#ggsave(paste0("browser_AS-HORs-vs-hg00733_censv0.3.png"), width=20, height=7)