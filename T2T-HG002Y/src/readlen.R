library(data.table)
library(ggplot2)
library(ggridges)
library(scales)

setwd("T2T-HG002Y/")

### .len generated with
#  for len in $(ls *.len.gz)
#  do
#    zcat $len | awk -v name=${len/.len.gz} 'NR%1000 {print name"\t"$0}' >> $out
#  done
#######################################

plot_density = function(dat = NULL) {
    ggplot(data = dat, aes(x=Len, colour = SMRT)) +
        geom_density(alpha=.2) +
        theme_bw() +
        theme(axis.title = element_text(size=6, face = "bold"),
              axis.text  = element_text(size=5),
              legend.key.size = unit(4, "mm"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=5),
              legend.margin=margin(c(0,0,0,0)))
}

dat1 = fread("input/sequelII_20kb_20210812.len", header=F, nThread=4)
dat2 = fread("input/sequelII_20kb_20210812_WUSTL.len", header=F, nThread=4)
colnames(dat1) = c("SMRT", "Len")
colnames(dat2) = c("SMRT", "Len")
dat=rbind(dat1, dat2)
colnames(dat) = c("SMRT", "Len")

plot_density(dat1)
ggsave(paste("output/", "sequelII_20kb_20210812", ".len.png", sep=""),
       width = 4, height = 3)
plot_density(dat2)
ggsave(paste("output/", "sequelII_20kb_20210812_WUSTL", ".len.png", sep=""),
       width = 4, height = 3)
plot_density(dat)
