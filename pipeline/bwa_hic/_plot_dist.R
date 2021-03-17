library(ggplot2)
library("ggsci")

dat1=read.table("dovetail8_to_pacbio_falcon_pa_pilon_ponly.merged.dist", header=F, sep="\t")
dat2=read.table("dovetail9_to_pacbio_falcon_pa_pilon_ponly.merged.dist", header=F, sep="\t")
dat3=read.table("arima_to_pacbio_falcon_pa_pilon_ponly.merged.dist", header=F, sep="\t")
dat4=read.table("phase_to_pacbio_falcon_pa_pilon_ponly.merged.dist", header=F, sep="\t")

dat1=cbind(dat1, type=c("dovetail8"))
dat2=cbind(dat2, type=c("dovetail9"))
dat3=cbind(dat3, type=c("arima"))
dat4=cbind(dat4, type=c("phase"))


dat=rbind(dat1, dat2, dat3, dat4)

#head(dat)
remove(dat1, dat2, dat3, dat4)

colnames(dat) <- c("hiC", "Type")

png(filename="hiC_dist_density.png")
p <- ggplot(dat, aes(dat$hiC, colour=dat$Type)) + labs(x="Distance (Insert Size)") + scale_x_log10() + theme_bw() + geom_density() + scale_color_npg() + scale_fill_npg()
plot(p)
dev.off()


