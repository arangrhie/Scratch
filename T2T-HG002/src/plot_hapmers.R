library(ggplot2)

setwd("T2T-HG002/")

dat=read.table("input/hifi_read.pos_hapmer.srt.txt", header=FALSE)
dat=read.table("input/ont_read.pos_hapmer.srt.txt", header=FALSE)
head(dat)

names(dat) = c("ReadPos", "Frequency_in_Illumina", "Hapmer")
ggplot(data = dat, aes(ReadPos, Frequency_in_Illumina)) +
    geom_point(aes(colour = Hapmer), size = 0.5) +
    geom_line(y=31) +
    geom_line(y=64) +
    theme_classic()
