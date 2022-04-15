library(ggplot2)
library(scales)

dat=read.table("input/0904_fixes/lengths.hist.indel", header=FALSE)
head(dat)
names(dat) <- c("Size", "Frequency")
summary(dat)

ggplot(dat, aes(x = Size, y=Frequency)) +
    geom_col() +
    theme_bw() +
    scale_x_continuous() +
    scale_y_continuous()

ggsave("output/0904_fixes/indel_size.png", width = 3, height = 3)


ggplot(dat, aes(x = Size, y=Frequency)) +
    geom_col() +
    theme_bw() +
    scale_x_continuous() +
    scale_y_log10("Num. INDELs",
                  breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))

ggsave("output/0904_fixes/indel_size_log.png", width = 3, height = 3)
