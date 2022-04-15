library(ggplot2)
library(scales)

dat=read.table("input/0904/Het_Size.len", header=T)
head(dat)

ggplot(dat, aes(x = Size)) +
    geom_histogram(binwidth = 100) +
    theme_bw() +
    scale_x_continuous() +
    scale_y_continuous()
ggsave("output/0904.het_size.png", width=3, height = 3)
