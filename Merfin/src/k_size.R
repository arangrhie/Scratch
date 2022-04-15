library(ggplot2)

setwd("Merfin")

# x = g
eq_0 = function(x) {
    e = 0.00000001
    log((1000000000 * x * (1-e)/e), base = 4)
}

eq_1 = function(x) {
    e = 0.001
    log((1000000000 * x * (1-e)/e), base = 4)
}

eq_2 = function(x) {
    e = 0.01
    log((1000000000 * x * (1-e)/e), base = 4)
}

from=0.1
to=10

ggplot(data.frame(x=c(from, to)), aes(x=x)) + 
    #stat_function(fun=eq_0, linetype = "dotdash", aes(colour = "0.00000001"), size = 0.3) +
    stat_function(fun=eq_1, linetype = "solid", aes(colour = "0.001"), size = 0.3) +
    stat_function(fun=eq_2, linetype = "longdash", aes(colour = "0.01"), size = 0.3) +
    #scale_colour_manual("Tolerable Collision Rate", values = c("darkgray", "black", "gray")) +
    scale_colour_manual("Tolerable Collision Rate", values = c("black", "gray")) +
    geom_vline(xintercept = 3.2, size = 0.2) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    xlab("Genome Size (Gb)") +
    ylab("k") +
    theme(legend.position='bottom') +
    theme(axis.title = element_text(size=6, face = "bold"),
          axis.text  = element_text(size=6),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          legend.margin=margin(c(0,0,0,0)))
ggsave("output/k_size.pdf", units = "mm", width = 60, height = 45)
ggsave("output/k_size.png", units = "mm", width = 60, height = 45)
