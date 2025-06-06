# Set preference -> packages -> CDN

install.packages("rgeos", repos="http://R-Forge.R-project.org", type="source")
install.packages("ggplot2", dependencies = TRUE, type="source")
install.packages("argparser")
install.packages("ggridges", dependencies = TRUE, type="source")
install.packages("ggpubr", dependencies = TRUE)
library(cowplot)

# KaryoploteR
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = TRUE, type="source")

BiocManager::install("karyoploteR")
BiocManager::install("BiocFileCache")

.libPaths(c(Sys.getenv('R_LIBS'), .libPaths()))
# This does not work. Permission issue?
.libPaths()

Sys.getenv('R_LIBS')

# for loading large numbers
install.packages('bit64')
library(bit64)
