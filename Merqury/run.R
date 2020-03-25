###################################
### Test merqury plotting scripts
###################################


# Load plot_blob.R
setwd("./Merqury/")
source("./src/plot_blob.R")

blob_plot(dat = "input/canu.hapmers.count", out = "output/canu")
