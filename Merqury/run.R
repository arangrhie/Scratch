###################################
### Test merqury plotting scripts
###################################


# Load plot_blob.R
setwd("./Merqury/")
source("./src/plot_blob.R")

blob_plot(dat = "input/triocanu_clr.hapmers.count", out = "output/triocanu", type = "scale")
blob_plot(dat = "input/canu.hapmers.count", out = "output/canu", type = "scale")
blob_plot(dat = "input/falcon_unzip.hapmers.count", out = "output/falcon_unzip", type = "scale")

blob_plot(dat = "input/triocanu_clr.hapmers.count", type = "scale", plot_log = TRUE, out = "output/triocanu.log")
blob_plot(dat = "input/canu.hapmers.count", type = "scale", plot_log = TRUE, out = "output/canu.log")
blob_plot(dat = "input/falcon_unzip.hapmers.count", type = "scale", plot_log = TRUE, out = "output/falcon_unzip.log")

dat = "input/canu.hapmers.count"
dat=read.table(dat, header=TRUE)
plot_blob_exact(dat)
plot_blob_scaled(dat)
