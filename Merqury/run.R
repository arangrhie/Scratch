###################################
### Test merqury plotting scripts
###################################


# Load plot_blob.R
setwd("./Merqury/")
source("./src/plot_blob.R")

### Fig. 4 c-e, smaller fig to enlarge axis fonts
blob_plot(dat = "input/triocanu_clr.hapmers.count", out = "output/triocanu.blob", w = 6.5, h = 4)
blob_plot(dat = "input/falcon_unzip.hapmers.count", out = "output/falcon_unzip.blob", w = 6.5, h = 4)
blob_plot(dat = "input/canu.hapmers.count", out = "output/canu.blob", w = 6, h = 4)

# Experimental scaled, log
blob_plot(dat = "input/triocanu_clr.hapmers.count", type = "scale", plot_log = TRUE, out = "output/triocanu.log")
blob_plot(dat = "input/falcon_unzip.hapmers.count", type = "scale", plot_log = TRUE, out = "output/falcon_unzip.log")
blob_plot(dat = "input/canu.hapmers.count", type = "scale", plot_log = TRUE, out = "output/canu.log")

dat = "input/canu.hapmers.count"
dat=read.table(dat, header=TRUE)
plot_blob_exact(dat)
plot_blob_scaled(dat)


### Fig. 4 f-k, fixed Y max
# Phase block and continuity, fix the YMAX
gsize = 130000000
ymax = 14000000

block = "input/triocanu.Col.100_20000.phased_block.sizes"
contig = "input/triocanu.Col.contig.sizes"
out = "output/triocanu.Col.100_20000.phased_block"

block_n(block = block, contig = contig, out = out, gsize = gsize, ymax = ymax, w = 3.5, h = 3)

block = "input/triocanu.Cvi.100_20000.phased_block.sizes"
contig = "input/triocanu.Cvi.contig.sizes"
out = "output/triocanu.Cvi.100_20000.phased_block"

block_n(block = block, contig = contig, out = out, gsize = gsize, ymax = ymax, w = 3.5, h = 3)

block = "input/falcon_unzip.pri.100_20000.phased_block.sizes"
contig = "input/falcon_unzip.pri.contig.sizes"
out = "output/falcon_unzip.pri.100_20000.phased_block"

block_n(block = block, contig = contig, out = out, gsize = gsize, ymax = ymax, w = 3.5, h = 3)

block = "input/falcon_unzip.alt.100_20000.phased_block.sizes"
contig = "input/falcon_unzip.alt.contig.sizes"
out = "output/falcon_unzip.alt.100_20000.phased_block"

block_n(block = block, contig = contig, out = out, gsize = gsize, ymax = ymax, w = 3.5, h = 3)

block = "input/canu.100_20000.phased_block.sizes"
contig = "input/canu.contig.sizes"
out = "output/canu.100_20000.phased_block"

block_n(block = block, contig = contig, out = out, gsize = gsize, ymax = ymax, w = 3.5, h = 3)








