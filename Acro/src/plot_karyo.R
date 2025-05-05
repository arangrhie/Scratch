library(karyoploteR)
library(BiocFileCache)

getwd()
setwd("Acro")

chm13.genome <- toGRanges("input/chm13v2.0.sizes.bed")
chm13.cytobands <- toGRanges("input/chm13v2.0_cytobands_allchrs.txt")

# Issues remaining in GRCh38
issues_hg38_all = toGRanges("input/chm13v2.0_GRCh38issues_lifted.20230315.fmt.bed")

# Non-syntenic to GRCh38 (Unique in CHM13 compared to GRCh38)
unique_to_hg38 = toGRanges("input/chm13v2-unique_to_hg38.bed")
unique_to_hg19 = toGRanges("input/chm13v2-unique_to_hg19.bed")

# Short-read accessibility mask
access_short = toGRanges("input/hs1.combined_mask.mrg_10kb.bed")

bfc <- BiocFileCache(ask=FALSE)
censat.file <- bfcrpath(bfc, "input/chm13v2.0.cenSatv2.1.bed")
chm13.censat = toGRanges(censat.file)

# For using cytobands
# kp <- plotKaryotype(genome = chm13.genome, cytobands = chm13.cytobands,
# mask = NA, plot.type = 1, non.overlapping = FALSE)

# Using autotrack; doesn't allow to overlay tracks
# at <- autotrack(current.track = 1, total.tracks = 2, r0 = 0, r1 = 0.5)
# at <- autotrack(current.track = 2, total.tracks = 2, r0 = 0.5, r1 = 1)

pdf("output/chm13v2.0_accessibility.pdf", width = 3, height = 6.5)

plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$data1inmargin <- 0
plot.params$data1outmargin <- 0

track1_r0 <- 0.05
track1_r1 <- 0.25
track2_r0 <- 0.3
track2_r1 <- 0.5
track3_r0 <- 0.55
track3_r1 <- 0.7

kp <- plotKaryotype(genome = chm13.genome, plot.params = plot.params,
  plot.type = 1, non.overlapping = FALSE)
kpPlotRegions(kp, data=censat, col = censat$itemRgb, border = NULL, data.panel = "ideogram")
kpPlotRegions(kp, data=access_short, col = "#666666", r0 = track1_r0, r1 = track1_r1, border = NA)
kpPlotRegions(kp, data=unique_to_hg19, col = "#CCCCCC", r0 = track2_r0, r1 = track2_r1, border = NA)
kpPlotRegions(kp, data=unique_to_hg38, col = "#AAAAAA", r0 = track2_r0, r1 = track2_r1, border = NA)
kpPlotRegions(kp, data=issues_hg38_all, col = "#ff9200", r0 = track3_r0, r1 = track3_r1)

dev.off()

plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$data1inmargin <- 0
plot.params$data1outmargin <- 0
plot.params$data1height <- 0.5

pdf("output/chr13p.pdf", width = 4, height = 0.5)
# Acro pattern
kp <- plotKaryotype(genome = chm13.genome, cytobands = chm13.cytobands, plot.params = plot.params,
					chromosomes = "chr13",
          zoom = "chr13:1-17514208",
          plot.type = 1, non.overlapping = FALSE)
kpPlotRegions(kp, data=censat, col = censat$itemRgb, border = NULL, data.panel = "ideogram")
kpAddBaseNumbers(kp)
dev.off()

pdf("output/chr14p.pdf", width = 4, height = 0.5)
kp <- plotKaryotype(genome = chm13.genome, cytobands = chm13.cytobands, plot.params = plot.params,
					chromosomes = "chr14",
          zoom = "chr14:1-12740172",
          plot.type = 1, non.overlapping = FALSE)
kpPlotRegions(kp, data=censat, col = censat$itemRgb, border = NULL, data.panel = "ideogram")
kpAddBaseNumbers(kp)
dev.off()

pdf("output/chr15p.pdf", width = 4, height = 0.5)
kp <- plotKaryotype(genome = chm13.genome, cytobands = chm13.cytobands, plot.params = plot.params,
					chromosomes = "chr15",
          zoom = "chr15:1-17694128",
          plot.type = 1, non.overlapping = FALSE)
kpPlotRegions(kp, data=censat, col = censat$itemRgb, border = NULL, data.panel = "ideogram")
kpAddBaseNumbers(kp)
dev.off()

pdf("output/chr21p.pdf", width = 4, height = 0.5)
kp <- plotKaryotype(genome = chm13.genome, cytobands = chm13.cytobands, plot.params = plot.params,
					chromosomes = "chr21",
          zoom = "chr21:1-11310801",
          plot.type = 1, non.overlapping = FALSE)
kpPlotRegions(kp, data=censat, col = censat$itemRgb, border = NULL, data.panel = "ideogram")
kpAddBaseNumbers(kp)
dev.off()

pdf("output/chr22p.pdf", width = 4, height = 0.5)
kp <- plotKaryotype(genome = chm13.genome, cytobands = chm13.cytobands, plot.params = plot.params,
					chromosomes = "chr22",
          zoom = "chr22:1-16860439",
          plot.type = 1, non.overlapping = FALSE)
kpPlotRegions(kp, data=censat, col = censat$itemRgb, border = NULL, data.panel = "ideogram")
kpAddBaseNumbers(kp)
dev.off()
