library(data.table)
library(ggplot2)
library(ggridges)
library(scales)

setwd("../VGP/")
#setwd("/data/rhiea/genome10k/analysis/readlen")

### Number format
fancy_scientific <- function(d) {
  # turn in to character string in scientific notation
  d <- format(d, scientific = TRUE)
  # quote the part before the exponent to keep all the digits and turn the 'e+' into 10^ format
  d <- gsub("^(.*)e\\+", "'\\1'%*%10^", d)
  # convert 0x10^00 to 0
  d <- gsub("\\'0[\\.0]*\\'(.*)", "'0'", d)
  # convert 1x10^00 to 10^00
  d <- gsub("\\'1\\'\\%\\*\\%10\\^", "10\\^", d)
  # make 10^00 1
  d <- gsub("10\\^00", "\\'1\\'", d)
  # return this as an expression
  parse(text=d)
}

# Color: Set1
red = "#E41A1C"
green = "#4DAF4A"
blue = "#377EB8" # light blue = "#56B4E9"
purple = "#984EA3"  # purple = "#CC79A7"

my_col = c(red, green, blue, purple)

### Get genomes order
#genomes=fread("input/ED_Fig2/genome.list.srt", header = F)
genomes=fread("genome.list.srt", header = F)

colnames(genomes) = c("Genome")
head(genomes)

# Prepare 10X data
#dat.10x=fread("input/ED_Fig2/linked.len", header=T)
dat.10x=fread("linked.len", header=T)
head(dat.10x)
result <- data.frame()
for (g in genomes$Genome) {
  print(g)
  dat_genome1=dat.10x[dat.10x$genome==g,]
  dat_genome1 <- data.frame(Length=rexp(100000, rate=2/dat_genome1$molecule_length_mean))
  totalbp=sum(as.numeric(dat_genome1$Length))
  print(totalbp)
  head(dat_genome1)
  dat_genome1$Weight=dat_genome1$Length/totalbp
  dat_genome1$Length=dat_genome1$Length/1000
  dat_genome1$Platform="Linked."
  dat_genome1$Genome=c(g)
  result <- rbind(result, dat_genome1)
}
head(result)

# Fix the order
result$Genome = factor(result$Genome, levels = rev(unique(result$Genome)))

## Test on all genomes 10XG data
ggplot(result, aes(x = Length, y = Genome, fill = Platform, color = Platform, weight = Weight)) +
  geom_density_ridges(aes(height=..density..),
                      scale = 0.90, alpha = 0.6,
                      fill = green, color = green,
                      rel_min_height = 0.01, stat="density") +
  scale_y_discrete(expand = c(0, 0.5)) +   # will generally have to set the `expand` option
  coord_cartesian(clip = "off") +
  scale_x_continuous(trans="log10", limits=c(1,5000000), expand = c(0, 0), labels=fancy_scientific) +
  annotation_logticks(sides = "b") +
  theme_ridges(center_axis_labels = TRUE) +
  theme(axis.title.y = element_blank(), legend.position='top') +
  xlab("Molecule Length (kbp)")

# Load a test set of fGouWil2
dat_all = fread("all.1kb.weighted", header=F, nThread=4)
#dat_all = fread("input/ED_Fig2/fGouWil2.1kb.weighted", header=F, nThread=4)
#dat_all = cbind(dat_all, "fGouWil2")

# Bind the 10X results and the rest
colnames(dat_all) <- c("Length", "Weight", "Platform", "Genome")
dat_all=rbind(result, dat_all)
head(dat_all)
tail(dat_all)

dat_all$Platform = factor(dat_all$Platform, levels = c("CLR", "Linked.", "Opt", "Hi-C"))

# Plot density
pdf("ED_Fig2_readlen_wt.h10.pdf", width=5, height=10, pointsize=12, useDingbats=F)
ggplot(dat_all, aes(x = Length, y = Genome, fill = Platform, color = Platform, weight = Weight)) +
  scale_fill_manual(values = my_col) +
  scale_color_manual(values = my_col) +
  geom_density_ridges(aes(height=..density..),
                      scale = 0.95, alpha = 0.6,
                      rel_min_height = 0.01, stat="density") +
  scale_y_discrete(expand = c(0, 0.5)) +   # will generally have to set the `expand` option
  coord_cartesian(clip = "off") +
  scale_x_continuous(trans="log10", limits=c(1,5000000), expand = c(0, 0), labels=fancy_scientific) +
  annotation_logticks(sides = "b") +
  theme_ridges(center_axis_labels = TRUE) +
  theme(axis.title.y = element_blank(), legend.position='top', legend.title = element_blank()) +
  xlab("Molecule Length (kb)")
dev.off()
