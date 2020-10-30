library(ggplot2)
library(RColorBrewer)
library(stringr)

fancy_scientific <- function(d) {
  # turn in to character string in scientific notation
  d <- format(d, scientific = TRUE)
  # quote the part before the exponent to keep all the digits and turn the 'e+' into 10^ format
  d <- gsub("^(.*)e\\+", "'\\1'%*%10^", d)
  # convert 0x10^00 to 0
  d <- gsub("\\'0[\\.0]*\\'(.*)", "'0'", d)
  # return this as an expression
  parse(text=d)
}

shared_theme <- function() {
  list(
    theme_bw(),
    theme(plot.title = element_blank(),
          #axis.line = element_line(size = 0.1),
          #legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.margin=margin(c(0,0,0,0)),
          axis.text = element_text(size = 5),
          plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
          axis.ticks = element_line(size=0.1),
          axis.title = element_text(size=6))
  )
}

shared_theme2 <- function() {
  list(
    theme(legend.position="none"),
    scale_y_continuous(labels = fancy_scientific)
  )
}

setwd("VGP")

#
# Test data that is ignored on the command line
#
allRepeats <- read.table("input/ExtD4/all_species.orig.tab")
head(allRepeats)
reptab <- read.table("input/ExtD4/all_repeat_summary.with_wm.orig.tsv")
head(reptab)
asmtab <- read.table("input/ExtD4/assembly_size.txt")
head(asmtab)
collapsedRepeatContentHighCopy <- "output/ExtD4_collapsed.dist.hc.pdf"
collapsedRepeatContentLowCopy <- "output/ExtD4_collapsed.dist.lc.pdf"
allSpeciesPlot <- "output/ExtD4_all_dups.orig.pdf"
cl <- "orig"

#
# Running as a script
#
args <- commandArgs(trailingOnly=T)

reptab <- read.table(args[1])
asmtab <- read.table(args[2])
collapsedRepeatContentHighCopy <- args[3]
collapsedRepeatContentLowCopy <- args[4]
allSpeciesTab <- args[5]
allSpeciesPlot <- args[6]
cl <- args[7]
allRepeats <- read.table(allSpeciesTab)

asmNames <- asmtab$V1
p        <- rainbow(length(asmNames))
names(p) <- asmNames
p


f <- as.numeric(as.factor(allRepeats$V1))
u <- unique(allRepeats$V1)

allRepeats$order <- 0
allRepeats$length <- allRepeats$V4-allRepeats$V3

idx <- sapply(u, function(i) which(allRepeats$V1 == i))

for(i in seq(1,length(idx))) {
  allRepeats$order[idx[[i]]] <- order(allRepeats$V9[idx[[i]]]);
}

allRepeats$orderedLength <- 0
allRepeats$pos <- 0
allRepeats$orderedIdentity <- 0
allRepeats$cumulativeLength <- 0
allRepeats$cumulativeLengthNorm <- 0
asmtab$normLength <- asmtab$V2/1E9
for(i in seq(1,length(idx))) {
  allRepeats$orderedIdentity[idx[[i]]] <- sort(allRepeats$V9[idx[[i]]]);
  allRepeats$pos[idx[[i]]] <- seq(1,length(idx[[i]]))
  o <- order(allRepeats$V9[idx[[i]]])
  allRepeats$orderedLength[idx[[i]]] <- allRepeats$V4[idx[[i]][o]]-allRepeats$V3[idx[[i]][o]]
  allRepeats$cumulativeLength[idx[[i]]] <- cumsum(allRepeats$orderedLength[idx[[i]]])
  allRepeats$cumulativeLengthNorm[idx[[i]]] <- cumsum(allRepeats$orderedLength[idx[[i]]])/asmtab$normLength[which(asmtab$V1 == allRepeats$V1[idx[[i]][1]])[1]]
}

pdf(allSpeciesPlot, width = 3.3, height = 2.9)
plt <- ggplot(allRepeats,aes(x=cumulativeLengthNorm, y=orderedIdentity, group=V1, colour=V1, size=orderedLength)) +
  geom_point(alpha = 0.5) +
  geom_line(lwd = 0.3) +
  geom_hline(yintercept = 0.9, lwd = 0.2) +
  ylab("% Repeat Masked") + 
  xlab("Cumulative Collapsed Bases / Genome Length (Gb) ") +
  scale_colour_discrete(name="Assembly") +
  scale_size(name="Collapse Length", range = c(1, 4), labels = fancy_scientific) +
  scale_color_manual(name = "Species", values = p) +
  scale_x_continuous(labels = fancy_scientific) +
  shared_theme() +
  guides(linetype=guide_legend(keyheight = unit(2, "mm")),
         shape=guide_legend(keyheight = unit(2, "mm")),
         color = guide_legend(keyheight = unit(2, "mm")),
         size = guide_legend(keyheight = unit(3, "mm")))
plt
dev.off()


#
# Pack a couple of types together
reptab$V2 <- as.character(reptab$V2)
idx <- which(reptab$V2 == "WindowMasker")
reptab$V2[idx] <- "Unknown"
idx <- which(reptab$V2 == "Low_complexity")
reptab$V2[idx] <- "Simple_repeat"
head(reptab)

#
# Remove underscores, etc. to make names readable.
#
underscore   <- str_replace(reptab$V2, "_", " ")
noSlash      <- str_replace(underscore, "\\?", " ")
reptab$names <- noSlash
ag           <- aggregate(x=reptab$V3, by=list(reptab$names),FUN=sum)
ago          <- order(ag$x, decreasing=T)

#
# convert to factor for ggplot help
#

reptab$rep <- factor(reptab$names, levels=ag$Group.1[ago])

#
# "highCopy" - duplications with many bases, not necessarily high copy.
highCopy <- which(reptab$V2 == "Unknown" | reptab$V2 == "Simple_repeat" | reptab$V2=="Satellite" | reptab$V2=="Low_complexity")
lowCopy  <- setdiff(seq(1,length(reptab$V2)), highCopy)
highCP <- reptab[highCopy,]
lowCP <- reptab[lowCopy,]
tail(reptab[highCopy,])
head(reptab[lowCopy,])

#
# gather lowCopy to one bar in highCopy
others=data.frame()
for (g in unique(lowCP$V1) ) {
  gDat <- lowCP[lowCP$V1 == g,]
  sum <- sum(gDat$V3)
  if (sum > 0) {
    gDatNew <- data.frame("V1" = c(g), "V2" = c("Others"), "V3" = c(sum), names = c("Others"), rep = c("Others"))
    others <- rbind(others, gDatNew)
  }
}
highCP <- rbind(highCP, others)
tail(highCP)

pdf(collapsedRepeatContentHighCopy, width = 2.4, height = 1)
ggplot(data=highCP, aes(x=rep, y=V3, fill=V1)) +
  geom_bar(stat="identity") +
  shared_theme() +
  shared_theme2() +
  xlab("Repeat Type") +
  ylab("Total Bases") +
  scale_fill_manual(name="Species", values=p)
dev.off()

pdf(collapsedRepeatContentLowCopy, width = 2.4, height = 2)
ggplot(data=lowCP, aes(x=rep, y=V3, fill=V1)) +
  geom_bar(stat="identity") +
  shared_theme() +
  shared_theme2() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) +
  xlab("Repeat Type") +
  ylab("Total Bases") +
  scale_fill_manual(name="Species", values=p)
dev.off()       

                    
theme(axis.title=element_text(size=12),
      axis.text.y=element_text(size=12),
      axis.text=element_text(size=12),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
      panel.grid.major = element_line(colour="lightgrey", size=0.5, linetype="solid"),
      panel.background=element_blank(), panel.grid.minor = element_blank() )