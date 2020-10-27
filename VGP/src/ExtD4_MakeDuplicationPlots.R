library(ggplot2)
library(RColorBrewer)
library(stringr)

setwd("VGP")
#
# Test data that is ignored on the command line
#
reptab <- read.table("input/ExtD4/all_repeat_summary.with_wm.orig.tsv")
head(reptab)
asmtab <- read.table("input/ExtD4/assembly_size.txt")
head(asmtab)
allRepeats <- read.table("input/ExtD4/all_species.orig.tab")
head(allRepeats)
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

pdf(collapsedRepeatContentHighCopy)
ggplot(data=reptab[highCopy,], aes(x=rep, y=V3, fill=V1)) +
  geom_bar(stat="identity") +
  theme(axis.title=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
        panel.grid.major = element_line(colour="lightgrey", size=0.5, linetype="solid"),
        panel.background=element_blank(), panel.grid.minor = element_blank() ) +
  xlab("Repeat type") +
  ylab("Total bases") +
  labs(fill="assembly") +
  scale_fill_manual(name="Species", values=p)
dev.off()

pdf(collapsedRepeatContentLowCopy)
ggplot(data=reptab[lowCopy,], aes(x=rep, y=V3, fill=V1)) +
  geom_bar(stat="identity") +
  theme(axis.title=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
        panel.grid.major = element_line(colour="lightgrey", size=0.5, linetype="solid"),
        panel.background=element_blank(), panel.grid.minor = element_blank() ) +
  xlab("Repeat type") +
  ylab("Total bases") +
  labs(fill="assembly") +
  scale_fill_manual(name="Species", values=p)
dev.off()       

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

plt <- ggplot(allRepeats,aes(x=cumulativeLengthNorm, y=orderedIdentity, group=V1, colour=V1, size=orderedLength))
pdf(allSpeciesPlot)
plt + geom_point() +
  geom_hline(yintercept=0.9) +
  ylab("Fraction repeat") + 
  xlab("Cumulative collapsed duplication / genome length (Gb) ") +
  scale_colour_discrete(name="Assembly") + scale_size(name="Collapse length") +
  geom_line(lwd=0.5) +
  scale_colour_manual(name="Species", values=p)
dev.off()
                                       
