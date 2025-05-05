chr = "chr22"

xl = c(0,17500000)
yl = c(80,100)
zl = c(80,100)
smooth = 100000
pal = hcl.colors(100, palette="YlOrRd", rev = TRUE)
linecol = "darkgray"

a <- read.table(paste(chr,"chr13.idy",sep="."))
b <- read.table(paste(chr,"chr14.idy",sep="."))
c <- read.table(paste(chr,"chr15.idy",sep="."))
d <- read.table(paste(chr,"chr21.idy",sep="."))
e <- read.table(paste(chr,"chr22.idy",sep="."))

a$V3 <- ksmooth(a$V1, a$V2, kernel="box", bandwidth=smooth)$y
b$V3 <- ksmooth(b$V1, b$V2, kernel="box", bandwidth=smooth)$y
c$V3 <- ksmooth(c$V1, c$V2, kernel="box", bandwidth=smooth)$y
d$V3 <- ksmooth(d$V1, d$V2, kernel="box", bandwidth=smooth)$y
e$V3 <- ksmooth(e$V1, e$V2, kernel="box", bandwidth=smooth)$y

pdf(file=paste(chr,"pdf",sep="."), width=20, height=1.5)

par(mfrow=c(4,1))
par(mar = c(0, 2, 0, 0))
par(oma = c(0, 0, 0, 0))
par(xaxt="n")
par(yaxt="s")
par(las=0)
par(bty="n")
par(tck=NA)

if ( chr != "chr13" ) {
  plot(1, type="n", xlim=xl, ylim=yl)
  image(a$V1, 1, as.matrix(a$V3), zlim=zl, col = pal, add=TRUE, useRaster=TRUE)
  abline(h=100, lty=3, col="gray")
  lines(a$V1, a$V3, lwd=2, pch=".", col=linecol)
}

if ( chr != "chr14" ) {
  plot(1, type="n", xlim=xl, ylim=yl)
  image(b$V1, 1, as.matrix(b$V3), zlim=zl, col = pal, add=TRUE, useRaster=TRUE)
  abline(h=100, lty=3, col="gray")
  lines(b$V1, b$V3, lwd=2, pch=".", col=linecol)
}

if ( chr != "chr15" ) {
  plot(1, type="n", xlim=xl, ylim=yl)
  image(c$V1, 1, as.matrix(c$V3), zlim=zl, col = pal, add=TRUE, useRaster=TRUE)
  abline(h=100, lty=3, col="gray")
  lines(c$V1, c$V3, lwd=2, pch=".", col=linecol)
}

if ( chr != "chr21" ) {
  plot(1, type="n", xlim=xl, ylim=yl)
  image(d$V1, 1, as.matrix(d$V3), zlim=zl, col = pal, add=TRUE, useRaster=TRUE)
  abline(h=100, lty=3, col="gray")
  lines(d$V1, d$V3, lwd=2, pch=".", col=linecol)
}

if ( chr != "chr22" ) {
  plot(1, type="n", xlim=xl, ylim=yl)
  image(e$V1, 1, as.matrix(e$V3), zlim=zl, col = pal, add=TRUE, useRaster=TRUE)
  abline(h=100, lty=3, col="gray")
  lines(e$V1, e$V3, lwd=2, pch=".", col=linecol)
}

dev.off()