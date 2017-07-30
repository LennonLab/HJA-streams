source("./analysis/InitialSetup.R")
source("./analysis/HJA-Functions.R")

#### PCoA 

# All HJA Catchment
hja.pcoa <- run.pcoa(comm = OTUsREL.hel, method = "euclidean")

# Sediments Only
sed.pcoa <- run.pcoa(comm = OTUsREL.hel[which(design$habitat == "sediment"),])

# Water Only
water.pcoa <- run.pcoa(comm = OTUsREL.hel[which(design$habitat == "water"),])


# Lookout Creek Watershed Only
lookout.pcoa <- run.pcoa(comm = OTUsREL.hel[which(design$watershed == "LC"),])

# Watershed 01 Only
ws01.pcoa <- run.pcoa(comm = OTUsREL.hel[which(design$watershed == "WS01"),])

# Is habitat or order an important factor in community structure?
hja.permanova <- adonis(hja.pcoa$dist.matrix ~ design$habitat * design$watershed + design$order, permutations = 999)
capture.output(hja.permanova$aov.tab, file = "./tables/hja_permanova.txt")



ord <- hja.pcoa$pcoa
scrs <- scores(ord)
xlim <- extendrange(x = .1, r = range(scrs[,1]))
ylim <- extendrange(x = .1, r = range(scrs[,2]))

with(design, levels(habitat))
cols <- c("white", "grey")

plot.new()
plot.window(xlim = xlim, ylim = ylim, asp = 1)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
with(design, points(scrs, col = "black", pch = 21, cex = 2, bg = cols[habitat]))
with(design, legend("topright", legend = levels(habitat), bty = "n",
                      col = cols, pch = 21, pt.bg = cols))
axis(side = 1, lwd.ticks = 2, cex.axis = 1.2)
axis(side = 2, lwd.ticks = 2, cex.axis = 1.2)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2)
title(xlab = "PCoA 1", ylab = "PCoA 2")
box(lwd = 2)




#### Figure 4: Sediment vs. Water PCoA
png(filename = "./figures/HJA_PCoA_Bray.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 3, 2) + 0.1)
plot(hja.pcoa$points[ ,1], hja.pcoa$points[ ,2], ylim = c(-0.5, 0.5), xlim = c(-.4, .4),
     xlab = paste("PCoA 1 (", var1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", var2, "%)", sep = ""), 
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, at = c(-0.4,-.2,0,.2,.4), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, at = c(-0.4,-.2,0,.2,.4), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(hja.pcoa$points[which(design$habitat == "sediment"),1], 
       hja.pcoa$points[which(design$habitat == "sediment"),2],
       pch=21, cex=2, bg="grey")
points(hja.pcoa$points[which(design$habitat == "water"),1], 
       hja.pcoa$points[which(design$habitat == "water"),2],
       pch=24, cex=2, bg="white")
legend("topleft", c("Water", "Sediment"),
       pt.bg = c("white", "grey"), pch = c(24,21), cex = 1.5, bty = "n")
ordiellipse(hja.pcoa, design$habitat, conf = 0.95)
dev.off()
graphics.off()
img <- readPNG("./figures/HJA_PCoA_Bray.png")
grid.raster(img)

png(filename = "./figures/HJA_PCoA_Sorensen.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 3, 2) + 0.1)
plot(hja.pcoa.sorensen$points[ ,1], hja.pcoa.sorensen$points[ ,2], ylim = c(-0.5, 0.5), xlim = c(-.5, .5),
     xlab = paste("PCoA 1 (", var1.sor, "%)", sep = ""),
     ylab = paste("PCoA 2 (", var2.sor, "%)", sep = ""), 
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, at = c(-0.5,-.25,0,.25,.5), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, at = c(-0.5,-.25,0,.25,.5), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(hja.pcoa.sorensen$points[which(design$habitat == "sediment"),1], 
       hja.pcoa.sorensen$points[which(design$habitat == "sediment"),2],
       pch=21, cex=2, bg="grey")
points(hja.pcoa.sorensen$points[which(design$habitat == "water"),1], 
       hja.pcoa.sorensen$points[which(design$habitat == "water"),2],
       pch=24, cex=2, bg="white")
legend("topleft", c("Water", "Sediment"),
       pt.bg = c("white", "grey"), pch = c(24,21), cex = 1.5, bty = "n")
ordiellipse(hja.pcoa.sorensen, design$habitat, conf = .95)
dev.off()
graphics.off()
img <- readPNG("./figures/HJA_PCoA_Sorensen.png")
grid.raster(img)


png(filename = "./figures/HJA_PCoA_jac.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 3, 2) + 0.1)
plot(hja.pcoa.jaccard$points[ ,1], hja.pcoa.jaccard$points[ ,2], 
     ylim = c(-0.5, 0.5), xlim = c(-.5, .5),
     xlab = paste("PCoA 1 (", var1.jac, "%)", sep = ""),
     ylab = paste("PCoA 2 (", var2.jac, "%)", sep = ""), 
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(hja.pcoa.jaccard$points[which(design$habitat == "sediment"),1], 
       hja.pcoa.jaccard$points[which(design$habitat == "sediment"),2],
       pch=21, cex=2, bg="grey")
points(hja.pcoa.jaccard$points[which(design$habitat == "water"),1], 
       hja.pcoa.jaccard$points[which(design$habitat == "water"),2],
       pch=24, cex=2, bg="white")
legend("topleft", c("Water", "Sediment"),
       pt.bg = c("white", "grey"), pch = c(24,21), cex = 1.5, bty = "n")
ordiellipse(hja.pcoa.jaccard, design$habitat, conf = .95)
dev.off()
graphics.off()
grid.raster(png::readPNG("./figures/HJA_PCoA_jac.png"))

