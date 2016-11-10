# source("./analysis/InitialSetup.R")

#### PCoA 

# All HJA Catchment
hja.db <- vegdist(OTUsREL, method = "bray")
hja.d.sorensen <- vegdist(OTUsREL, method = "bray", binary = T)
hja.d.jaccard <- vegdist(OTUsREL, method = "jaccard")
hja.pcoa <- cmdscale(hja.db, eig=TRUE)
hja.pcoa.sorensen <- cmdscale(hja.d.sorensen, eig = TRUE)
hja.pcoa.jaccard <- cmdscale(hja.d.jaccard, eig = TRUE)
var1 <- round(hja.pcoa$eig[1] / sum(hja.pcoa$eig),3) * 100
var2 <- round(hja.pcoa$eig[2] / sum(hja.pcoa$eig),3) * 100
var3 <- round(hja.pcoa$eig[3] / sum(hja.pcoa$eig),3) * 100
var1.sor <- round(hja.pcoa.sorensen$eig[1] / sum(hja.pcoa.sorensen$eig),3) * 100
var2.sor <- round(hja.pcoa.sorensen$eig[2] / sum(hja.pcoa.sorensen$eig),3) * 100
var3.sor <- round(hja.pcoa.sorensen$eig[3] / sum(hja.pcoa.sorensen$eig),3) * 100
var1.jac <- round(hja.pcoa.jaccard$eig[1] / sum(hja.pcoa.jaccard$eig),3) * 100
var2.jac <- round(hja.pcoa.jaccard$eig[2] / sum(hja.pcoa.jaccard$eig),3) * 100
var3.jac <- round(hja.pcoa.jaccard$eig[3] / sum(hja.pcoa.jaccard$eig),3) * 100


# Sediments Only
sediment <- OTUsREL[which(design$habitat == "sediment"),]
sediment.db <- vegdist(sediment, method = "bray")
sediment.dsorensen <- vegdist(decostand(sediment, method = "pa"), method = "bray")
sediment.d.jac <- vegdist(sediment, method = "jaccard")
sediment.pcoa <- cmdscale(sediment.db, eig=TRUE)
sed.design <- design[which(design$habitat == "sediment"),]
s.var1 <- round(sediment.pcoa$eig[1] / sum(sediment.pcoa$eig),3) * 100
s.var2 <- round(sediment.pcoa$eig[2] / sum(sediment.pcoa$eig),3) * 100
s.var3 <- round(sediment.pcoa$eig[3] / sum(sediment.pcoa$eig),3) * 100

# Water Only
water <- OTUsREL[which(design$habitat == "water"),]
water.db <- vegdist(water, method = "bray")
water.dsorensen <- vegdist(decostand(water, method = "pa"), method = "bray")
water.d.jac <- vegdist(water, method = "jaccard")
water.pcoa <- cmdscale(water.db, eig=TRUE)
water.design <- design[which(design$habitat == "water"),]
w.var1 <- round(water.pcoa$eig[1] / sum(water.pcoa$eig),3) * 100
w.var2 <- round(water.pcoa$eig[2] / sum(water.pcoa$eig),3) * 100
w.var3 <- round(water.pcoa$eig[3] / sum(water.pcoa$eig),3) * 100

# Lookout Creek Watershed Only
lc <- OTUsREL[which(design$watershed == "LC"),]
lc.db <- vegdist(lc, method = "bray")
lc.pcoa <- cmdscale(lc.db, eig=TRUE)
lc.design <- design[which(design$watershed == "LC"),]
lc.var1 <- round(lc.pcoa$eig[1] / sum(lc.pcoa$eig),3) * 100
lc.var2 <- round(lc.pcoa$eig[2] / sum(lc.pcoa$eig),3) * 100
lc.var3 <- round(lc.pcoa$eig[3] / sum(lc.pcoa$eig),3) * 100

# Watershed 01 Only
w1 <- OTUsREL[which(design$watershed == "WS01"),]
w1.db <- vegdist(w1, method = "bray")
w1.pcoa <- cmdscale(w1.db, eig=TRUE)
w1.design <- design[which(design$watershed == "WS01"),]
w1.var1 <- round(w1.pcoa$eig[1] / sum(w1.pcoa$eig),3) * 100
w1.var2 <- round(w1.pcoa$eig[2] / sum(w1.pcoa$eig),3) * 100
w1.var3 <- round(w1.pcoa$eig[3] / sum(w1.pcoa$eig),3) * 100

# Is habitat or order an important factor in community structure?
hja.permanova <- adonis(hja.db ~ design$habitat * design$watershed + design$order, permutations = 999)
capture.output(hja.permanova$aov.tab, file = "./tables/hja_permanova.txt")


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

