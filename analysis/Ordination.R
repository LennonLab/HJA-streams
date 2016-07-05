# source("./InitialSetup.R")

#### PCoA 

# All HJA Catchment
hja.db <- vegdist(OTUsREL, method = "bray", upper = TRUE, diag = TRUE)
hja.pcoa <- cmdscale(hja.db, eig=TRUE)
var1 <- round(hja.pcoa$eig[1] / sum(hja.pcoa$eig),3) * 100
var2 <- round(hja.pcoa$eig[2] / sum(hja.pcoa$eig),3) * 100
var3 <- round(hja.pcoa$eig[3] / sum(hja.pcoa$eig),3) * 100

# Sediments Only
sediment <- OTUsREL[which(design$habitat == "sediment"),]
sediment.db <- vegdist(sediment, method = "bray", diag = T)
sediment.pcoa <- cmdscale(sediment.db, eig=TRUE)
sed.design <- design[which(design$habitat == "sediment"),]
s.var1 <- round(sediment.pcoa$eig[1] / sum(sediment.pcoa$eig),3) * 100
s.var2 <- round(sediment.pcoa$eig[2] / sum(sediment.pcoa$eig),3) * 100
s.var3 <- round(sediment.pcoa$eig[3] / sum(sediment.pcoa$eig),3) * 100

# Water Only
water <- OTUsREL[which(design$habitat == "water"),]
water.db <- vegdist(water, method = "bray", diag = T)
water.pcoa <- cmdscale(water.db, eig=TRUE)
water.design <- design[which(design$habitat == "water"),]
w.var1 <- round(water.pcoa$eig[1] / sum(water.pcoa$eig),3) * 100
w.var2 <- round(water.pcoa$eig[2] / sum(water.pcoa$eig),3) * 100
w.var3 <- round(water.pcoa$eig[3] / sum(water.pcoa$eig),3) * 100

# Lookout Creek Watershed Only
lc <- OTUsREL[which(design$watershed == "LC"),]
lc.db <- vegdist(lc, method = "bray", diag = T)
lc.pcoa <- cmdscale(lc.db, eig=TRUE)
lc.design <- design[which(design$watershed == "LC"),]
lc.var1 <- round(lc.pcoa$eig[1] / sum(lc.pcoa$eig),3) * 100
lc.var2 <- round(lc.pcoa$eig[2] / sum(lc.pcoa$eig),3) * 100
lc.var3 <- round(lc.pcoa$eig[3] / sum(lc.pcoa$eig),3) * 100

# Watershed 01 Only
w1 <- OTUsREL[which(design$watershed == "WS01"),]
w1.db <- vegdist(w1, method = "bray", diag = T)
w1.pcoa <- cmdscale(w1.db, eig=TRUE)
w1.design <- design[which(design$watershed == "WS01"),]
w1.var1 <- round(w1.pcoa$eig[1] / sum(w1.pcoa$eig),3) * 100
w1.var2 <- round(w1.pcoa$eig[2] / sum(w1.pcoa$eig),3) * 100
w1.var3 <- round(w1.pcoa$eig[3] / sum(w1.pcoa$eig),3) * 100

# Is habitat or order an important factor in community structure?
adonis(hja.db ~ design$order, permutations = 999)
hja.permanova <- adonis(hja.db ~ design$habitat * design$watershed + design$order, permutations = 999)
capture.output(hja.permanova$aov.tab, file = "../tables/hja_permanova.txt")



#### Figure 4: Sediment vs. Water PCoA
png(filename = "../figures/HJA_PCoA.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 3, 2) + 0.1)
plot(hja.pcoa$points[ ,1], hja.pcoa$points[ ,2], ylim = c(-0.5, 0.4), xlim = c(-.35, .35),
     xlab = paste("PCoA 1 (", var1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", var2, "%)", sep = ""), 
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, at = c(-0.3,-.15,0,.15,.3), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, at = c(-0.3,-.15,0,.15,.3), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(hja.pcoa$points[which(design$habitat == "sediment"),1], 
       hja.pcoa$points[which(design$habitat == "sediment"),2],
       pch=21, cex=2, bg="grey")
points(hja.pcoa$points[which(design$habitat == "water"),1], 
       hja.pcoa$points[which(design$habitat == "water"),2],
       pch=24, cex=2, bg="white")
legend("topright", c("Water", "Sediment"),
       pt.bg = c("white", "grey"), pch = c(24,21), cex = 1.5, bty = "n")
dev.off()
graphics.off()

img <- readPNG("../figures/HJA_PCoA.png")
grid.raster(img)
