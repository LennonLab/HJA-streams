# Catchment scale elevational beta diveristy

hja.bray <- 1 - (vegdist(OTUsREL, method = "bray"))
hja.bray.list <- liste(hja.bray, entry = "struc")
hja.elev.dist <- dist(env[,10])
hja.elev.dist.list <- liste(hja.elev.dist, entry = "elev")
hja.model <- lm(log(hja.bray.list$struc) ~ hja.elev.dist.list$elev)
summary(hja.model)

png(filename = "../figures/FigureS5.png", 
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 1, 3) + 0.1)
plot(log(hja.bray.list$struc) ~ hja.elev.dist.list$elev, 
     xlim = c(0, 800), xlab="",
     ylab = "", xaxt = "n", yaxt = "n")
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.2)
mtext("Difference in Elevation (m)", side = 1, line = 3, cex = 1.5)
box(lwd=2)
abline(hja.model, lwd = 2)
dev.off()
graphics.off()
img <- readPNG("../figures/FigureS5.png")
grid.raster(img)



# Elevation gradient analysis
water.bray <- 1 - (vegdist(OTUsREL.log[which(design$habitat=="water"),], method = "bray"))
water.bray.list <- liste(water.bray, entry = "struc")
water.elev.dist <- dist(env[which(design$habitat == "water"),10])
water.elev.dist.list <- liste(water.elev.dist, entry = "elev")
water.elev.lm <- lm(log(water.bray.list$struc) ~ water.elev.dist.list$elev)
summary(water.elev.lm)

sed.bray <- 1 - (vegdist(OTUsREL.log[which(design$habitat=="sediment"),], method = "bray"))
sed.bray.list <- liste(sed.bray, entry = "struc")
sed.elev.dist <- dist(env[which(design$habitat == "sediment"),10])
sed.elev.dist.list <- liste(sed.elev.dist, entry = "elev")
sed.elev.lm <- lm(log(sed.bray.list$struc) ~ sed.elev.dist.list$elev)
summary(sed.elev.lm)


### Figure S3: Bray-Curtis Similarity by Elevation
png(filename = "../figures/Figure10.png", 
    width = 1200, height = 1200, res = 96*2)

par(mfrow = c(2, 1))
par(mar = c(1, 5, 4, 3) + 0.3)

plot(log(water.bray.list$struc) ~ water.elev.dist.list$elev, 
     xlim = c(0,800), xlab="",
     ylab = "", xaxt = "n", yaxt = "n")
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
mtext("log Community Similarity\nSurface Water", side = 2, line = 3, cex = 1.2)
box(lwd=2)
abline(water.elev.lm, lty = 2, lwd = 2)

par(mar = c(5, 5, 1, 3) + 0.3)


plot(log(sed.bray.list$struc) ~ sed.elev.dist.list$elev,
     xlim = c(0,800), xlab="",
     ylab = "", xaxt = "n", yaxt = "n")
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
mtext("log Community Similarity\nSediments", side = 2, line = 3, cex = 1.2)
mtext("Elevation difference (m)", side = 1, line = 3, cex = 1.5)
box(lwd=2)
abline(sed.elev.lm, lwd = 2)


dev.off()
graphics.off()

img <- readPNG("../figures/Figure10.png")
grid.raster(img)
