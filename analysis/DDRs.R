# source("./analysis/InitialSetup.R")
# source("./analysis/DistanceCalcs.R")
# source("./analysis/Ordination.R")

# Water Distance Decay
water.xy <- xy[which(design$habitat == "water"),]
water.den.dist <- den.dists[which(design$habitat == "water"),which(design$habitat == "water")]
water.geo.dist.ls <- na.omit(liste(dist.mat[which(design$habitat == "water"),
                                            which(design$habitat == "water")], 
                                   entry = "geo.dist")[,3])
water.struc.dist <- 1 - water.db
water.env.dist <- vegdist(env.2[which(design$habitat == "water"),], 
                          method = "gower", upper = F, diag = F)
water.env.dist.ls <- liste(water.env.dist, entry = "env")[,3]
water.struc.dist.ls <- liste(water.struc.dist, entry = "struc")[,3]
water.den.dist.ls <- na.omit(liste(water.den.dist, entry = "den.dist")[,3])
water.dists <- data.frame(water.den.dist.ls, water.geo.dist.ls, 
                              water.struc.dist.ls, water.env.dist.ls)
names(water.dists) <- c("den", "geo", "comm.struc", "env")

water.env.lm <- (lm(log(water.dists$comm.struc) ~ water.dists$env))
water.geo.lm <- (lm(log(water.dists$comm.struc) ~ water.dists$geo))
capture.output(summary(water.env.lm), file = "./tables/DDR_water-env.txt")
capture.output(summary(water.geo.lm), file = "./tables/DDR_water-space.txt")

# Sediment Distance Decay
sed.xy <- xy[which(design$habitat == "sediment"),]
sed.den.dist <- den.dists[which(design$habitat == "sediment"),which(design$habitat == "sediment")]
sed.geo.dist.ls <- na.omit(liste(dist.mat[which(design$habitat == "sediment"),
                                          which(design$habitat == "sediment")], 
                                 entry = "geo.dist")[,3])
sed.struc.dist <- 1 - sediment.db
sed.env.dist <- vegdist(env.2[which(design$habitat == "sediment"),], 
                        method = "gower", upper = F, diag = F)
sed.env.dist.ls <- liste(sed.env.dist, entry = "env")[,3]
sed.struc.dist.ls <- liste(sed.struc.dist, entry = "struc")[,3]
sed.den.dist.ls <- na.omit(liste(sed.den.dist, entry = "den.dist")[,3])
sed.dists <- data.frame(sed.den.dist.ls, sed.geo.dist.ls, 
                            sed.struc.dist.ls, sed.env.dist.ls)
names(sed.dists) <- c("den", "geo", "comm.struc", "env")

sed.env.lm <- (lm(log(sed.dists$comm.struc) ~ sed.dists$env))
sed.geo.lm <- (lm(log(sed.dists$comm.struc) ~ sed.dists$geo))
capture.output(summary(sed.env.lm), file = "./tables/DDR_sed-env.txt")
capture.output(summary(sed.geo.lm), file = "./tables/DDR_sed-space.txt")

# Headwaters DDRs
headwater.env.dists <- vegdist(env.2[which(design$order < 2),], method = "gower")
headwater.env.dists <- liste(headwater.env.dists, entry = "env")[,3]
headwater.geo.dists <- na.omit(liste(dist.mat[which(design$order < 2), which(design$order < 2)],
                              entry = "geo.dist"))
headwater.db <- vegdist(OTUsREL[which(design$order < 2),])
headwater.dists <- cbind(liste((1 - headwater.db), entry = "comm.struc"), 
                         headwater.env.dists, headwater.geo.dists)
headwater.env.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$headwater.env.dists))
headwater.geo.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$geo.dist))
capture.output(summary(headwater.env.lm), file = "./tables/DDR_headwater-env.txt")
capture.output(summary(headwater.geo.lm), file = "./tables/DDR_headwater-space.txt")

# Higher Order DDRs
downstream.env.dists <- vegdist(env.2[which(design$order >= 2),], method = "gower")
downstream.env.dists <- liste(downstream.env.dists, entry = "env")[,3]
downstream.geo.dists <- na.omit(liste(dist.mat[which(design$order >= 2), which(design$order >= 2)],
                               entry = "geo.dist"))
downstream.db <- vegdist(OTUsREL[which(design$order >= 2),])
downstream.dists <- cbind(liste((1 - downstream.db), entry = "comm.struc"), 
                        downstream.env.dists, downstream.geo.dists)
downstream.env.lm <- (lm(log(downstream.dists$comm.struc) ~ downstream.dists$downstream.env.dists))
downstream.geo.lm <- (lm(log(downstream.dists$comm.struc) ~ downstream.dists$geo.dist))
capture.output(summary(downstream.env.lm), file = "./tables/DDR_downstream-env.txt")
capture.output(summary(downstream.geo.lm), file = "./tables/DDR_downstream-space.txt")

# Downstream seds
downstream.sed.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "sediment"),], method = "gower")
downstream.sed.env.dists <- liste(downstream.sed.env.dists, entry = "env")[,3]
downstream.sed.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "sediment"), 
                                                  which(design$order >1 & design$habitat == "sediment")],
                                     entry = "geo.dist"))
downstream.sed.db <- vegdist(OTUsREL[which(design$order >1 & design$habitat == "sediment"),])
downstream.sed.dists <- cbind(liste((1 - downstream.sed.db), entry = "comm.struc"), 
                         downstream.sed.env.dists, downstream.sed.geo.dists)
downstream.sed.env.lm <- (lm(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$downstream.sed.env.dists))
downstream.sed.geo.lm <- (lm(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo.dist))
summary(downstream.sed.env.lm)
summary(downstream.sed.geo.lm)

plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$downstream.sed.env.dists)
plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo.dist)

# Downstream water
downstream.water.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "water"),], method = "gower")
downstream.water.env.dists <- liste(downstream.water.env.dists, entry = "env")[,3]
downstream.water.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "water"), 
                                                   which(design$order >1 & design$habitat == "water")],
                                          entry = "geo.dist"))
downstream.water.db <- vegdist(OTUsREL[which(design$order >1 & design$habitat == "water"),])
downstream.water.dists <- cbind(liste((1 - downstream.water.db), entry = "comm.struc"), 
                              downstream.water.env.dists, downstream.water.geo.dists)
downstream.water.env.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$downstream.water.env.dists))
downstream.water.geo.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo.dist))
summary(downstream.water.env.lm)
summary(downstream.water.geo.lm)

plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$downstream.water.env.dists)
plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo.dist)

### Catchment Scale DDRs

hja.env.dists <- vegdist(env.2, method = "gower")
hja.env.dists <- liste(hja.env.dists, entry = "env")[,3]
hja.geo.dists <- na.omit(liste(dist.mat, entry = "geo.dist"))
hja.den.dists <- na.omit(liste(den.dists, entry = "dend.dist"))
hja.dists <- cbind(liste(hja.db, entry = "comm.struc"), hja.env.dists)
hja.dists <- cbind(hja.dists, hja.den.dists[,3], hja.geo.dists[,3])
names(hja.dists)[c(4,5,6)] <- c("env.dists", "den.dists", "geo.dists")

hja.env.lm <- (lm(log(1-hja.dists$comm.struc) ~ hja.dists$env.dists))
hja.den.lm <- (lm(log(1-hja.dists$comm.struc) ~ hja.dists$den.dists))
hja.geo.lm <- (lm(log(1-hja.dists$comm.struc) ~ hja.dists$geo.dists))
hja.geo_env.lm <- lm(log(1-hja.dists$comm.struc) ~ hja.dists$geo.dists * hja.dists$env.dists)

AIC(hja.geo_env.lm)
AIC(hja.env.lm)
AIC(hja.geo.lm)

summary(hja.geo_env.lm)
summary(hja.env.lm)
summary(hja.geo.lm)
summary(hja.den.lm)

################# FIGURES

##### Figure: Headwater vs. Mainstem DDRs

png(filename = "./figures/DDR_HeadwaterDownstream.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(headwater.dists$headwater.env.dists, 
     log(headwater.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.8))
# abline(headwater.env.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.dists$downstream.env.dists, 
     log(downstream.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.8))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
abline(downstream.env.lm, lwd = 2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(headwater.dists$geo.dist, 
     log(headwater.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
# abline(headwater.geo.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Headwaters", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.dists$geo.dist, 
     log(downstream.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(downstream.geo.lm, lwd = 2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Downstream", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_HeadwaterDownstream.png")
grid.raster(img)





#### Figure: Water vs Sediment DDRs
png(filename = "./figures/DDR_WaterSed.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(water.dists$env, 
     log(water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
abline(water.env.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(sed.dists$env, 
     log(sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
abline(sed.env.lm, lwd = 2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(water.dists$geo, 
     log(water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(water.geo.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(sed.dists$geo, 
     log(sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(sed.geo.lm, lwd = 2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment Bacteria", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_WaterSed.png")
grid.raster(img)





# Catchment-Scale DDRs
png(filename = "./figures/DDR_HJA_com-env.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$env.dists, log(1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.env.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.5)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_com-den.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$den.dists, log(1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.den.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.5)
mtext("Dendritic distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_com-geo.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$geo.dists, log(1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.geo.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.5)
mtext("Geographic Distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_env-geo.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$geo.dists, hja.dists$env.dists, xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(lm(hja.dists$env.dists ~ hja.dists$geo.dists), lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Environmental Distance", side = 2, line = 3, cex = 1.5)
mtext("Geographic Distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_env-den.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$den.dists, hja.dists$env.dists, xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(lm(hja.dists$env.dists ~ hja.dists$den.dists), lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Environmental Distance", side = 2, line = 3, cex = 1.5)
mtext("Dendritic Distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()


# Downstream sed vs. water

png(filename = "./figures/DDR_DownstreamSedWater.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(downstream.water.dists$downstream.water.env.dists, 
     log(downstream.water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
abline(downstream.water.env.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.sed.dists$downstream.sed.env.dists, 
     log(downstream.sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
abline(downstream.sed.env.lm, lwd = 2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(downstream.water.dists$geo.dist, 
     log(downstream.water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(downstream.water.geo.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.sed.dists$geo.dist, 
     log(downstream.sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(downstream.sed.geo.lm, lwd = 2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_DownstreamSedWater.png")
grid.raster(img)



