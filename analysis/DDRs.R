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
water.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$habitat == "water"),
                                         which(design$habitat == "water")]),
                             entry = "unifrac")[,3]
water.dists <- data.frame(water.den.dist.ls, water.geo.dist.ls, 
                          water.struc.dist.ls, water.env.dist.ls,
                          water.phylo.dist.ls)
names(water.dists) <- c("den", "geo", "comm.struc", "env", "unifrac")

water.env.lm <- (lm(log(water.dists$comm.struc) ~ water.dists$env))
water.geo.lm <- (lm(log(water.dists$comm.struc) ~ water.dists$geo))
summary(lm(log(water.dists$unifrac) ~ water.dists$env))
plot(water.dists$unifrac ~ water.dists$geo)
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
sed.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$habitat == "sediment"),
                                                 which(design$habitat == "sediment")]),
                             entry = "unifrac")[,3]
sed.dists <- data.frame(sed.den.dist.ls, sed.geo.dist.ls, 
                          sed.struc.dist.ls, sed.env.dist.ls,
                          sed.phylo.dist.ls)
names(sed.dists) <- c("den", "geo", "comm.struc", "env", "unifrac")

sed.env.lm <- (lm(log(sed.dists$comm.struc) ~ sed.dists$env))
sed.geo.lm <- (lm(log(sed.dists$comm.struc) ~ sed.dists$geo))
plot(sed.dists$unifrac ~ sed.dists$env)
plot(log10(sed.dists$comm.struc) ~ sed.dists$unifrac)
capture.output(summary(sed.env.lm), file = "./tables/DDR_sed-env.txt")
capture.output(summary(sed.geo.lm), file = "./tables/DDR_sed-space.txt")

# Headwaters DDRs
headwater.env.dists <- vegdist(env.2[which(design$order < 2),], method = "gower")
headwater.env.dists <- liste(headwater.env.dists, entry = "env")[,3]
headwater.den.dist <- den.dists[which(design$order < 2),which(design$order < 2)]
headwater.den.dists <- na.omit(liste(headwater.den.dist, entry = "den.dist")[,3])
headwater.geo.dists <- na.omit(liste(dist.mat[which(design$order < 2), which(design$order < 2)],
                              entry = "geo.dist")[,3])
headwater.db <- vegdist(OTUsREL[which(design$order < 2),])
headwater.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$order < 2),
                                                 which(design$order < 2)]),
                             entry = "unifrac")[,3]

headwater.dists <- data.frame(liste((1 - headwater.db), entry = "comm.struc")[,3], 
                         headwater.env.dists, headwater.geo.dists, headwater.den.dists,
                         headwater.phylo.dist.ls)
colnames(headwater.dists) <- c("comm.struc", "env", "geo", "den", "unifrac")
headwater.env.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$env))
headwater.geo.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$geo))
headwater.den.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$den))
headwater.phy.env.lm <- (lm(log(headwater.dists$unifrac) ~ headwater.dists$env))
headwater.phy.geo.lm <- (lm(log(headwater.dists$unifrac) ~ headwater.dists$geo))
capture.output(summary(headwater.env.lm), file = "./tables/DDR_headwater-env.txt")
capture.output(summary(headwater.geo.lm), file = "./tables/DDR_headwater-space.txt")

# Downstream DDRs
downstream.env.dists <- vegdist(env.2[which(design$order >= 2),], method = "gower")
downstream.env.dists <- liste(downstream.env.dists, entry = "env")[,3]
downstream.geo.dists <- na.omit(liste(dist.mat[which(design$order >= 2), which(design$order >= 2)],
                               entry = "geo.dist"))[,3]
downstream.den.dist <- den.dists[which(design$order > 1),which(design$order > 1)]
downstream.den.dists <- na.omit(liste(downstream.den.dist, entry = "den.dist")[,3])

downstream.db <- vegdist(OTUsREL[which(design$order >= 2),])
downstream.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$order >= 2),
                                                     which(design$order >= 2)]),
                                 entry = "unifrac")[,3]
downstream.dists <- data.frame(liste((1 - downstream.db), entry = "comm.struc")[,3], 
                        downstream.env.dists, downstream.geo.dists, downstream.den.dists,
                        downstream.phylo.dist.ls)
colnames(downstream.dists) <- c("comm.struc", "env", "geo", "den", "unifrac")
downstream.env.lm <- (lm(log(downstream.dists$comm.struc) ~ downstream.dists$env))
downstream.geo.lm <- (lm(log(downstream.dists$comm.struc) ~ downstream.dists$geo))
downstream.den.lm <- (lm(log(downstream.dists$comm.struc) ~ downstream.dists$den))
downstream.phy.env.lm <- lm(log(downstream.dists$unifrac) ~ downstream.dists$env)
downstream.phy.geo.lm <- lm(log(downstream.dists$unifrac) ~ downstream.dists$geo)
capture.output(summary(downstream.env.lm), file = "./tables/DDR_downstream-env.txt")
capture.output(summary(downstream.geo.lm), file = "./tables/DDR_downstream-space.txt")

# Downstream seds
downstream.sed.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "sediment"),], method = "gower")
downstream.sed.env.dists <- liste(downstream.sed.env.dists, entry = "env")[,3]
downstream.sed.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "sediment"), 
                                                  which(design$order >1 & design$habitat == "sediment")],
                                     entry = "geo.dist"))[,3]
downstream.sed.db <- vegdist(OTUsREL[which(design$order >1 & design$habitat == "sediment"),])
downstream.sed.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$order >= 2 & design$habitat == "sediment"),
                                                      which(design$order >= 2 & design$habitat == "sediment")]),
                                  entry = "unifrac")[,3]
downstream.sed.dists <- data.frame(liste((1 - downstream.sed.db), entry = "comm.struc")[,3], 
                         downstream.sed.env.dists, downstream.sed.geo.dists,
                         downstream.sed.phylo.dist.ls)
colnames(downstream.sed.dists) <- c("comm.struc", "env", "geo", "unifrac")
downstream.sed.env.lm <- (lm(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$env))
downstream.sed.geo.lm <- (lm(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo))
downstream.sed.phy.env.lm <- (lm(log(downstream.sed.dists$unifrac) ~ downstream.sed.dists$env))
downstream.sed.phy.geo.lm <- (lm(log(downstream.sed.dists$unifrac) ~ downstream.sed.dists$geo))
summary(downstream.sed.env.lm)
summary(downstream.sed.phy.env.lm)
summary(downstream.sed.geo.lm)
summary(downstream.sed.phy.geo.lm)


plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$env)
abline(downstream.sed.env.lm)
plot(log(downstream.sed.dists$unifrac) ~ downstream.sed.dists$env)
abline(downstream.sed.phy.env.lm)
plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo)
abline(downstream.sed.geo.lm)
plot(log(downstream.sed.dists$unifrac) ~ downstream.sed.dists$geo)
abline(downstream.sed.phy.geo.lm)

# Downstream water
downstream.water.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "water"),], method = "gower")
downstream.water.env.dists <- liste(downstream.water.env.dists, entry = "env")[,3]
downstream.water.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "water"), 
                                                   which(design$order >1 & design$habitat == "water")],
                                          entry = "geo.dist"))[,3]
downstream.water.db <- vegdist(OTUsREL[which(design$order >1 & design$habitat == "water"),])
downstream.water.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$order >= 2 & design$habitat == "water"),
                                                          which(design$order >= 2 & design$habitat == "water")]),
                                      entry = "unifrac")[,3]
downstream.water.dists <- data.frame(liste((1 - downstream.water.db), entry = "comm.struc")[,3], 
                              downstream.water.env.dists, downstream.water.geo.dists,
                              downstream.water.phylo.dist.ls)
colnames(downstream.water.dists) <- c("comm.struc", "env", "geo", "unifrac")
downstream.water.env.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$env))
downstream.water.geo.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo))
downstream.water.phy.env.lm <- lm(log(downstream.water.dists$unifrac) ~ downstream.water.dists$env)
downstream.water.phy.geo.lm <- lm(log(downstream.water.dists$unifrac) ~ downstream.water.dists$geo)
summary(downstream.water.env.lm)
summary(downstream.water.phy.env.lm)
summary(downstream.water.geo.lm)
summary(downstream.water.phy.geo.lm)
plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$env)
abline(downstream.water.env.lm)
plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo)
abline(downstream.water.geo.lm)
plot(log(downstream.water.dists$unifrac) ~ downstream.water.dists$env)
abline(downstream.water.phy.env.lm)
plot(log(downstream.water.dists$unifrac) ~ downstream.water.dists$geo)
abline(downstream.water.phy.geo.lm)

### Catchment Scale DDRs

hja.env.dists <- vegdist(env.2, method = "gower")
hja.env.dists <- liste(hja.env.dists, entry = "env")[,3]
hja.geo.dists <- na.omit(liste(dist.mat, entry = "geo"))[,3]
hja.den.dists <- na.omit(liste(den.dists, entry = "den"))[,3]
hja.com.dists <- 1-liste(hja.db, entry = "comm.struc")[,3]
hja.phy.dists <- 1-liste(hja.unifrac.dist, entry = "unifrac")[,3]
hja.dists <- data.frame(hja.com.dists, hja.env.dists, hja.den.dists, hja.geo.dists, hja.phy.dists)
colnames(hja.dists) <- c("comm.struc", "env", "den", "geo", "unifrac")

hja.env.lm <- (lm(log(hja.dists$comm.struc) ~ hja.dists$env))
hja.den.lm <- (lm(log(hja.dists$comm.struc) ~ hja.dists$den))
hja.geo.lm <- (lm(log(hja.dists$comm.struc) ~ hja.dists$geo))
hja.geo_env.lm <- lm(log(hja.dists$comm.struc) ~ hja.dists$geo * hja.dists$env)

AIC(hja.geo_env.lm)
AIC(hja.env.lm)
AIC(hja.geo.lm)

summary(hja.geo_env.lm)
summary(hja.env.lm)
summary(hja.geo.lm)
summary(hja.den.lm)
plot(log(hja.dists$comm.struc) ~ hja.dists$geo)
plot(log(hja.dists$comm.struc) ~ hja.dists$env)
plot(log(hja.dists$unifrac) ~ hja.dists$env)
plot(log(hja.dists$unifrac) ~ hja.dists$geo)

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

png(filename = "./figures/DDR_HJA_geo-den.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$den.dists, hja.dists$geo.dists, xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(lm(hja.dists$env.dists ~ hja.dists$den.dists), lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Geo Distance", side = 2, line = 3, cex = 1.5)
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



