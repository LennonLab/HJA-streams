# source("./InitialSetup.R")
# source("./DistanceCalcs.R")

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

water.lm <- (lm(log(water.dists$comm.struc) ~ water.dists$env))
water.dd <- (lm(log(water.dists$comm.struc) ~ water.dists$geo))
summary(water.lm)
summary(water.dd)

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

sed.lm <- (lm(log(sed.dists$comm.struc) ~ sed.dists$env))
sed.dd <- (lm(log(sed.dists$comm.struc) ~ sed.dists$geo))
summary(sed.lm)
summary(sed.dd)

# Headwaters DDRs
headwater.env.dists <- vegdist(env.2[which(design$order < 2),], method = "bray")
headwater.env.dists <- liste(headwater.env.dists, entry = "env")[,3]
headwater.geo.dists <- na.omit(liste(dist.mat[which(design$order < 2), which(design$order < 2)],
                              entry = "geo.dist"))
headwater.db <- vegdist(OTUsREL[which(design$order < 2),])
headwater.dists <- cbind(liste((1 - headwater.db), entry = "comm.struc"), 
                         headwater.env.dists, headwater.geo.dists)
headwater.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$headwater.env.dists))
headwater.dd <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$geo.dist))
summary(headwater.lm)
summary(headwater.dd)

# Higher Order DDRs
mainstem.env.dists <- vegdist(env.2[which(design$order >= 2),], method = "bray")
mainstem.env.dists <- liste(mainstem.env.dists, entry = "env")[,3]
mainstem.geo.dists <- na.omit(liste(dist.mat[which(design$order >= 2), which(design$order >= 2)],
                               entry = "geo.dist"))
mainstem.db <- vegdist(OTUsREL[which(design$order >= 2),])
mainstem.dists <- cbind(liste((1 - mainstem.db), entry = "comm.struc"), 
                        mainstem.env.dists, mainstem.geo.dists)
mainstem.lm <- (lm(log(mainstem.dists$comm.struc) ~ mainstem.dists$mainstem.env.dists))
mainstem.dd <- (lm(log(mainstem.dists$comm.struc) ~ mainstem.dists$geo.dist))
summary(mainstem.lm)
summary(mainstem.dd)


### Catchment Scale DDRs

# # Principal Components Analysis on HJA Environment
# hja.pca <- princomp(env.mat)
# summary(hja.pca)
# plot(hja.pca, type = "l")
# biplot(hja.pca)
# pc1 <- hja.pca$scores[,1]
# pc2 <- hja.pca$scores[,2]
# pc3 <- hja.pca$scores[,3]
# 
# cor(pc1, env.mat)
# cor(pc2, env.mat)

env.dists <- vegdist(env.2, method = "gower")
env.dists <- liste(env.dists, entry = "env")[,3]
geo.dists <- na.omit(liste(dist.mat, entry = "geo.dist"))
#env.dists <- liste(dist(pc1), entry = "env")[,3]
den.dists <- na.omit(liste(den.dists, entry = "dend.dist"))
hja.dists <- cbind(liste(hja.db, entry = "comm.struc"), env.dists)
hja.dists <- cbind(hja.dists, den.dists[,3], geo.dists[,3])
names(hja.dists)[c(5,6)] <- c("den.dists", "geo.dists")

hja.env.lm <- (lm(log(1-hja.dists$comm.struc) ~ hja.dists$env.dists))
hja.den.lm <- (lm(log(1-hja.dists$comm.struc) ~ hja.dists$den.dists))
hja.geo.lm <- (lm(log(1-hja.dists$comm.struc) ~ hja.dists$geo.dists))
hja.geo_env.lm <- lm(log(1-hja.dists$comm.struc) ~ hja.dists$geo.dists * hja.dists$env.dists)

summary(hja.geo_env.lm)
summary(hja.env.lm)
summary(hja.geo.lm)
summary(hja.den.lm)

# water.pca <- princomp(env.mat[which(design$habitat == "water"),2:6])
# summary(water.pca)
# water.pc1 <- water.pca$scores[,1]
# cor(water.pc1, env.mat[which(design$habitat == "water"),2:6])
# water.env.dists <- vegdist(env.2[which(design$habitat == "water"),], method = "gower")
# water.env.dists <- liste(water.env.dists, entry = "env")[,3]
# water.geo.dists <- na.omit(liste(dist.mat[which(design$habitat == "water"), which(design$habitat == "water")],
#                                  entry = "geo.dist"))
# #water.env.dists <- liste(dist(water.pc1), entry = "env")[,3]
# water.dists <- cbind(liste((1 - water.db), entry = "comm.struc"), water.env.dists, water.geo.dists)
# plot(water.dists$water.env.dists, log(water.dists$comm.struc))

# sed.env.dists <- vegdist(env.2[which(design$habitat == "sediment"),], method = "gower")
# sed.env.dists <- liste(sed.env.dists, entry = "env")[,3]
# sed.pca <- princomp(env.mat[which(design$habitat == "sediment"),2:6])
# summary(sed.pca)
# sed.pc1 <- sed.pca$scores[,3]
# cor(sed.pc1, env.mat[which(design$habitat == "sediment"),2:6])
# sed.geo.dists <- na.omit(liste(dist.mat[which(design$habitat == "sediment"), which(design$habitat == "sediment")],
#                                entry = "geo.dist"))
# #sed.env.dists <- liste(dist(sed.pc1), entry = "env")[,3]
# sed.dists <- cbind(liste((1-sediment.db), entry = "comm.struc"), sed.env.dists, sed.geo.dists)
# plot(sed.dists$sed.env.dists, log(sed.dists$comm.struc))


plot(sed.dists$geo.dist, log(sed.dists$comm.struc))
plot(water.dists$geo.dist, log(water.dists$comm.struc))

summary(sed.dd)
summary(water.dd)

################# FIGURES

##### Figure: Headwater vs. Mainstem DDRs

png(filename = "../figures/DDR_HeadwaterMainstem.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(headwater.dists$headwater.env.dists, 
     log(headwater.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.5))
abline(headwater.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(mainstem.dists$mainstem.env.dists, 
     log(mainstem.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.5))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.2)
abline(mainstem.lm, lwd = 2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(headwater.dists$geo.dist, 
     log(headwater.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(headwater.dd, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Headwaters", side = 4, line = 1, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(mainstem.dists$geo.dist, 
     log(mainstem.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(mainstem.dd, lwd = 2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Higher Orders", side = 4, line = 1, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("../figures/DDR_HeadwaterMainstem.png")
grid.raster(img)





#### Figure: Water vs Sediment DDRs
png(filename = "../figures/DDR_WaterSed.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(water.dists$env, 
     log(water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
abline(water.lm, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(sed.dists$env, 
     log(sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.2)
abline(sed.lm, lwd = 2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(water.dists$geo, 
     log(water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(water.dd, lty = 1, lwd = 2)
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Waters", side = 4, line = 1, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(sed.dists$geo, 
     log(sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(sed.dd, lwd = 2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment", side = 4, line = 1, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("../figures/DDR_WaterSed.png")
grid.raster(img)





# Catchment-Scale DDRs
png(filename = "../figures/Figure10.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$env.dists, log(1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.env.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.5)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "../figures/Figure11.png",
    width = 1200, height = 1200, res = 96*2)
grid.raster(readPNG("../figures/Figure10.png"))

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$den.dists, log(1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.den.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.5)
mtext("Dendritic distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "../figures/Figure12.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$geo.dists, log(1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.geo.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("log Community Similarity", side = 2, line = 3, cex = 1.5)
mtext("Geographic Distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "../figures/Figure13.png",
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

png(filename = "../figures/Figure14.png",
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

