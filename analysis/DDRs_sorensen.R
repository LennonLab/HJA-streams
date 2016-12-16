source("./analysis/InitialSetup.R")
source("./analysis/DistanceCalcs.R")
source("./analysis/Ordination.R")

# Water Distance Decay
water.xy <- xy[which(design$habitat == "water"),]
water.den.dist <- den.dists[which(design$habitat == "water"),which(design$habitat == "water")]
water.geo.dist.ls <- na.omit(liste(dist.mat[which(design$habitat == "water"),
                                            which(design$habitat == "water")], 
                                   entry = "geo.dist")[,3])
water.struc.dist <- 1 - water.dsorensen
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
capture.output(summary(water.env.lm), file = "./tables/DDR_water-env-sor.txt")
capture.output(summary(water.geo.lm), file = "./tables/DDR_water-space-sor.txt")

# Sediment Distance Decay
sed.xy <- xy[which(design$habitat == "sediment"),]
sed.den.dist <- den.dists[which(design$habitat == "sediment"),which(design$habitat == "sediment")]
sed.geo.dist.ls <- na.omit(liste(dist.mat[which(design$habitat == "sediment"),
                                          which(design$habitat == "sediment")], 
                                 entry = "geo.dist")[,3])
sed.struc.dist <- 1 - sediment.dsorensen
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
capture.output(summary(sed.env.lm), file = "./tables/DDR_sed-env-sor.txt")
capture.output(summary(sed.geo.lm), file = "./tables/DDR_sed-space-sor.txt")

# Headwaters DDRs
headwater.env.dists <- vegdist(env.2[which(design$order < 2),], method = "gower")
headwater.env.dists <- liste(headwater.env.dists, entry = "env")[,3]
headwater.geo.dists <- na.omit(liste(dist.mat[which(design$order < 2), which(design$order < 2)],
                              entry = "geo.dist"))[,3]
headwater.den.dist <- den.dists[which(design$order < 2),which(design$order < 2)]
headwater.den.dists <- na.omit(liste(headwater.den.dist, entry = "den.dist")[,3])
headwater.db <- vegdist(decostand(OTUsREL[which(design$order < 2),], method = "pa"))
headwater.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$order < 2),
                                                       which(design$order < 2)]),
                                   entry = "unifrac")[,3]

headwater.dists <- data.frame(liste((1 - headwater.db), entry = "comm.struc")[,3], 
                              headwater.env.dists, headwater.geo.dists, headwater.den.dists,
                              headwater.phylo.dist.ls)
colnames(headwater.dists) <- c("comm.struc", "env", "geo", "den", "unifrac")
headwater.env.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$env))
headwater.geo.lm <- (lm(log(headwater.dists$comm.struc) ~ headwater.dists$geo))
capture.output(summary(headwater.env.lm), file = "./tables/DDR_headwater-env-sor.txt")
capture.output(summary(headwater.geo.lm), file = "./tables/DDR_headwater-space-sor.txt")

# Higher Order DDRs
downstream.env.dists <- vegdist(env.2[which(design$order >= 2),], method = "gower")
downstream.env.dists <- liste(downstream.env.dists, entry = "env")[,3]
downstream.geo.dists <- na.omit(liste(dist.mat[which(design$order >= 2), which(design$order >= 2)],
                               entry = "geo.dist"))[,3]
downstream.db <- vegdist(decostand(OTUsREL[which(design$order >= 2),],method = "pa"))
downstream.den.dist <- den.dists[which(design$order > 1),which(design$order > 1)]
downstream.den.dists <- na.omit(liste(downstream.den.dist, entry = "den.dist")[,3])
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
capture.output(summary(downstream.env.lm), file = "./tables/DDR_downstream-env-sor.txt")
capture.output(summary(downstream.geo.lm), file = "./tables/DDR_downstream-space-sor.txt")

# Downstream seds
downstream.sed.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "sediment"),], method = "gower")
downstream.sed.env.dists <- liste(downstream.sed.env.dists, entry = "env")[,3]
downstream.sed.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "sediment"), 
                                                  which(design$order >1 & design$habitat == "sediment")],
                                     entry = "geo.dist"))[,3]
downstream.sed.db <- vegdist(decostand(OTUsREL[which(design$order >1 & design$habitat == "sediment"),], method = "pa"))
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

capture.output(summary(downstream.sed.env.lm), file = "./tables/DDR_downstream-seds-env-sor.txt")
capture.output(summary(downstream.sed.geo.lm), file = "./tables/DDR_downstream-seds-space-sor.txt")

plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$env)
plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo)

# Downstream water
downstream.water.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "water"),], method = "gower")
downstream.water.env.dists <- liste(downstream.water.env.dists, entry = "env")[,3]
downstream.water.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "water"), 
                                                   which(design$order >1 & design$habitat == "water")],
                                          entry = "geo.dist"))[,3]
downstream.water.db <- vegdist(decostand(OTUsREL[which(design$order >1 & design$habitat == "water"),], method = "pa"))
downstream.water.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$order >= 2 & design$habitat == "water"),
                                                              which(design$order >= 2 & design$habitat == "water")]),
                                          entry = "unifrac")[,3]
downstream.water.dists <- data.frame(liste((1 - downstream.water.db), entry = "comm.struc")[,3], 
                                     downstream.water.env.dists, downstream.water.geo.dists,
                                     downstream.water.phylo.dist.ls)
colnames(downstream.water.dists) <- c("comm.struc", "env", "geo", "unifrac")
downstream.water.env.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$env))
downstream.water.geo.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo))
capture.output(summary(downstream.water.env.lm), file = "./tables/DDR_downstream-water-env-sor.txt")
capture.output(summary(downstream.water.geo.lm), file = "./tables/DDR_downstream-water-space-sor.txt")

plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$env)
plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo.dist)

# # Downstream seds
# downstream.sed.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "sediment"),], method = "gower")
# downstream.sed.env.dists <- liste(downstream.sed.env.dists, entry = "env")[,3]
# downstream.sed.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "sediment"), 
#                                                    which(design$order >1 & design$habitat == "sediment")],
#                                           entry = "geo.dist"))
# downstream.sed.db <- vegdist(decostand(OTUsREL[which(design$order >1 & design$habitat == "sediment"),], method = "pa"))
# downstream.sed.dists <- cbind(liste((1 - downstream.sed.db), entry = "comm.struc"), 
#                               downstream.sed.env.dists, downstream.sed.geo.dists)
# downstream.sed.env.lm <- (lm(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$downstream.sed.env.dists))
# downstream.sed.geo.lm <- (lm(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo.dist))
# capture.output(summary(downstream.sed.env.lm), file = "./tables/DDR_downstream-seds-env-sor.txt")
# capture.output(summary(downstream.sed.geo.lm), file = "./tables/DDR_downstream-seds-space-sor.txt")
# 
# plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$downstream.sed.env.dists)
# plot(log(downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo.dist)
# 
# # Downstream water
# downstream.water.env.dists <- vegdist(env.2[which(design$order > 1 & design$habitat == "water"),], method = "gower")
# downstream.water.env.dists <- liste(downstream.water.env.dists, entry = "env")[,3]
# downstream.water.geo.dists <- na.omit(liste(dist.mat[which(design$order >1  & design$habitat == "water"), 
#                                                      which(design$order >1 & design$habitat == "water")],
#                                             entry = "geo.dist"))
# downstream.water.db <- vegdist(decostand(OTUsREL[which(design$order >1 & design$habitat == "water"),], method = "pa"))
# downstream.water.dists <- cbind(liste((1 - downstream.water.db), entry = "comm.struc"), 
#                                 downstream.water.env.dists, downstream.water.geo.dists)
# downstream.water.env.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$downstream.water.env.dists))
# downstream.water.geo.lm <- (lm(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo.dist))
# capture.output(summary(downstream.water.env.lm), file = "./tables/DDR_downstream-water-env-sor.txt")
# capture.output(summary(downstream.water.geo.lm), file = "./tables/DDR_downstream-water-space-sor.txt")
# 
# plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$downstream.water.env.dists)
# plot(log(downstream.water.dists$comm.struc) ~ downstream.water.dists$geo.dist)
# 

### Catchment Scale DDRs

hja.env.dists <- vegdist(env.2, method = "gower")
hja.env.dists <- liste(hja.env.dists, entry = "env")[,3]
hja.geo.dists <- na.omit(liste(dist.mat, entry = "geo.dist"))
hja.den.dists <- na.omit(liste(den.dists, entry = "dend.dist"))
hja.dists <- cbind(liste(hja.d.sorensen, entry = "comm.struc"), hja.env.dists)
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

png(filename = "./figures/DDR_HeadwaterDownstream_sorensen.png",
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
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.2)

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
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.2)
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
img <- readPNG("./figures/DDR_HeadwaterDownstream_sorensen.png")
grid.raster(img)





#### Figure: Water vs Sediment DDRs
png(filename = "./figures/DDR_WaterSed_sorensen.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(water.dists$env, 
     log(water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
abline(water.env.lm, lwd = 2)
r2 <- round(summary(water.env.lm)$r.squared, 3)
my.p <- round(summary(water.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

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
r2 <- round(summary(sed.env.lm)$r.squared, 3)
my.p <- round(summary(sed.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(water.dists$geo, 
     log(water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(water.geo.lm, lty = 1, lwd = 2)
r2 <- round(summary(water.geo.lm)$r.squared, 3)
my.p <- round(summary(water.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

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
r2 <- round(summary(sed.geo.lm)$r.squared, 3)
my.p <- round(summary(sed.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment Bacteria", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_WaterSed_sorensen.png")
grid.raster(img)





# Catchment-Scale DDRs
png(filename = "./figures/DDR_HJA_com-env_sorensen.png",
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
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.5)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_com-den_sorensen.png",
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
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.5)
mtext("Dendritic distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_com-geo_sorensen.png",
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
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.5)
mtext("Geographic Distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

# Downstream sed vs. water

png(filename = "./figures/DDR_DownstreamSedWater_sorensen.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(downstream.water.dists$env, 
     log(downstream.water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
abline(downstream.water.env.lm, lwd = 2)
r2 <- round(summary(downstream.water.env.lm)$r.squared, 3)
my.p <- round(summary(downstream.water.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Compositional Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.sed.dists$env, 
     log(downstream.sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.2)
abline(downstream.sed.env.lm, lwd = 2)
r2 <- round(summary(downstream.sed.env.lm)$r.squared, 3)
my.p <- round(summary(downstream.sed.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(downstream.water.dists$geo, 
     log(downstream.water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(downstream.water.geo.lm, lty = 1, lwd = 2)
r2 <- round(summary(downstream.water.geo.lm)$r.squared, 3)
my.p <- round(summary(downstream.water.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.sed.dists$geo, 
     log(downstream.sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(downstream.sed.geo.lm, lwd = 2)
r2 <- round(summary(downstream.sed.geo.lm)$r.squared, 3)
my.p <- round(summary(downstream.sed.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_DownstreamSedWater_sorensen.png")
grid.raster(img)


################################################################
png(filename = "./figures/DDR_HeadwaterDownstream_sor_dt.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

dt.headwater.env.lm<-(lm(residuals(headwater.geo.lm) + coefficients(summary(headwater.env.lm))[1,1] ~ headwater.dists$env))
dt.downstream.env.lm<-(lm(residuals(downstream.geo.lm) + coefficients(summary(downstream.env.lm))[1,1] ~ downstream.dists$env))
dt.headwater.geo.lm<-(lm(residuals(headwater.env.lm) + coefficients(summary(headwater.geo.lm))[1,1] ~ headwater.dists$geo))
dt.downstream.geo.lm<-(lm(residuals(downstream.env.lm) + coefficients(summary(downstream.geo.lm))[1,1] ~ downstream.dists$geo))

yrange.headwater <- c(
  min(residuals(headwater.geo.lm)+coefficients(summary(headwater.env.lm))[1,1],
      residuals(headwater.env.lm)+coefficients(summary(headwater.geo.lm))[1,1]),
  max(residuals(headwater.geo.lm)+coefficients(summary(headwater.env.lm))[1,1],
      residuals(headwater.env.lm)+coefficients(summary(headwater.geo.lm))[1,1]))
yrange.downstream <- c(
  min(residuals(downstream.env.lm)+coefficients(summary(downstream.geo.lm))[1,1],
      residuals(downstream.geo.lm)+coefficients(summary(downstream.env.lm))[1,1]),
  max(residuals(downstream.env.lm)+coefficients(summary(downstream.geo.lm))[1,1],
      residuals(downstream.geo.lm)+coefficients(summary(downstream.env.lm))[1,1]))


par(mar = c(1, 5, 3, 0) + 0.4)
plot(headwater.dists$env, 
     residuals(headwater.geo.lm)+ coefficients(summary(headwater.env.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,max(hja.env.dists)),
     ylim = yrange.headwater)
abline(dt.headwater.env.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.headwater.env.lm)$r.squared, 3)
my.p <- round(summary(dt.headwater.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n')

axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.dists$env, 
     residuals(downstream.geo.lm)+ coefficients(summary(downstream.env.lm))[1,1] , xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,max(hja.env.dists)),
     ylim = yrange.downstream)

abline(dt.downstream.env.lm, lwd=2)
r2 <- round(summary(dt.downstream.env.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')


axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(headwater.dists$geo, 
     residuals(headwater.env.lm)+ coefficients(summary(headwater.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000), ylim = yrange.headwater)
abline(dt.headwater.geo.lm, lty = 1, lwd = 2)

r2 <- round(summary(dt.headwater.geo.lm)$r.squared, 3)
my.p <- round(summary(dt.headwater.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n')


axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Headwaters", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.dists$geo, 
     residuals(downstream.env.lm) + coefficients(summary(downstream.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000), ylim = yrange.downstream)

abline(dt.downstream.geo.lm, lwd=2)
r2 <- round(summary(dt.downstream.geo.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Downstream", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_HeadwaterDownstream_sor_dt.png")
grid.raster(img)


#-------------------------------------------------------
png(filename = "./figures/DDR_WaterSed_sor_dt.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

#models
dt.water.env.lm<-(lm(residuals(water.geo.lm) + coefficients(summary(water.env.lm))[1,1] ~ water.dists$env))
dt.sed.env.lm<-(lm(residuals(sed.geo.lm) + coefficients(summary(sed.env.lm))[1,1] ~ sed.dists$env))
dt.water.geo.lm<-(lm(residuals(water.env.lm) + coefficients(summary(water.geo.lm))[1,1] ~ water.dists$geo))
dt.sed.geo.lm<-(lm(residuals(sed.env.lm) + coefficients(summary(sed.geo.lm))[1,1] ~ sed.dists$geo))

yrange.water <- c(
  min(residuals(water.geo.lm)+coefficients(summary(water.env.lm))[1,1],
      residuals(water.env.lm)+coefficients(summary(water.geo.lm))[1,1]),
  max(residuals(water.geo.lm)+coefficients(summary(water.env.lm))[1,1],
      residuals(water.env.lm)+coefficients(summary(water.geo.lm))[1,1]))
yrange.sed <- c(
  min(residuals(sed.env.lm)+coefficients(summary(sed.geo.lm))[1,1],
      residuals(sed.geo.lm)+coefficients(summary(sed.env.lm))[1,1]),
  max(residuals(sed.env.lm)+coefficients(summary(sed.geo.lm))[1,1],
      residuals(sed.geo.lm)+coefficients(summary(sed.env.lm))[1,1]))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(water.dists$env, 
     residuals(water.geo.lm)+coefficients(summary(water.env.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6),
     ylim = yrange.water)
abline(dt.water.env.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.water.env.lm)$r.squared, 3)
my.p <- round(summary(dt.water.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n')

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
     residuals(sed.geo.lm)+ coefficients(summary(sed.env.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6),
     ylim = yrange.sed)
abline(dt.sed.env.lm, lty = 1, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
r2 <- round(summary(dt.sed.env.lm)$r.squared, 3)
my.p <- round(summary(dt.sed.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

par(mar = c(1, 1, 3, 4) + 0.4)
plot(water.dists$geo, 
     residuals(water.env.lm) + coefficients(summary(water.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000), ylim = yrange.water)
abline(dt.water.geo.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.water.geo.lm)$r.squared, 3)
my.p <- round(summary(dt.water.geo.lm)$coefficients[2,4],3)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n')

axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(sed.dists$geo, 
     residuals(sed.env.lm)+ coefficients(summary(sed.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000), ylim = yrange.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(dt.sed.geo.lm, lwd = 2)
r2 <- round(summary(dt.sed.geo.lm)$r.squared, 3)
my.p <- round(summary(dt.sed.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment Bacteria", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_WaterSed_sor_dt.png")
grid.raster(img)



# Downstream sed vs. water
png(filename = "./figures/DDR_DownstreamSedWater_sor_dt.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

dt.downstream.water.env.lm<-(lm(residuals(downstream.water.geo.lm) + coefficients(summary(downstream.water.env.lm))[1,1] ~ downstream.water.dists$env))
dt.downstream.sed.env.lm<-(lm(residuals(downstream.sed.geo.lm) + coefficients(summary(downstream.sed.env.lm))[1,1] ~ downstream.sed.dists$env))
dt.downstream.water.geo.lm<-(lm(residuals(downstream.water.env.lm) + coefficients(summary(downstream.water.geo.lm))[1,1] ~ downstream.water.dists$geo))
dt.downstream.sed.geo.lm<-(lm(residuals(downstream.sed.env.lm) + coefficients(summary(downstream.sed.geo.lm))[1,1] ~ downstream.sed.dists$geo))

yrange.ds.water <- c(
  min(residuals(downstream.water.geo.lm)+coefficients(summary(downstream.water.env.lm))[1,1],
      residuals(downstream.water.env.lm)+coefficients(summary(downstream.water.geo.lm))[1,1]),
  max(residuals(downstream.water.geo.lm)+coefficients(summary(downstream.water.env.lm))[1,1],
      residuals(downstream.water.env.lm)+coefficients(summary(downstream.water.geo.lm))[1,1]))
yrange.ds.sed <- c(
  min(residuals(downstream.sed.env.lm)+coefficients(summary(downstream.sed.geo.lm))[1,1],
      residuals(downstream.sed.geo.lm)+coefficients(summary(downstream.sed.env.lm))[1,1]),
  max(residuals(downstream.sed.env.lm)+coefficients(summary(downstream.sed.geo.lm))[1,1],
      residuals(downstream.sed.geo.lm)+coefficients(summary(downstream.sed.env.lm))[1,1]))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(downstream.water.dists$env, 
     residuals(downstream.water.geo.lm)+ coefficients(summary(downstream.water.env.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6), ylim = yrange.ds.water)
abline(dt.downstream.water.env.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.downstream.water.env.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.water.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.sed.dists$env, 
     residuals(downstream.sed.geo.lm) + coefficients(summary(downstream.sed.env.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6), ylim = yrange.ds.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
abline(dt.downstream.sed.env.lm, lwd = 2)
r2 <- round(summary(dt.downstream.sed.env.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.sed.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(downstream.water.dists$geo, 
     residuals(downstream.water.env.lm)+ coefficients(summary(downstream.water.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000), ylim = yrange.ds.water)
abline(dt.downstream.water.geo.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.downstream.water.geo.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.water.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.sed.dists$geo, 
     residuals(downstream.sed.env.lm)+ coefficients(summary(downstream.sed.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000), ylim = yrange.ds.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=log(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(dt.downstream.sed.geo.lm, lwd = 2)
r2 <- round(summary(dt.downstream.sed.geo.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.sed.geo.lm)$coefficients[2,4],3)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Sediment", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
img <- readPNG("./figures/DDR_DownstreamSedWater_sor_dt.png")
grid.raster(img)

