source("DDRs.R")

################# FIGURES

##### Figure: Headwater vs. Mainstem DDRs

png(filename = "./figures/DDR_HeadwaterDownstream_bray.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(headwater.dists$env, 
     (headwater.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.8))
abline(headwater.env.lm, lty = 1, lwd = 2)
r2 <- round(summary(headwater.env.lm)$r.squared, 3)
my.p <- round(summary(headwater.env.lm)$coefficients[2,4],3)
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.dists$env, 
     (downstream.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.8))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.2)
abline(downstream.env.lm, lwd = 2)
r2 <- round(summary(downstream.env.lm)$r.squared, 3)
my.p <- round(summary(downstream.env.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(headwater.dists$geo, 
     (headwater.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,12000))
abline(headwater.geo.lm, lty = 1, lwd = 2)
r2 <- round(summary(headwater.geo.lm)$r.squared, 3)
my.p <- round(summary(headwater.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Headwaters", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.dists$geo, 
     (downstream.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(downstream.geo.lm, lwd = 2)
r2 <- round(summary(downstream.geo.lm)$r.squared, 3)
my.p <- round(summary(downstream.geo.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Downstream", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
grid.raster(readPNG("./figures/DDR_HeadwaterDownstream_bray.png"))






#### Figure: Water vs Sediment DDRs
png(filename = "./figures/DDR_WaterSed_bray.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(water.dists$env, 
     (water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
     (sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
     (water.dists$comm.struc), xlab="", 
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
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(sed.dists$geo, 
     (sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
grid.raster(readPNG("./figures/DDR_WaterSed_bray.png"))





# Catchment-Scale DDRs
png(filename = "./figures/DDR_HJA_com-env_bray.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$env, (1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.env.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.5)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_com-den_bray.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$den, (1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.den.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.5)
mtext("Dendritic distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DDR_HJA_com-geo_bray.png",
    width = 1200, height = 1200, res = 96*2)

par(mar = c(5, 5, 3, 3) + 0.4)
plot(hja.dists$geo, (1 - hja.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n")
abline(hja.geo.lm, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Compositional Similarity", side = 2, line = 3, cex = 1.5)
mtext("Geographic Distance (m)", side = 1, line = 3, cex = 1.5)
dev.off()
graphics.off()

# Downstream sed vs. water

png(filename = "./figures/DDR_DownstreamSedWater_geo.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(downstream.water.dists$env, 
     (downstream.water.dists$comm.struc), xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,.6))
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
     (downstream.sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
     (downstream.water.dists$comm.struc), xlab="", 
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
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.sed.dists$geo, 
     (downstream.sed.dists$comm.struc), xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000))
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
grid.raster(readPNG("./figures/DDR_DownstreamSedWater_geo.png"))



################################################################
png(filename = "./figures/DDR_HeadwaterDownstream_bray_dt.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))


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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Geographic Distance", side = 1, line = 3, cex = 1.5)
mtext("Downstream", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
grid.raster(readPNG("./figures/DDR_HeadwaterDownstream_bray_dt.png"))



#-------------------------------------------------------
png(filename = "./figures/DDR_WaterSed_bray_dt.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))


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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(sed.dists$geo, 
     residuals(sed.env.lm)+ coefficients(summary(sed.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000), ylim = yrange.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
grid.raster(readPNG("./figures/DDR_WaterSed_bray_dt.png"))



# Downstream sed vs. water
png(filename = "./figures/DDR_DownstreamSedWater_bray_dt.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.sed.dists$env, 
     residuals(downstream.sed.geo.lm) + coefficients(summary(downstream.sed.env.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6), ylim = yrange.ds.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.sed.dists$geo, 
     residuals(downstream.sed.env.lm)+ coefficients(summary(downstream.sed.geo.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,12000), ylim = yrange.ds.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
grid.raster(readPNG("./figures/DDR_DownstreamSedWater_bray_dt.png"))


###### Dendritic distance decay relationships
################################################################
png(filename = "./figures/DDR_HeadwaterDownstream_bray_dt_dend.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

yrange.headwater <- c(
  min(residuals(headwater.den.lm)+coefficients(summary(headwater.env.lm))[1,1],
      residuals(headwater.env.lm)+coefficients(summary(headwater.den.lm))[1,1]),
  max(residuals(headwater.den.lm)+coefficients(summary(headwater.env.lm))[1,1],
      residuals(headwater.env.lm)+coefficients(summary(headwater.den.lm))[1,1]))
yrange.downstream <- c(
  min(residuals(downstream.env.lm)+coefficients(summary(downstream.den.lm))[1,1],
      residuals(downstream.den.lm)+coefficients(summary(downstream.env.lm))[1,1]),
  max(residuals(downstream.env.lm)+coefficients(summary(downstream.den.lm))[1,1],
      residuals(downstream.den.lm)+coefficients(summary(downstream.env.lm))[1,1]))


par(mar = c(1, 5, 3, 0) + 0.4)
plot(headwater.dists$env, 
     residuals(headwater.den.lm)+ coefficients(summary(headwater.env.lm))[1,1], xlab="", 
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.dists$env, 
     residuals(downstream.den.lm)+ coefficients(summary(downstream.env.lm))[1,1] , xlab="", 
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)
mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)

par(mar = c(1, 1, 3, 4) + 0.4)
plot(headwater.dists$den, 
     residuals(headwater.env.lm)+ coefficients(summary(headwater.den.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,20000), ylim = yrange.headwater)
abline(dt.headwater.den.lm, lty = 1, lwd = 2)

r2 <- round(summary(dt.headwater.den.lm)$r.squared, 3)
my.p <- round(summary(dt.headwater.den.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomleft', legend = rp, bty = 'n')


axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Headwaters", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.dists$den, 
     residuals(downstream.env.lm) + coefficients(summary(downstream.den.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,20000), ylim = yrange.downstream)

abline(dt.downstream.den.lm, lwd=2)
r2 <- round(summary(dt.downstream.den.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.den.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext("Dendritic Distance (m)", side = 1, line = 3, cex = 1.5)
mtext("Downstream", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
grid.raster(readPNG("./figures/DDR_HeadwaterDownstream_bray_dt_dend.png"))



#-------------------------------------------------------
png(filename = "./figures/DDR_WaterSed_bray_dt_dend.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

yrange.water <- c(
  min(residuals(water.den.lm)+coefficients(summary(water.env.lm))[1,1],
      residuals(water.env.lm)+coefficients(summary(water.den.lm))[1,1]),
  max(residuals(water.den.lm)+coefficients(summary(water.env.lm))[1,1],
      residuals(water.env.lm)+coefficients(summary(water.den.lm))[1,1]))
yrange.sed <- c(
  min(residuals(sed.env.lm)+coefficients(summary(sed.den.lm))[1,1],
      residuals(sed.den.lm)+coefficients(summary(sed.env.lm))[1,1]),
  max(residuals(sed.env.lm)+coefficients(summary(sed.den.lm))[1,1],
      residuals(sed.den.lm)+coefficients(summary(sed.env.lm))[1,1]))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(water.dists$env, 
     residuals(water.den.lm)+coefficients(summary(water.env.lm))[1,1], xlab="", 
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(sed.dists$env, 
     residuals(sed.den.lm)+ coefficients(summary(sed.env.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6),
     ylim = yrange.sed)
abline(dt.sed.env.lm, lty = 1, lwd = 2)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
plot(water.dists$den, 
     residuals(water.env.lm) + coefficients(summary(water.den.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,20000), ylim = yrange.water)
abline(dt.water.den.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.water.den.lm)$r.squared, 3)
my.p <- round(summary(dt.water.den.lm)$coefficients[2,4],3)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n')

axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(sed.dists$den, 
     residuals(sed.env.lm)+ coefficients(summary(sed.den.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,20000), ylim = yrange.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(dt.sed.den.lm, lwd = 2)
r2 <- round(summary(dt.sed.den.lm)$r.squared, 3)
my.p <- round(summary(dt.sed.den.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Dendritic Distance (m)", side = 1, line = 3, cex = 1.5)
mtext("Sediment Bacteria", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
grid.raster(readPNG("./figures/DDR_WaterSed_bray_dt_dend.png"))



# Downstream sed vs. water
png(filename = "./figures/DDR_DownstreamSedWater_bray_dt_dend.png",
    width = 1600, height = 1600, res = 96*2)
par(mfcol = c(2, 2))

yrange.ds.water <- c(
  min(residuals(downstream.water.den.lm)+coefficients(summary(downstream.water.env.lm))[1,1],
      residuals(downstream.water.env.lm)+coefficients(summary(downstream.water.den.lm))[1,1]),
  max(residuals(downstream.water.den.lm)+coefficients(summary(downstream.water.env.lm))[1,1],
      residuals(downstream.water.env.lm)+coefficients(summary(downstream.water.den.lm))[1,1]))
yrange.ds.sed <- c(
  min(residuals(downstream.sed.env.lm)+coefficients(summary(downstream.sed.den.lm))[1,1],
      residuals(downstream.sed.den.lm)+coefficients(summary(downstream.sed.env.lm))[1,1]),
  max(residuals(downstream.sed.env.lm)+coefficients(summary(downstream.sed.den.lm))[1,1],
      residuals(downstream.sed.den.lm)+coefficients(summary(downstream.sed.env.lm))[1,1]))

par(mar = c(1, 5, 3, 0) + 0.4)
plot(downstream.water.dists$env, 
     residuals(downstream.water.den.lm)+ coefficients(summary(downstream.water.env.lm))[1,1], xlab="", 
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
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd = 2)
mtext("Community Similarity", side = 2, line = 3, cex = 1.2)

par(mar = c(4, 5, 1, 0) + 0.4)
plot(downstream.sed.dists$env, 
     residuals(downstream.sed.den.lm) + coefficients(summary(downstream.sed.env.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,.6), ylim = yrange.ds.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, 
     labels=c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.0"), 
     at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
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
plot(downstream.water.dists$den, 
     residuals(downstream.water.env.lm)+ coefficients(summary(downstream.water.den.lm))[1,1], xlab="", 
     ylab = "", xaxt="n", yaxt="n", xlim = c(0,20000), ylim = yrange.ds.water)
abline(dt.downstream.water.den.lm, lty = 1, lwd = 2)
r2 <- round(summary(dt.downstream.water.den.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.water.den.lm)$coefficients[2,4],3)
if(my.p < 0.001) my.p <- "< 0.001"
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')
axis(side=1, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
mtext("Bacterioplankton", side = 4, line = 1.5, cex = 1.2)
box(lwd = 2)

par(mar = c(4, 1, 1, 4) + 0.4)
plot(downstream.sed.dists$den, 
     residuals(downstream.sed.env.lm)+ coefficients(summary(downstream.sed.den.lm))[1,1], xlab="", 
     ylab = "", xaxt = "n", yaxt = "n", xlim = c(0,20000), ylim = yrange.ds.sed)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, at=(seq(0.1:1, by = 0.1)), lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
abline(dt.downstream.sed.den.lm, lwd = 2)
r2 <- round(summary(dt.downstream.sed.den.lm)$r.squared, 3)
my.p <- round(summary(dt.downstream.sed.den.lm)$coefficients[2,4],3)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topright', legend = rp, bty = 'n')

mtext("Dendritic Distance (m)", side = 1, line = 3, cex = 1.5)
mtext("Sediment", side = 4, line = 1.5, cex = 1.2)

dev.off()
graphics.off()
grid.raster(readPNG("./figures/DDR_DownstreamSedWater_bray_dt_dend.png"))
