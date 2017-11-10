library(gdm)
geo.dists <- geoXY(env$latitude, env$longitude)
hja.biol.dat <- cbind.data.frame(sample = rownames(design), 
                                 as.matrix(vegdist(OTUsREL, method = 'bray', binary = F)))
#hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat, order = design[,5])
hja.sitepair <- formatsitepair(bioData = hja.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = hja.env.dat, 
                               siteColumn = "sample", 
                               YColumn = "Y", XColumn = "X")

head(hja.sitepair)
hja.gdm <- gdm(hja.sitepair, geo = TRUE)
summary(hja.gdm)
coef(hja.gdm)
length(hja.gdm$predictors)
plot(hja.gdm, plot.layout = c(3,4))
hja.pred <- predict(hja.gdm, hja.sitepair)
graphics.off()
plot(hja.pred, hja.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", xlim=c(0.4,1), ylim=c(0.3,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
#text("1:1", x = .98, y = .93)
mtext(paste("Deviance Explained: ",round(hja.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)



water.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == "water"),]
water.biol.dat <- cbind.data.frame(sample = rownames(design[which(design$habitat == "water"),]), 
                                   as.matrix(vegdist(OTUsREL[which(design$habitat == "water"),], method = 'bray', binary = F)))
water.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat[,-1])[which(design$habitat == "water"),]
water.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat[,-1], order = design[,5])[which(design$habitat == "water"),]
water.sitepair <- formatsitepair(bioData = water.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = water.env.dat, 
                               siteColumn = "sample", 
                               YColumn = "Y", XColumn = "X")

head(water.sitepair)
water.gdm <- gdm(water.sitepair, geo = TRUE)
summary(water.gdm)
length(water.gdm$predictors)
plot(water.gdm, plot.layout = c(3,4))
water.pred <- predict(water.gdm, water.sitepair)
graphics.off()
plot(water.pred, water.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", xlim=c(0.4,1), ylim=c(0.3,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
mtext(paste("Deviance Explained: ",round(water.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)

sediment.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == "sediment"),]
sediment.biol.dat <- cbind.data.frame(sample = rownames(design[which(design$habitat == "sediment"),]), 
                                   as.matrix(vegdist(OTUsREL[which(design$habitat == "sediment"),], method = 'bray', binary = F)))
sediment.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat[,-1])[which(design$habitat == "sediment"),]
sediment.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat[,-1], order = design[,5])[which(design$habitat == "sediment"),]
sediment.sitepair <- formatsitepair(bioData = sediment.biol.dat, 
                                 bioFormat = 3, abundance = T, 
                                 predData = sediment.env.dat, 
                                 siteColumn = "sample", 
                                 YColumn = "Y", XColumn = "X")

head(sediment.sitepair)
sediment.gdm <- gdm(sediment.sitepair, geo = TRUE)
summary(sediment.gdm)
length(sediment.gdm$predictors)
graphics.off()
plot.gdm(sediment.gdm, plot.layout = c(3,4))
sediment.pred <- predict(sediment.gdm, sediment.sitepair)
graphics.off()
plot(sediment.pred, sediment.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", xlim=c(0.4,1), ylim=c(0.3,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
mtext(paste("Deviance Explained: ",round(sediment.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)


# mntd
geo.dists <- geoXY(env$latitude, env$longitude)
hja.biol.dat <- cbind.data.frame(sample = rownames(design), as.matrix(mntd.hja))
#hja.biol.dat <- cbind.data.frame(sample = rownames(design), as.matrix(normalize.matrix(bNTI.dist)))
#hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat, order = design[,5])
hja.sitepair <- formatsitepair(bioData = hja.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = hja.env.dat, 
                               siteColumn = "sample",
                               YColumn = "Y", XColumn = "X")
head(hja.sitepair)
hja.gdm <- gdm(hja.sitepair, geo = TRUE)
summary(hja.gdm)
coef(hja.gdm)
length(hja.gdm$predictors)
plot(hja.gdm, plot.layout = c(3,4))
hja.pred <- predict(hja.gdm, hja.sitepair)
graphics.off()
plot(hja.pred, hja.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
text("1:1", x = (max(hja.pred) - .001), y = (max(hja.sitepair$distance) - 0.02))
mtext(paste("Deviance Explained: ",round(hja.gdm$explained,3),"%",sep = ""), 
     side = 3, line = 1.2, cex = 1.2)

# mntd - sed
geo.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == 'sediment'),]
hja.biol.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'sediment')], as.matrix(mntd.hja)[which(design$habitat == 'sediment'),which(design$habitat == 'sediment')])
#hja.biol.dat <- cbind.data.frame(sample = rownames(design), as.matrix(normalize.matrix(bNTI.dist)))
#hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
hja.env.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'sediment')], geo.dists, env.mat[which(design$habitat == 'sediment'),], order = design[which(design$habitat == 'sediment'),5])
hja.sitepair <- formatsitepair(bioData = hja.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = hja.env.dat, 
                               siteColumn = "sample",
                               YColumn = "Y", XColumn = "X")
head(hja.sitepair)
hja.gdm <- gdm(hja.sitepair, geo = TRUE, splines = rep(5, 10))
summary(hja.gdm)
coef(hja.gdm)
length(hja.gdm$predictors)
plot(hja.gdm, plot.layout = c(3,4))
hja.pred <- predict(hja.gdm, hja.sitepair)
graphics.off()
plot(hja.pred, hja.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
text("1:1", x = (max(hja.pred) - .001), y = (max(hja.sitepair$distance) - 0.02))
mtext(paste("Deviance Explained: ",round(hja.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)

# mntd  - water
geo.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == 'water'),]
hja.biol.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'water')], as.matrix(mntd.hja)[which(design$habitat == 'water'),which(design$habitat == 'water')])
#hja.biol.dat <- cbind.data.frame(sample = rownames(design), as.matrix(normalize.matrix(bNTI.dist)))
#hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
hja.env.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'water')], geo.dists, env.mat[which(design$habitat == 'water'),], order = design[which(design$habitat == 'water'),5])
hja.sitepair <- formatsitepair(bioData = hja.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = hja.env.dat, 
                               siteColumn = "sample",
                               YColumn = "Y", XColumn = "X")
head(hja.sitepair)
hja.gdm <- gdm(hja.sitepair, geo = TRUE, splines = rep(5, 10))
summary(hja.gdm)
coef(hja.gdm)
length(hja.gdm$predictors)
plot(hja.gdm, plot.layout = c(3,4))
hja.pred <- predict(hja.gdm, hja.sitepair)
graphics.off()
plot(hja.pred, hja.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
text("1:1", x = (max(hja.pred) - .001), y = (max(hja.sitepair$distance) - 0.02))
mtext(paste("Deviance Explained: ",round(hja.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)
