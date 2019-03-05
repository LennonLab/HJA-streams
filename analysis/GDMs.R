source("analysis/InitialSetup.R")
library(gdm)
geo.dists <- SoDA::geoXY(env$latitude, env$longitude)

hja.biol.dat <- cbind.data.frame(sample = rownames(design), 
                                 as.matrix(1/sqrt(2)*vegdist(OTUsREL.hel, method = 'euclidean', binary = F)))
#hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
hja.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.subs[,-1], order = design[,5])
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
                                   as.matrix(1/sqrt(2)*vegdist(OTUsREL.hel[which(design$habitat == "water"),], method = 'euclid', binary = F)))
water.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.subs[,-1])[which(design$habitat == "water"),]
water.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.subs[,-1], order = design[,5])[which(design$habitat == "water"),]
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
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))
mtext(paste("Deviance Explained: ",round(water.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)

water.splines <- isplineExtract(water.gdm)
matplot(x = water.splines$x[,-1], y = water.splines$y[,-1], type = 'l', xlim = c(-1,1))

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
sed.biol.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'sediment')], as.matrix(mntd.hja)[which(design$habitat == 'sediment'),which(design$habitat == 'sediment')])
#sed.biol.dat <- cbind.data.frame(sample = rownames(design), as.matrix(normalize.matrix(bNTI.dist)))
#sed.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
sed.env.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'sediment')], geo.dists, env.mat[which(design$habitat == 'sediment'),], order = design[which(design$habitat == 'sediment'),5])
sed.sitepair <- formatsitepair(bioData = sed.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = sed.env.dat, 
                               siteColumn = "sample",
                               YColumn = "Y", XColumn = "X")
head(sed.sitepair)
sed.gdm <- gdm(sed.sitepair, geo = TRUE, splines = rep(5,10))
summary(sed.gdm)
coef(sed.gdm)
length(sed.gdm$predictors)
plot(sed.gdm, plot.layout = c(3,4))
sed.splines <- isplineExtract(sed.gdm)
sed.pred <- predict(sed.gdm, sed.sitepair)
graphics.off()
plot(sed.pred, sed.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
text("1:1", x = (max(sed.pred) - .001), y = (max(sed.sitepair$distance) - 0.02))
mtext(paste("Deviance Explained: ",round(sed.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)

# mntd  - water
geo.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == 'water'),]
water.biol.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'water')], as.matrix(mntd.hja)[which(design$habitat == 'water'),which(design$habitat == 'water')])
#water.biol.dat <- cbind.data.frame(sample = rownames(design), as.matrix(normalize.matrix(bNTI.dist)))
#water.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat)
water.env.dat <- cbind.data.frame(sample = rownames(design)[which(design$habitat == 'water')], geo.dists, env.mat[which(design$habitat == 'water'),], order = design[which(design$habitat == 'water'),5])
#water.env.dat <- water.env.dat[,c("sample", "X","Y", "temperature", "TN", "order")]
water.sitepair <- formatsitepair(bioData = water.biol.dat, 
                               bioFormat = 3, abundance = T, 
                               predData = water.env.dat, 
                               siteColumn = "sample",
                               YColumn = "Y", XColumn = "X")
head(water.sitepair)
water.gdm <- gdm(water.sitepair, geo = TRUE, splines = rep(5, 10))
summary(water.gdm)
coef(water.gdm)
length(water.gdm$predictors)
plot(water.gdm, plot.layout = c(3,4))
water.splines <- isplineExtract(water.gdm)
water.pred <- predict(water.gdm, water.sitepair)
graphics.off()
plot(water.pred, water.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
text("1:1", x = (max(water.pred) - .001), y = (max(water.sitepair$distance) - 0.02))
mtext(paste("Deviance Explained: ",round(water.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)




# comparison






# turnover
geo.dists <- geoXY(env$latitude, env$longitude)
hja.biol.dat <- cbind.data.frame(sample = rownames(design), 
                                 as.matrix(beta.pair(OTUs.PA, index.family = 'sorensen')$beta.sim))
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
lines(c(-1,2), c(-1,2))
#text("1:1", x = .98, y = .93)
mtext(paste("Deviance Explained: ",round(hja.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)


geo.dists <- geoXY(env$latitude, env$longitude)
water.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == "water"),]
water.biol.dat <- cbind.data.frame(sample = rownames(design[which(design$habitat == "water"),]), 
                                   as.matrix(beta.pair(OTUs.PA[which(design$habitat == "water"),], index.family = "sorensen")$beta.sim))
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
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
mtext(paste("Deviance Explained: ",round(water.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)

geo.dists <- geoXY(env$latitude, env$longitude)
sediment.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == "sediment"),]
sediment.biol.dat <- cbind.data.frame(sample = rownames(design[which(design$habitat == "sediment"),]), 
                                   as.matrix(beta.pair(OTUs.PA[which(design$habitat == "sediment"),], index.family = "sorensen")$beta.sim))
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
plot(sediment.gdm, plot.layout = c(3,4))
sediment.pred <- predict(sediment.gdm, sediment.sitepair)
graphics.off()
plot(sediment.pred, sediment.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
mtext(paste("Deviance Explained: ",round(sediment.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)

# downstream
geo.dists <- geoXY(env$latitude, env$longitude)
water.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == "water" & design$order > 1),]
water.biol.dat <- cbind.data.frame(sample = rownames(design[which(design$habitat == "water" & design$order > 1),]), 
                                   as.matrix(vegdist(OTUsREL[which(design$habitat == "water" & design$order > 1),], method = 'bray')))
water.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat[,-1])[which(design$habitat == "water" & design$order > 1),]
water.sitepair <- formatsitepair(bioData = water.biol.dat, 
                                 bioFormat = 3, abundance = T, 
                                 predData = water.env.dat, 
                                 siteColumn = "sample", 
                                 YColumn = "Y", XColumn = "X")

head(water.sitepair)
water.gdm <- gdm(water.sitepair, geo = TRUE)
summary(water.gdm)
length(water.gdm$predictors)
water.pred <- water.gdm$predicted
water.splines <- isplineExtract(water.gdm)
water.splines.x <- as.data.frame(water.splines$x)
water.splines.y <- as.data.frame(water.splines$y)
water.splines.x.long <- water.splines.x %>% rowid_to_column() %>% gather(key = spline, value = xpos, -rowid)
water.splines.y.long <- water.splines.y %>% rowid_to_column() %>% gather(key = spline, value = ypos, -rowid)
full_join(water.splines.x.long, water.splines.y.long) %>% filter(spline == "Geographic") %>%
  ggplot(aes(x = xpos, y = ypos)) + 
  geom_line(size = 2) +
  ylab("F(Geographic Distance)\n") +
  xlab("\nGeographic Distance (m)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        aspect.ratio = 1) + 
  ggsave(filename = "figures/gdm-water-splines-geo.pdf", bg = "white", width = 6, height = 6)
full_join(water.splines.x.long, water.splines.y.long) %>% 
  filter(spline != "Geographic", spline != "DOC", spline != "elevation", spline != "ph") %>%
  ggplot(aes(x = xpos, y = ypos, col = spline)) + 
  geom_line(size = 2) +
  ylab("F(Env Variable)") +
  xlab("Env Variable (z-score)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), aspect.ratio = 1) +
  ggsave(filename = "figures/gdm-water-splines-env.pdf", bg = "white", width = 6, height = 6)

pdf("figures/gdm-water-splines.pdf", bg = "white", height = 10, width = 10)
plot(water.gdm, plot.layout = c(3,3), cex = 1.2)
dev.off()
water.pred <- predict(water.gdm, water.sitepair)
graphics.off()

pdf("figures/gdm-water.pdf", bg = "white", height = 6, width = 6)
plot(water.pred, water.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
mtext(paste("Deviance Explained: ",round(water.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)
dev.off()

geo.dists <- geoXY(env$latitude, env$longitude)
sediment.dists <- geoXY(env$latitude, env$longitude)[which(design$habitat == "sediment" & design$order > 1),]
sediment.biol.dat <- cbind.data.frame(sample = rownames(design[which(design$habitat == "sediment" & design$order > 1),]), 
                                      as.matrix(vegdist(OTUsREL[which(design$habitat == "sediment" & design$order > 1),],method = 'bray')))
sediment.env.dat <- cbind.data.frame(sample = rownames(design), geo.dists, env.mat[,-1])[which(design$habitat == "sediment" & design$order > 1),]
sediment.sitepair <- formatsitepair(bioData = sediment.biol.dat, 
                                    bioFormat = 3, abundance = T, 
                                    predData = sediment.env.dat, 
                                    siteColumn = "sample", 
                                    YColumn = "Y", XColumn = "X")

head(sediment.sitepair)
sediment.gdm <- gdm(sediment.sitepair, geo = TRUE)
summary(sediment.gdm)
length(sediment.gdm$predictors)
sediment.splines <- isplineExtract(sediment.gdm)
sediment.splines.x <- as.data.frame(sediment.splines$x)
sediment.splines.y <- as.data.frame(sediment.splines$y)
sediment.splines.x.long <- sediment.splines.x %>% rowid_to_column() %>% gather(key = spline, value = xpos, -rowid)
sediment.splines.y.long <- sediment.splines.y %>% rowid_to_column() %>% gather(key = spline, value = ypos, -rowid)
full_join(sediment.splines.x.long, sediment.splines.y.long) %>% filter(spline == "Geographic") %>%
  ggplot(aes(x = xpos, y = ypos)) + 
  geom_line(size = 2) +
  ylab("F(Geographic Distance)\n") +
  xlab("\nGeographic Distance (m)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        aspect.ratio = 1) + 
  ggsave(filename = "figures/gdm-sediment-splines-geo.pdf", bg = "white", width = 6, height = 6)
full_join(sediment.splines.x.long, sediment.splines.y.long) %>% 
  filter(spline != "Geographic") %>%
  ggplot(aes(x = xpos, y = ypos, col = spline)) + 
  geom_line(size = 2) +
  ylab("F(Env Variable)") +
  xlab("Env Variable (z-score)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), aspect.ratio = 1) +
  ggsave(filename = "figures/gdm-sediment-splines-env.pdf", bg = "white", width = 6, height = 4)


pdf("figures/gdm-sed-splines.pdf", bg = "white", height = 10, width = 10)
plot(sediment.gdm, plot.layout = c(3,3))
dev.off()
sediment.pred <- predict(sediment.gdm, sediment.sitepair)
graphics.off()
pdf("figures/gdm-sed.pdf", bg = "white", height = 6, width = 6)
plot(sediment.pred, sediment.sitepair$distance, ylab="Observed dissimilarity",
     xlab="Predicted dissimilarity", pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2), lwd = 2)
mtext(paste("Deviance Explained: ",round(sediment.gdm$explained,3),"%",sep = ""), 
      side = 3, line = 1.2, cex = 1.2)
dev.off()



gdm.w.geo.plt <- full_join(water.splines.x.long, water.splines.y.long) %>% filter(spline == "Geographic") %>%
  ggplot(aes(x = xpos, y = ypos)) + 
  geom_line(size = 2) +
  ylab("F(Geographic Distance)\n") +
  xlab("\nGeographic Distance (m)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        aspect.ratio = 1)
gdm.w.env.plt <- full_join(water.splines.x.long, water.splines.y.long) %>% 
  filter(spline != "Geographic", spline != "DOC", spline != "elevation", spline != "ph") %>%
  ggplot(aes(x = xpos, y = ypos, col = spline)) + 
  geom_line(size = 2) +
  ylab("F(Env Variable)") +
  xlab("Env Variable (z-score)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), aspect.ratio = 1)
gdm.s.geo.plt <- full_join(sediment.splines.x.long, sediment.splines.y.long) %>% filter(spline == "Geographic") %>%
  ggplot(aes(x = xpos, y = ypos)) + 
  geom_line(size = 2) +
  ylab("F(Geographic Distance)\n") +
  xlab("\nGeographic Distance (m)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        aspect.ratio = 1)
gdm.s.env.plt <- full_join(sediment.splines.x.long, sediment.splines.y.long) %>% 
  filter(spline != "Geographic") %>%
  ggplot(aes(x = xpos, y = ypos, col = spline)) + 
  geom_line(size = 2) +
  ylab("F(Env Variable)") +
  xlab("Env Variable (z-score)") +
  theme_cowplot() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), aspect.ratio = 1)
plot_grid(gdm.w.env.plt, gdm.w.geo.plt, gdm.s.env.plt, gdm.s.geo.plt, 
          labels = c("A", "B", "C", "D")) + 
  ggsave("figures/gdm-grid.pdf", bg = "white", width = 10, height = 10)
