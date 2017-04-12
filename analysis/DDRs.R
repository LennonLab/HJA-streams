source("./analysis/InitialSetup.R")
source("./analysis/DistanceCalcs.R")
source("./analysis/Ordination.R")

RC.bray.dist <- readRDS(file = "./data/RCbraydist.csv")

# Water Distance Decay
water.xy <- xy[which(design$habitat == "water"),]
water.den.dist <- as.dist(
  den.dists[which(design$habitat == "water"),
            which(design$habitat == "water")])
water.geo.dist.ls <- na.omit(liste(as.dist(
  as.matrix(dist.mat)[which(design$habitat == "water"), 
           which(design$habitat == "water")]),
  entry = "geo.dist")[,3])
water.struc.dist <- 1 - vegdist(OTUsREL[which(design$habitat == "water"), 
                                            which(design$habitat == "water")],
                                method = dist.met)
water.env.dist <- vegdist(env.2[which(design$habitat == "water"),], 
                          method = "gower", upper = F, diag = F)
rc.water <- as.dist(RC.bray[which(design$habitat == "water"), which(design$habitat == "water")])
water.rc.dist.ls <- liste(rc.water, entry = "rc.bray")[,3]
water.env.dist.ls <- liste(water.env.dist, entry = "env")[,3]
water.struc.dist.ls <- liste(water.struc.dist, entry = "struc")[,3]
water.den.dist.ls <- na.omit(liste(water.den.dist, entry = "den.dist")[,3])
water.phylo.dist.ls <- 1-liste(as.dist(hja.unifrac[which(design$habitat == "water"),
                                         which(design$habitat == "water")]),
                             entry = "unifrac")[,3]
water.dists <- data.frame(water.den.dist.ls, water.geo.dist.ls, 
                          water.struc.dist.ls, water.env.dist.ls,
                          water.phylo.dist.ls, water.rc.dist.ls,
                          water.bnti.dist.ls[,3])
names(water.dists) <- c("den", "geo", "comm.struc", "env", "unifrac", "rc.bray", "bNTI")

water.env.lm <- (lm((water.dists$comm.struc) ~ water.dists$env))
water.geo.lm <- (lm((water.dists$comm.struc) ~ water.dists$geo))
water.den.lm <- (lm((water.dists$comm.struc) ~ water.dists$den))
summary(water.env.lm)
summary(water.geo.lm)
capture.output(summary(water.env.lm), file = "./tables/DDR_water-env.txt")
capture.output(summary(water.geo.lm), file = "./tables/DDR_water-space.txt")

# Sediment Distance Decay
sed.xy <- xy[which(design$habitat == "sediment"),]
sed.den.dist <- as.dist(
  den.dists[which(design$habitat == "sediment"),
            which(design$habitat == "sediment")])
sed.geo.dist.ls <- na.omit(liste(as.dist(
  as.matrix(dist.mat)[which(design$habitat == "sediment"),
                 which(design$habitat == "sediment")]), 
        entry = "geo.dist")[,3])
sed.struc.dist <- 1- vegdist(OTUsREL[which(design$habitat == "sediment"),
                                          which(design$habitat == "sediment")],
                              method = dist.met)
sed.env.dist <- vegdist(env.2[which(design$habitat == "sediment"),], 
                        method = "gower", upper = F, diag = F)
rc.sed <- as.dist(RC.bray[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sed.rc.dist.ls <- liste(rc.sed, entry = "rc.bray")[,3]
sed.env.dist.ls <- liste(sed.env.dist, entry = "env")[,3]
sed.struc.dist.ls <- liste(sed.struc.dist, entry = "struc")[,3]
sed.den.dist.ls <- na.omit(liste(sed.den.dist, entry = "den.dist")[,3])
sed.phylo.dist.ls <- 1-liste(
  as.dist(hja.unifrac[which(design$habitat == "sediment"),
                      which(design$habitat == "sediment")]),
  entry = "unifrac")[,3]
sed.dists <- data.frame(sed.den.dist.ls, sed.geo.dist.ls, 
                          sed.struc.dist.ls, sed.env.dist.ls,
                          sed.phylo.dist.ls, sed.rc.dist.ls,
                        sed.bnti.dist.ls[,3])
names(sed.dists) <- c("den", "geo", "comm.struc", "env", "unifrac", "rc.bray", "bNTI")
plot(sed.dists$comm.struc~sed.dists$env)
sed.env.lm <- (lm((sed.dists$comm.struc) ~ sed.dists$env))
sed.geo.lm <- (lm((sed.dists$comm.struc) ~ sed.dists$geo))
sed.den.lm <- (lm((sed.dists$comm.struc) ~ sed.dists$den))
summary(sed.env.lm)
summary(sed.geo.lm)
summary(sed.den.lm)
capture.output(summary(sed.env.lm), file = "./tables/DDR_sed-env.txt")
capture.output(summary(sed.geo.lm), file = "./tables/DDR_sed-space.txt")

# Headwaters DDRs
headwater.env.dists <- vegdist(env.2[which(design$order < 2),],
                               method = "gower")
headwater.env.dists <- liste(headwater.env.dists, entry = "env")[,3]
headwater.den.dist <- as.dist(den.dists[which(design$order < 2),
                                        which(design$order < 2)])
headwater.den.dists <- na.omit(
  liste(headwater.den.dist, entry = "den.dist")[,3])
headwater.geo.dists <- na.omit(
  liste(as.dist(
    as.matrix(dist.mat)[which(design$order < 2), which(design$order < 2)]),
        entry = "geo.dist")[,3])
headwater.db <- 1 - vegdist(OTUsREL[which(design$order < 2),], 
                        method = dist.met)
headwater.phylo.dist.ls <- 1-liste(
  as.dist(hja.unifrac[which(design$order < 2),
                      which(design$order < 2)]),
  entry = "unifrac")[,3]
headwater.rc <- liste(as.dist(RC.bray[which(design$order < 2), which(design$order < 2)]))[,3]
headwater.dists <- data.frame(
  liste((headwater.db), entry = "comm.struc")[,3],
  headwater.env.dists, headwater.geo.dists, headwater.den.dists,
  headwater.phylo.dist.ls, headwater.rc)
colnames(headwater.dists) <- c("comm.struc", "env", "geo", "den", "unifrac", "rc.bray")
headwater.env.lm <- (lm((headwater.dists$comm.struc) ~ headwater.dists$env))
headwater.geo.lm <- (lm((headwater.dists$comm.struc) ~ headwater.dists$geo))
headwater.den.lm <- (lm((headwater.dists$comm.struc) ~ headwater.dists$den))
headwater.rc.lm <- (lm((headwater.dists$rc.bray) ~ headwater.dists$den))
headwater.phy.env.lm <- (lm((headwater.dists$unifrac) ~ headwater.dists$env))
headwater.phy.geo.lm <- (lm((headwater.dists$unifrac) ~ headwater.dists$geo))
capture.output(summary(headwater.env.lm), file = "./tables/DDR_headwater-env.txt")
capture.output(summary(headwater.geo.lm), file = "./tables/DDR_headwater-space.txt")
capture.output(summary(headwater.den.lm), file = "./tables/DDR_headwater-den.txt")
capture.output(summary(headwater.phy.env.lm), file = "./tables/DDR_headwater-phy-env.txt")
capture.output(summary(headwater.phy.geo.lm), file = "./tables/DDR_headwater-phy-geo.txt")

# Downstream DDRs
downstream.env.dists <- vegdist(env.2[which(design$order >= 2),], 
                                method = "gower")
downstream.env.dists <- liste(downstream.env.dists, entry = "env")[,3]
downstream.geo.dists <- na.omit(
  liste(as.dist(as.matrix(dist.mat)[which(design$order >= 2), which(design$order >= 2)]),
        entry = "geo.dist"))[,3]
downstream.den.dist <- as.dist(
  den.dists[which(design$order > 1),which(design$order > 1)])
downstream.den.dists <- na.omit(
  liste(downstream.den.dist, entry = "den.dist")[,3])

downstream.db <- 1 - vegdist(OTUsREL[which(design$order >= 2),],
                         method = dist.met)
downstream.phylo.dist.ls <- 1-liste(
  as.dist(hja.unifrac[which(design$order >= 2),
                      which(design$order >= 2)]),
  entry = "unifrac")[,3]
downstream.rc <- liste(as.dist(RC.bray[which(design$order >= 2), which(design$order >= 2)]))[,3]

downstream.dists <- data.frame(
  liste((downstream.db), entry = "comm.struc")[,3],
  downstream.env.dists, downstream.geo.dists, downstream.den.dists,
  downstream.phylo.dist.ls, downstream.rc)
colnames(downstream.dists) <- c("comm.struc", "env", "geo", "den", "unifrac", "rc.bray")
downstream.env.lm <- (lm(
  (downstream.dists$comm.struc) ~ downstream.dists$env))
downstream.geo.lm <- (lm(
  (downstream.dists$comm.struc) ~ downstream.dists$geo))
downstream.den.lm <- (lm(
  (downstream.dists$comm.struc) ~ downstream.dists$den))
downstream.phy.env.lm <- lm(
  (downstream.dists$unifrac) ~ downstream.dists$env)
downstream.phy.geo.lm <- lm(
  (downstream.dists$unifrac) ~ downstream.dists$geo)
downstream.rc.lm <- lm(
  downstream.dists$rc.bray ~ downstream.dists$den)

capture.output(summary(downstream.env.lm), 
               file = "./tables/DDR_downstream-env.txt")
capture.output(summary(downstream.geo.lm), 
               file = "./tables/DDR_downstream-space.txt")
capture.output(summary(downstream.den.lm), 
               file = "./tables/DDR_downstream-den.txt")
capture.output(summary(downstream.phy.env.lm), 
               file = "./tables/DDR_downstream-phy-env.txt")
capture.output(summary(downstream.phy.geo.lm), 
               file = "./tables/DDR_downstream-phy-geo.txt")



# Downstream seds
downstream.sed.env.dists <- vegdist(
  env.2[which(design$order > 1 & design$habitat == "sediment"),], 
  method = "gower")
downstream.sed.env.dists <- liste(
  downstream.sed.env.dists, entry = "env")[,3]
downstream.sed.geo.dists <- na.omit(
  liste(as.dist(as.matrix(dist.mat)[which(design$order >1  & design$habitat == "sediment"),
                 which(design$order >1 & design$habitat == "sediment")]),
        entry = "geo.dist"))[,3]
downstream.sed.den.dists <- na.omit(liste(
  as.dist(den.dists[which(design$order >1 & design$habitat == "sediment"),
                    which(design$order >1 & design$habitat == "sediment")]),
  entry = "den.dist"))[,3]
downstream.sed.db <- 1 - vegdist(
  OTUsREL[which(design$order >1 & design$habitat == "sediment"),], 
  method = dist.met)
downstream.sed.phylo.dist.ls <- 1-liste(
  as.dist(hja.unifrac[which(design$order >= 2 & design$habitat == "sediment"),
                      which(design$order >= 2 & design$habitat == "sediment")]),
  entry = "unifrac")[,3]
downstream.sed.rc <- liste(as.dist(RC.bray[which(design$order > 1 & design$habitat == "sediment"),
                                           which(design$order > 1 & design$habitat == "sediment")]))[,3]
downstream.sed.dists <- data.frame(
  liste((downstream.sed.db), entry = "comm.struc")[,3],
  downstream.sed.env.dists, downstream.sed.geo.dists,
  downstream.sed.phylo.dist.ls, downstream.sed.den.dists,
  downstream.sed.rc)
colnames(downstream.sed.dists) <- c("comm.struc", "env", "geo", "unifrac", "den", "rc.bray")
downstream.sed.env.lm <- (lm(
  (downstream.sed.dists$comm.struc) ~ downstream.sed.dists$env))
downstream.sed.geo.lm <- (lm(
  (downstream.sed.dists$comm.struc) ~ downstream.sed.dists$geo))
downstream.sed.phy.env.lm <- (lm(
  (downstream.sed.dists$unifrac) ~ downstream.sed.dists$env))
downstream.sed.phy.geo.lm <- (lm(
  (downstream.sed.dists$unifrac) ~ downstream.sed.dists$geo))
summary(downstream.sed.env.lm)

summary(downstream.sed.phy.env.lm)
summary(downstream.sed.geo.lm)
summary(downstream.sed.phy.geo.lm)


# Downstream water
downstream.water.env.dists <- vegdist(
  env.2[which(design$order > 1 & design$habitat == "water"),], 
  method = "gower")
downstream.water.env.dists <- liste(
  downstream.water.env.dists, entry = "env")[,3]
downstream.water.geo.dists <- na.omit(
  liste(as.dist(as.matrix(dist.mat)[which(design$order >1  & design$habitat == "water"),
                 which(design$order >1 & design$habitat == "water")]),
        entry = "geo.dist"))[,3]
downstream.water.db <- 1 - vegdist(
  OTUsREL[which(design$order >1 & design$habitat == "water"),],
  method = dist.met)
downstream.water.phylo.dist.ls <- 1-liste(
  as.dist(hja.unifrac[which(design$order >= 2 & design$habitat == "water"),
                      which(design$order >= 2 & design$habitat == "water")]),
  entry = "unifrac")[,3]
downstream.water.den.dists <- na.omit(liste(
  as.dist(den.dists[which(design$order >1  & design$habitat == "water"),
                    which(design$order >1 & design$habitat == "water")]),
  entry = "den.dist"))[,3]
downstream.water.rc <- liste(as.dist(RC.bray[which(design$order > 1 & design$habitat == "water"),
                                             which(design$order > 1 & design$habitat == "water")]))[,3]

downstream.water.dists <- data.frame(
  liste((downstream.water.db), entry = "comm.struc")[,3], 
  downstream.water.env.dists, downstream.water.geo.dists,
  downstream.water.phylo.dist.ls, downstream.water.den.dists,
  downstream.water.rc)
colnames(downstream.water.dists) <- c("comm.struc", "env", "geo", "unifrac", "den", "rc.bray")
downstream.water.env.lm <- (lm(
  (downstream.water.dists$comm.struc) ~ downstream.water.dists$env))
downstream.water.geo.lm <- (lm(
  (downstream.water.dists$comm.struc) ~ downstream.water.dists$geo))
downstream.water.phy.env.lm <- lm(
  (downstream.water.dists$unifrac) ~ downstream.water.dists$env)
downstream.water.phy.geo.lm <- lm(
  (downstream.water.dists$unifrac) ~ downstream.water.dists$geo)
summary(downstream.water.env.lm)
summary(downstream.water.phy.env.lm)
summary(downstream.water.geo.lm)
summary(downstream.water.phy.geo.lm)

### Catchment Scale DDRs

hja.env.dists <- vegdist(env.2, method = "gower")
hja.env.dists <- liste(hja.env.dists, entry = "env")[,3]
hja.geo.dists <- na.omit(liste(dist.mat, entry = "geo"))[,3]
hja.den.dists <- na.omit(liste(as.dist(den.dists), entry = "den"))[,3]
hja.com.dists <- 1-liste(hja.db, entry = "comm.struc")[,3]
hja.phy.dists <- 1-liste(hja.unifrac.dist, entry = "unifrac")[,3]
hja.rc.dists <- liste(RC.bray.dist)[,3]
hja.dists <- data.frame(
  hja.com.dists, hja.env.dists, hja.den.dists, hja.geo.dists, hja.phy.dists, hja.rc.dists)
colnames(hja.dists) <- c("comm.struc", "env", "den", "geo", "unifrac", "rc.bray")

hja.env.lm <- (lm((hja.dists$comm.struc) ~ hja.dists$env))
hja.den.lm <- (lm((hja.dists$comm.struc) ~ hja.dists$den))
hja.geo.lm <- (lm((hja.dists$comm.struc) ~ hja.dists$geo))
hja.geo_env.lm <- lm((hja.dists$comm.struc) ~ hja.dists$geo * hja.dists$env)

summary(hja.geo_env.lm)
summary(hja.env.lm)
summary(hja.geo.lm)
summary(hja.den.lm)



###### DETRENDED REGRESSIONS

# Headwaters vs Downstream
dt.headwater.env.lm <- (lm(
  residuals(headwater.den.lm) ~ headwater.dists$env))
dt.downstream.env.lm <- (lm(
  residuals(downstream.den.lm) ~ downstream.dists$env))
dt.headwater.den.lm <- (lm(
  residuals(headwater.env.lm) ~ headwater.dists$den))
dt.downstream.den.lm <- (lm(
  residuals(downstream.env.lm) ~ downstream.dists$den))

summary(dt.headwater.env.lm)
summary(dt.headwater.den.lm)
summary(dt.downstream.env.lm)
summary(dt.downstream.den.lm)

# dt.headwater.env.lm<-(lm(residuals(headwater.geo.lm) + coefficients(summary(headwater.env.lm))[1,1] ~ headwater.dists$env))
# dt.downstream.env.lm<-(lm(residuals(downstream.geo.lm) + coefficients(summary(downstream.env.lm))[1,1] ~ downstream.dists$env))
# dt.headwater.geo.lm<-(lm(residuals(headwater.env.lm) + coefficients(summary(headwater.geo.lm))[1,1] ~ headwater.dists$geo))
# dt.downstream.geo.lm<-(lm(residuals(downstream.env.lm) + coefficients(summary(downstream.geo.lm))[1,1] ~ downstream.dists$geo))

# dt.headwater.env.lm<-(lm(residuals(headwater.den.lm) + coefficients(summary(headwater.env.lm))[1,1] ~ headwater.dists$env))
# dt.downstream.env.lm<-(lm(residuals(downstream.den.lm) + coefficients(summary(downstream.env.lm))[1,1] ~ downstream.dists$env))
# dt.headwater.den.lm<-(lm(residuals(headwater.env.lm) + coefficients(summary(headwater.den.lm))[1,1] ~ headwater.dists$den))
# dt.downstream.den.lm<-(lm(residuals(downstream.env.lm) + coefficients(summary(downstream.den.lm))[1,1] ~ downstream.dists$den))



# Water vs Sediment
dt.water.env.lm <- (lm(
  residuals(water.den.lm) ~ water.dists$env))
dt.sed.env.lm <- (lm(
  residuals(sed.den.lm) ~ sed.dists$env))
dt.water.den.lm <- (lm(
  residuals(water.env.lm) ~ water.dists$den))
dt.sed.den.lm <- (lm(
  residuals(sed.env.lm) ~ sed.dists$den))

summary(dt.downstream.env.lm)
summary(dt.downstream.den.lm)
summary(dt.downstream.env.lm)
summary(dt.downstream.den.lm)


# dt.water.env.lm<-(lm(residuals(water.geo.lm) + coefficients(summary(water.env.lm))[1,1] ~ water.dists$env))
# dt.sed.env.lm<-(lm(residuals(sed.geo.lm) + coefficients(summary(sed.env.lm))[1,1] ~ sed.dists$env))
# dt.water.geo.lm<-(lm(residuals(water.env.lm) + coefficients(summary(water.geo.lm))[1,1] ~ water.dists$geo))
# dt.sed.geo.lm<-(lm(residuals(sed.env.lm) + coefficients(summary(sed.geo.lm))[1,1] ~ sed.dists$geo))

# water.den.lm <- lm((water.dists$comm.struc) ~ water.dists$den)
# sed.den.lm <- lm((sed.dists$comm.struc) ~ sed.dists$den)
# dt.water.env.lm<-(lm(residuals(water.den.lm) + coefficients(summary(water.env.lm))[1,1] ~ water.dists$env))
# dt.sed.env.lm<-(lm(residuals(sed.den.lm) + coefficients(summary(sed.env.lm))[1,1] ~ sed.dists$env))
# dt.water.den.lm<-(lm(residuals(water.env.lm) + coefficients(summary(water.den.lm))[1,1] ~ water.dists$den))
# dt.sed.den.lm<-(lm(residuals(sed.env.lm) + coefficients(summary(sed.den.lm))[1,1] ~ sed.dists$den))


# Downstream Sediment vs Water
downstream.water.den.lm <- lm((downstream.water.dists$comm.struc) ~ downstream.water.dists$den)
downstream.sed.den.lm <- lm((downstream.sed.dists$comm.struc) ~ downstream.sed.dists$den)
dt.downstream.water.env.lm <- (lm(
  residuals(downstream.water.den.lm) ~ downstream.water.dists$env))
dt.downstream.sed.env.lm <- (lm(
  residuals(downstream.sed.den.lm) ~ downstream.sed.dists$env))
dt.downstream.water.den.lm <- (lm(
  residuals(downstream.water.env.lm) ~ downstream.water.dists$den))
dt.downstream.sed.den.lm <- (lm(
  residuals(downstream.sed.env.lm) ~ downstream.sed.dists$den))


# dt.downstream.water.env.lm<-(lm(residuals(downstream.water.geo.lm) + coefficients(summary(downstream.water.env.lm))[1,1] ~ downstream.water.dists$env))
# dt.downstream.sed.env.lm<-(lm(residuals(downstream.sed.geo.lm) + coefficients(summary(downstream.sed.env.lm))[1,1] ~ downstream.sed.dists$env))
# dt.downstream.water.geo.lm<-(lm(residuals(downstream.water.env.lm) + coefficients(summary(downstream.water.geo.lm))[1,1] ~ downstream.water.dists$geo))
# dt.downstream.sed.geo.lm<-(lm(residuals(downstream.sed.env.lm) + coefficients(summary(downstream.sed.geo.lm))[1,1] ~ downstream.sed.dists$geo))

# downstream.water.den.lm <- lm((downstream.water.dists$comm.struc) ~ downstream.water.dists$den)
# downstream.sed.den.lm <- lm((downstream.sed.dists$comm.struc) ~ downstream.sed.dists$den)
# dt.downstream.water.env.lm<-(lm(residuals(downstream.water.den.lm) + coefficients(summary(downstream.water.env.lm))[1,1] ~ downstream.water.dists$env))
# dt.downstream.sed.env.lm<-(lm(residuals(downstream.sed.den.lm) + coefficients(summary(downstream.sed.env.lm))[1,1] ~ downstream.sed.dists$env))
# dt.downstream.water.den.lm<-(lm(residuals(downstream.water.env.lm) + coefficients(summary(downstream.water.den.lm))[1,1] ~ downstream.water.dists$den))
# dt.downstream.sed.den.lm<-(lm(residuals(downstream.sed.env.lm) + coefficients(summary(downstream.sed.den.lm))[1,1] ~ downstream.sed.dists$den))

summary(dt.downstream.sed.env.lm)
summary(dt.downstream.sed.den.lm)
summary(dt.downstream.sed.env.lm)
summary(dt.downstream.sed.den.lm)
plot(residuals(downstream.water.den.lm) ~ downstream.water.dists$env)
plot(residuals(downstream.sed.den.lm) ~ downstream.sed.dists$env)
plot(residuals(downstream.water.env.lm) ~ downstream.water.dists$den)
plot(residuals(downstream.sed.env.lm) ~ downstream.sed.dists$den)


residuals(water.den.lm)
residuals(sed.den.lm)
residuals(water.env.lm)
residuals(sed.env.lm)

water.dists$habitat <- "water"
sed.dists$habitat <- "sediment"
water.dists$comm.struc.detrended.den <- residuals(water.den.lm)
water.dists$comm.struc.detrended.env <- residuals(water.env.lm)
sed.dists$comm.struc.detrended.den <- residuals(sed.den.lm)
sed.dists$comm.struc.detrended.env <- residuals(sed.env.lm)
ddrs <- as.data.frame(rbind(water.dists, sed.dists))
names(ddrs)

ggplot(data = ddrs, aes(x = env, y = comm.struc.detrended.den, color = habitat)) +
  geom_point() +
  stat_smooth(method = "lm")

downstream.dists$orders <- "downstream"
headwater.dists$orders <- "headwaters"
ddrs.updown <- as.data.frame(rbind(downstream.dists, headwater.dists))
names(ddrs.updown)
ggplot(data = ddrs.updown, aes(x = den, y = comm.struc, color = orders)) +
  geom_point()+
  stat_smooth(method = "lm")

downstream.water.dists$habitat <- "water"
downstream.sed.dists$habitat <- "sediment"
ds.ddrs <- as.data.frame(rbind(downstream.water.dists, downstream.sed.dists))
ggplot(data = ds.ddrs, aes(x = env, y = comm.struc, color = habitat)) +
  geom_point()+
  stat_smooth(method = "lm")


