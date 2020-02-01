# AEM 
library(spdep)

geoXY

design.sed <- subset(design, habitat == "sediment")
design.water <- subset(design, habitat == "water")
env %>% 
  ggplot(aes(x = longitude, y = latitude)) +
  geom_point(alpha = .2) + 
  geom_text_repel(aes(label = sample), point.padding = .1) +
  facet_wrap(~habitat) +
  theme_minimal() +
  ggsave("temp/sample_orientation.pdf")

ggplot(data = fortify(rgdal::readOGR("imagery/lidar_stream/lidar_stream.shp"))) + 
  geom_polygon(aes(x = long, y = lat, group = group),
               color = "gray", fill = "white", size = 0.2)
  

# # Construct side-by-edge matrix for sediments
# sed.E <- matrix(0, nrow = 24, ncol = 31)
# rownames(sed.E) <- rownames(design.sed)
# sed.E["LC_01_S",c(1:17)] <- 1
# sed.E["LC_02_S",c(2:17)] <- 1
# sed.E["LC_05_S",c(5,6,7,8,9,10,12,13,14,15)] <- 1
# sed.E["LC_08_S",17] <- 1
# sed.E["LC_09_S",13] <- 1
# sed.E["LC_10_S",12] <- 1
# sed.E["LC_11_S",16] <- 1
# sed.E["LC_13_S",c(8,9,10,14,15)] <- 1
# sed.E["LC_16_S",11] <- 1
# sed.E["LC_17_S",c(9,10,14,15)] <- 1
# sed.E["LC_18_S",15] <- 1
# sed.E["LC_19_S",14] <- 1
# sed.E["W1_01_S",c(18:31)] <- 1
# sed.E["W1_02_S",c(19:31)] <- 1
# sed.E["W1_05_S",c(22:31)] <- 1
# sed.E["W1_06_S",c(23:31)] <- 1
# sed.E["W1_08_S",25] <- 1
# sed.E["W1_10_S",c(26:31)] <- 1
# sed.E["W1_12_S",c(27:31)] <- 1
# sed.E["W1_18_S",c(29,30)] <- 1
# sed.E["W1_19_S",30] <- 1
# sed.E["W1_20_S",31] <- 1
# sed.E["W1_21_S",c(20:31)] <- 1
# sed.E["W1_22_S",c(21:31)] <- 1

sed.SBE <- read_csv("data/AEMsed_site-by-edge.csv")
sed.OTUsREL.hel <- OTUsREL.hel[which(design$habitat == "sediment"),]
sed.aem <- aem(binary.mat = as.matrix(sed.SBE[,-1]))
matplot(sed.aem$vectors[,1:9], type = 'l')
plot(sed.aem$values/sum(sed.aem$values))
sed.aem.vec <- as.data.frame(sed.aem$vectors)
sed.db <- vegdist(sed.OTUsREL.hel, method = "euclidean")
sed.aem.full <- rda(sed.OTUsREL.hel ~ ., sed.aem.vec[,1:11])
(sed.r2.aem <- RsquareAdj(sed.aem.full)$adj.r.squared)
sed.aem.fwd <- forward.sel(sed.OTUsREL.hel, sed.aem.vec)
sed.aem.sig.vec <- sed.aem.vec[,c(sort(sed.aem.fwd[,2]))]
sed.aem.sig <- rda(sed.OTUsREL.hel ~ ., sed.aem.sig.vec)
anova(sed.aem.sig)

sed.envs <- sed.SBE %>% rename(sample = "X1") %>% 
  select(sample) %>% 
  left_join(cbind.data.frame(env[,1], env.subs)) %>% 
  select(-sample, -habitat)
varpart(sed.OTUsREL.hel, sed.aem.sig.vec, sed.envs)

mod0 <- dbrda(sed.db ~ 1, data = as.data.frame(sed.aem$vectors))
mod1 <- dbrda(sed.db ~ ., data = as.data.frame(sed.aem$vectors))
envfit(mod1, as.data.frame(sed.aem$vectors))
anova.cca(mod1, by = "axis")
sed.aem.mod <- ordiR2step(mod0, mod1)

sed.aem.2 <- cmdscale(dist(sed.E), eig = T, k = 23)
plot(sed.aem.2$eig/sum(sed.aem.2$eig))
sed.aem.eigvec <- as.data.frame(sed.aem.2$points)
mod1 <- dbrda(sed.db ~ ., data = sed.aem.eigvec)
envfit(mod1, sed.aem.eigvec)
sed.eigs <- subset(sed.aem.eigvec, select = c(V2, V3, V7, V11, V12, V13, V19, V21, V23))
sed.spatial0 <- dbrda(sed.db ~ 1, sed.eigs)
sed.spatial1 <- dbrda(sed.db ~ ., sed.eigs[,c("V2", "V7")])
sed.spatial <- ordiR2step(sed.spatial0, sed.spatial1)
plot(sed.spatial)
envfit(sed.spatial1, sed.eigs)

sed.space.mod <- model.matrix(sed.spatial1)

sed.env <- as.data.frame(env.subs[which(design$habitat == "sediment"),-c(1,2)])
sed.env.mod <- model.matrix(~ .,sed.env)[,-1]
sed.space.mod <- model.matrix(~ ., sed.aem.vec[,1:11])[,-1]
varpart(sed.OTUsREL.hel, sed.env.mod, sed.space.mod)


# # construct water site-by-edge matrix
# water.E <- matrix(0, nrow = 32, ncol = 40)
# rownames(water.E) <- rownames(design.water)
# water.E["LC_02_W",c(1:24)] <- 1
# water.E["LC_03_W",c(2,5:9)] <- 1
# water.E["LC_04_W",c(5:9)] <- 1
# water.E["LC_05_W",c(10:15,17:24)] <- 1
# water.E["LC_07_W",c(8,9)] <- 1
# water.E["LC_08_W",6] <- 1
# water.E["LC_09_W",11] <- 1
# water.E["LC_10_W",13] <- 1
# water.E["LC_11_W",9] <- 1
# water.E["LC_13_W",c(18:24)] <- 1
# water.E["LC_14_W",15] <- 1
# water.E["LC_16_W",16] <- 1
# water.E["LC_18_W",c(23,24)] <- 1
# water.E["LC_20_W",24] <- 1
# water.E["LC_21_W",21] <- 1
# water.E["LC_22_W",19] <- 1
# water.E["W1_02_W",c(25:40)] <- 1
# water.E["W1_03_W",c(26:40)] <- 1
# water.E["W1_21_W",c(27:40)] <- 1
# water.E["W1_22_W",c(28:40)] <- 1
# water.E["W1_04_W",c(29:40)] <- 1
# water.E["W1_05_W",c(30:40)] <- 1
# water.E["W1_06_W",c(31:40)] <- 1
# water.E["W1_08_W",32] <- 1
# water.E["W1_11_W",c(34:40)] <- 1
# water.E["W1_12_W",c(35:40)] <- 1
# water.E["W1_13_W",c(36:40)] <- 1
# water.E["W1_14_W",c(37:40)] <- 1
# water.E["W1_15_W",c(38:40)] <- 1
# water.E["W1_17_W",c(39:40)] <- 1
# water.E["W1_20_W",40] <- 1

waterSBE <- read_csv("data/AEMwater_site-by-edge.csv")
waterOTUsREL.hel <- OTUsREL.hel[which(design$habitat == "water"),]
wateraem <- aem(binary.mat = as.matrix(waterSBE[,-1]))
matplot(wateraem$vectors, type = 'l')
plot(wateraem$values/sum(wateraem$values))
wateraem.vec <- as.data.frame(wateraem$vectors)
waterdb <- vegdist(waterOTUsREL.hel, method = "euclidean")
wateraem.full <- rda(waterOTUsREL.hel ~ ., wateraem.vec[,1:(floor(ncol(wateraem.vec)/2))])
(waterr2.aem <- RsquareAdj(wateraem.full)$adj.r.squared)
wateraem.fwd <- forward.sel(waterOTUsREL.hel, wateraem.vec)
wateraem.sig.vec <- wateraem.vec[,c(sort(wateraem.fwd[,2]))]
wateraem.sig <- rda(waterOTUsREL.hel ~ ., as.data.frame(wateraem.sig.vec))
anova(wateraem.sig)

water.space.mod <- model.matrix(~ ., wateraem.vec[,1:(floor(ncol(wateraem.vec)/2))])[,-1]
water.env <- as.data.frame(env.subs[which(design$habitat == "water"),-c(1,2)])
water.env.mod <- model.matrix(~ .,water.env)[,-1]
varpart(Y = waterOTUsREL.hel, as.data.frame(water.env.mod), as.data.frame(water.space.mod))
