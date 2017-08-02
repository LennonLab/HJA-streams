# AEM 

design.sed <- subset(design, habitat == "sediment")
design.water <- subset(design, habitat == "water")

# Construct side-by-edge matrix for sediments
sed.E <- matrix(0, nrow = 24, ncol = 31)
rownames(sed.E) <- rownames(design.sed)
sed.E["LC_01_S",c(1:17)] <- 1
sed.E["LC_02_S",c(2:17)] <- 1
sed.E["LC_05_S",c(5,6,7,8,9,10,12,13,14,15)] <- 1
sed.E["LC_08_S",17] <- 1
sed.E["LC_09_S",13] <- 1
sed.E["LC_10_S",12] <- 1
sed.E["LC_11_S",16] <- 1
sed.E["LC_13_S",c(8,9,10,14,15)] <- 1
sed.E["LC_16_S",11] <- 1
sed.E["LC_17_S",c(9,10,14,15)] <- 1
sed.E["LC_18_S",15] <- 1
sed.E["LC_19_S",14] <- 1
sed.E["W1_01_S",c(18:31)] <- 1
sed.E["W1_02_S",c(19:31)] <- 1
sed.E["W1_05_S",c(22:31)] <- 1
sed.E["W1_06_S",c(23:31)] <- 1
sed.E["W1_08_S",25] <- 1
sed.E["W1_10_S",c(26:31)] <- 1
sed.E["W1_12_S",c(27:31)] <- 1
sed.E["W1_18_S",c(29,30)] <- 1
sed.E["W1_19_S",30] <- 1
sed.E["W1_20_S",31] <- 1
sed.E["W1_21_S",c(20:31)] <- 1
sed.E["W1_22_S",c(21:31)] <- 1

sed.aem <- aem(binary.mat = sed.E)
plot(sed.aem$values/sum(sed.aem$values))
sed.db <- vegdist(OTUsREL[which(design$habitat == "sediment"),])
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

sed.env <- as.data.frame(env.mat[which(design$habitat == "sediment"),-c(1,2)])
sed.env.mod <- model.matrix(~ .,sed.env)[,-1]
varpart(sed.db, sed.env.mod, sed.space.mod)


# construct water site-by-edge matrix
water.E <- matrix(0, nrow = 32, ncol = 40)
rownames(water.E) <- rownames(design.water)
water.E["LC_02_W",c(1:24)] <- 1
water.E["LC_03_W",c(2,5:9)] <- 1
water.E["LC_04_W",c(5:9)] <- 1
water.E["LC_05_W",c(10:15,17:24)] <- 1
water.E["LC_07_W",c(8,9)] <- 1
water.E["LC_08_W",6] <- 1
water.E["LC_09_W",11] <- 1
water.E["LC_10_W",13] <- 1
water.E["LC_11_W",9] <- 1
water.E["LC_13_W",c(18:24)] <- 1
water.E["LC_14_W",15] <- 1
water.E["LC_16_W",16] <- 1
water.E["LC_18_W",c(23,24)] <- 1
water.E["LC_20_W",24] <- 1
water.E["LC_21_W",21] <- 1
water.E["LC_22_W",19] <- 1
water.E["W1_02_W",c(25:40)] <- 1
water.E["W1_03_W",c(26:40)] <- 1
water.E["W1_21_W",c(27:40)] <- 1
water.E["W1_22_W",c(28:40)] <- 1
water.E["W1_04_W",c(29:40)] <- 1
water.E["W1_05_W",c(30:40)] <- 1
water.E["W1_06_W",c(31:40)] <- 1
water.E["W1_08_W",32] <- 1
water.E["W1_11_W",c(34:40)] <- 1
water.E["W1_12_W",c(35:40)] <- 1
water.E["W1_13_W",c(36:40)] <- 1
water.E["W1_14_W",c(37:40)] <- 1
water.E["W1_15_W",c(38:40)] <- 1
water.E["W1_17_W",c(39:40)] <- 1
water.E["W1_20_W",40] <- 1


water.db <- vegdist(OTUsREL[which(design$habitat == "water"),])
water.aem <- cmdscale(dist(water.E), eig = T, k = 31)
plot(water.aem$eig/sum(water.aem$eig))
water.aem.eigvec <- as.data.frame(water.aem$points)
mod1 <- dbrda(water.db ~ ., data = water.aem.eigvec)
envfit(mod1, water.aem.eigvec)
water.eigs <- subset(water.aem.eigvec, select = c(V1, V2, V6, V10, V27))
water.spatial0 <- dbrda(water.db ~ 1, water.eigs)
water.spatial1 <- dbrda(water.db ~ ., water.eigs)
water.spatial <- ordiR2step(water.spatial0, water.spatial1)
plot(water.spatial1)
water.space.mod <- model.matrix(water.spatial1)

water.env <- as.data.frame(env.mat[which(design$habitat == "water"),-c(1,2)])
water.env.mod <- model.matrix(~ .,water.env)[,-1]
varpart(water.db, water.env.mod, water.space.mod)
?varpart
