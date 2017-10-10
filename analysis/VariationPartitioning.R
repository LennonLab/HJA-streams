# setwd("~/GitHub/HJA-streams/")
# source("./analysis/InitialSetup.R")
# source("./analysis/Ordination.R")
library(adespatial)
library(spdep)

# Variation Partitioning
hja.coords <- cbind(env$latitude, env$longitude)
hja.dist <- vegdist(OTUsREL.hel, method = "euclidean")
hja.dist <- vegdist(OTUsREL, method = "bray")

# Create env model
hja.dbrda.mod0 <- dbrda(hja.dist ~ 1, as.data.frame(env.mat))
hja.dbrda.mod1 <- dbrda(hja.dist ~ ., as.data.frame(env.mat))
hja.dbrda.env <- ordiR2step(hja.dbrda.mod0, hja.dbrda.mod1, perm.max = 200)
hja.dbrda.env$call
hja.env.mod <- model.matrix(~ sediment + elevation + TP,
                            data = as.data.frame(env.mat))[,-1]
permutest(hja.dbrda.env, permutations = 999)
envfit(hja.dbrda.env, hja.env.mod, perm = 999)




# Create spatial mod
rs <- rowSums(OTUs) / sum(OTUs)
hja.pcnm <- pcnm(den.dists, w = rs, dist.ret = T)
hja.pcnm.mod0 <- dbrda(hja.dist ~ 1, as.data.frame(hja.pcnm$vectors))
hja.pcnm.mod1 <- dbrda(hja.dist ~ ., as.data.frame(hja.pcnm$vectors))
hja.pcnm.mod <- ordiR2step(hja.pcnm.mod0, hja.pcnm.mod1)
hja.pcnm.mod$call
permutest(hja.pcnm.mod, permutations = 999)
hja.space.mod <- model.matrix(~ PCNM1 + PCNM13 + PCNM9, data = as.data.frame(hja.pcnm$vectors))[,-1]
hja.pcnm.mod <- dbrda(hja.dist ~ hja.space.mod)

hja.varpart <- varpart(hja.dist, hja.env.mod, hja.space.mod)
hja.varpart

capture.output(
  hja.varpart,
  file = "./tables/varpart_hja.txt")


hja.env.var <- vegan::capscale(hja.dist ~ hja.env.mod)
hja.spa.var <- vegan::capscale(hja.dist ~ hja.space.mod)
hja.env_spa.var <- vegan::capscale(hja.dist ~ hja.env.mod + Condition(hja.space.mod))
hja.spa_env.var <- vegan::capscale(hja.dist ~ hja.space.mod + Condition(hja.env.mod))
capture.output(
  permutest(hja.env.var, permutations = 999),
  permutest(hja.spa.var, permutations = 999),
  permutest(hja.env_spa.var, permutations = 999),
  permutest(hja.spa_env.var, permutations = 999),
  file = "./tables/varpart_hja_permutation_tests.txt")


# headwaters vs. downstream
headwaters.db <- vegdist(decostand(OTUsREL[which(design$order==1),], method = "hellinger"), method = "euclidean")
downstream.db <- vegdist(decostand(OTUsREL[which(design$order!=1),], method = "hellinger"), method = "euclidean")
headwaters.env <- hja.env.mod[which(design$order==1),] 
downstream.env <- hja.env.mod[which(design$order!=1),]
headwaters.coords <- hja.coords[which(design$order==1),] 
downstream.coords <- hja.coords[which(design$order!=1),]

anova.cca(vegan::capscale(headwaters.db ~ headwaters.coords), permutations = 999)
anova.cca(vegan::capscale(downstream.db ~ downstream.coords), permutations = 999)
# Significant linear trend in downstream, so will condition on coords

headwaters.dist <- as.dist(as.matrix(den.dists)[which(design$order==1),which(design$order==1)])
headwaters.pcnm <- pcnm(headwaters.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$order==1),]) / sum(OTUsREL[which(design$order==1),])
headwaters.pcnm <- pcnm(headwaters.dist, w = rs)
headwaters.pcnm <- scores(headwaters.pcnm)[,which(headwaters.pcnm$values>0)]
headwaters.pcnm <- as.data.frame(headwaters.pcnm)
headwaters.dbrda <- dbrda(headwaters.db ~ ., headwaters.pcnm)
envfit(headwaters.dbrda, headwaters.pcnm)
headwaters.pcnm <- subset(headwaters.pcnm, select = c("PCNM1", "PCNM3", "PCNM4")) 


downstream.dist <- as.dist(as.matrix(den.dists)[which(design$order!=1),which(design$order!=1)])
downstream.pcnm <- pcnm(downstream.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$order!=1),]) / sum(OTUsREL[which(design$order!=1),])
downstream.pcnm <- pcnm(downstream.dist, w = rs)
downstream.pcnm <- scores(downstream.pcnm)[,which(downstream.pcnm$values>0)]
downstream.pcnm <- as.data.frame(downstream.pcnm)
downstream.dbrda <- dbrda(downstream.db ~ ., downstream.pcnm)
envfit(downstream.dbrda, downstream.pcnm)
downstream.pcnm <- subset(downstream.pcnm, select = c("PCNM1", "PCNM7", "PCNM11")) 

headwaters.varpart <- varpart(headwaters.db, headwaters.env, headwaters.pcnm)
downstream.varpart <- varpart(downstream.db, downstream.env, downstream.pcnm)
capture.output(
  headwaters.varpart,
  file = "./tables/varpart_headwaters.txt")
capture.output(
  downstream.varpart,
  file = "./tables/varpart_downstream.txt")

headwaters.env.mat <- as.matrix(headwaters.env)
headwaters.pcnm.mat <- as.matrix(headwaters.pcnm)
headwaters.env.var <- vegan::capscale(headwaters.db ~ headwaters.env.mat)
headwaters.spa.var <- vegan::capscale(headwaters.db ~ headwaters.pcnm.mat)
headwaters.env_spa.var <- vegan::capscale(headwaters.db ~ headwaters.env.mat + Condition(headwaters.pcnm.mat))
headwaters.spa_env.var <- vegan::capscale(headwaters.db ~ headwaters.pcnm.mat + Condition(headwaters.env.mat))
capture.output(
  permutest(headwaters.env.var, permutations = 999),
  permutest(headwaters.spa.var, permutations = 999), # NS
  permutest(headwaters.env_spa.var, permutations = 999), # NS
  permutest(headwaters.spa_env.var, permutations = 999), # NS
  file = "./tables/varpart_headwaters_permutation_tests.txt")

downstream.env.mat <- as.matrix(downstream.env)
downstream.pcnm.mat <- as.matrix(downstream.pcnm)
downstream.env.var <- vegan::capscale(downstream.db ~ downstream.env.mat)
downstream.spa.var <- vegan::capscale(downstream.db ~ downstream.pcnm.mat)
downstream.env_spa.var <- vegan::capscale(downstream.db ~ downstream.env.mat + Condition(downstream.pcnm.mat))
downstream.spa_env.var <- vegan::capscale(downstream.db ~ downstream.pcnm.mat + Condition(downstream.env.mat))
capture.output(
  permutest(downstream.env.var, permutations = 999),
  permutest(downstream.spa.var, permutations = 999),
  permutest(downstream.env_spa.var, permutations = 999),
  permutest(downstream.spa_env.var, permutations = 999),
  file = "./tables/varpart_downstream_permutation_tests.txt")

# water vs. sediments
water.db <- vegdist(decostand(OTUsREL[which(design$habitat=="water"),], method = "hellinger"), method = "euclidean")
sed.db <- vegdist(decostand(OTUsREL[which(design$habitat=="sediment"),], method = "hellinger"), method = "euclidean")
water.env <- env.mat[which(design$habitat=="water"),-c(1,2)] 
sed.env <- env.mat[which(design$habitat=="sediment"),-c(1,2)]

water.dist <- as.dist(as.matrix(den.dists)[which(design$habitat=="water"),which(design$habitat=="water")])
water.pcnm <- pcnm(water.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$habitat=="water"),]) / sum(OTUsREL[which(design$habitat=="water"),])
water.pcnm <- pcnm(water.dist, w = rs)
water.pcnm <- scores(water.pcnm)[,which(water.pcnm$values>0)]
water.pcnm <- as.data.frame(water.pcnm)
water.dbrda <- dbrda(water.db ~ ., water.pcnm)
envfit(water.dbrda, water.pcnm)
water.pcnm <- subset(water.pcnm, select = c("PCNM1", "PCNM3", "PCNM7", "PCNM11")) 

sed.dist <- as.dist(as.matrix(den.dists)[which(design$habitat=="sediment"),which(design$habitat=="sediment")])
sed.pcnm <- pcnm(sed.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$habitat=="sediment"),]) / sum(OTUsREL[which(design$habitat=="sediment"),])
sed.pcnm <- pcnm(sed.dist, w = rs)
sed.pcnm <- scores(sed.pcnm)[,which(sed.pcnm$values>0)]
sed.pcnm <- as.data.frame(sed.pcnm)
sed.dbrda <- dbrda(sed.db ~ ., sed.pcnm)
envfit(sed.dbrda, sed.pcnm)
sed.pcnm <- subset(sed.pcnm, select = c("PCNM4", "PCNM7", "PCNM8", "PCNM13")) 

water.varpart <- varpart(water.db, water.env[,-c(1:2)], water.pcnm)
sed.varpart <- varpart(sed.db, sed.env[,-c(1:2)], sed.pcnm)

capture.output(
  water.varpart,
  file = "./tables/varpart_water.txt")
capture.output(
  sed.varpart,
  file = "./tables/varpart_sed.txt")

water.env.mat <- as.matrix(water.env)
water.pcnm.mat <- as.matrix(water.pcnm)
water.env.var <- vegan::capscale(water.db ~ water.env.mat)
water.spa.var <- vegan::capscale(water.db ~ water.pcnm.mat)
water.env_spa.var <- vegan::capscale(water.db ~ water.env.mat + Condition(water.pcnm.mat))
water.spa_env.var <- vegan::capscale(water.db ~ water.pcnm.mat + Condition(water.env.mat))
capture.output(
  permutest(water.env.var, permutations = 999),
  permutest(water.spa.var, permutations = 999),
  permutest(water.env_spa.var, permutations = 999), # NS
  permutest(water.spa_env.var, permutations = 999), # NS
  file = "./tables/varpart_water_permutation_tests.txt")

sed.env.mat <- as.matrix(sed.env)
sed.pcnm.mat <- as.matrix(sed.pcnm)
sed.env.var <- vegan::capscale(sed.db ~ sed.env.mat)
sed.spa.var <- vegan::capscale(sed.db ~ sed.pcnm.mat)
sed.env_spa.var <- vegan::capscale(sed.db ~ sed.env.mat + Condition(sed.pcnm.mat))
sed.spa_env.var <- vegan::capscale(sed.db ~ sed.pcnm.mat + Condition(sed.env.mat))
capture.output(
  permutest(sed.env.var, permutations = 999),
  permutest(sed.spa.var, permutations = 999), # NS
  permutest(sed.env_spa.var, permutations = 999), # NS
  permutest(sed.spa_env.var, permutations = 999), # NS
  file = "./tables/varpart_sed_permutation_tests.txt")







# water vs. sediments using MEM
water.coords <- xy[which(design$habitat=="water"),]
sed.coords <- xy[which(design$habitat=="sediment"),]

water.nb <- dnearneigh(water.coords, 0, .032)
plot(water.nb, water.coords)

water.lw <- nb2listw(water.nb)
wat.mem <- mem(water.lw)
barplot(attr(wat.mem, "values"))
plot(wat.mem)
plot(wat.mem, SpORcoords = water.coords, nb = water.nb)
mI <- moran.randtest(wat.mem, water.lw, 999)
wat.sig <- which(mI$pvalue < 0.05)
plot(wat.mem[,wat.sig], SpORcoords = water.coords, nb = water.nb)
water.mem.df <- as.data.frame(wat.mem[,wat.sig])

sed.nb <- dnearneigh(sed.coords, 0, .038)
plot(sed.nb, sed.coords)

sed.lw <- nb2listw(sed.nb)
sed.mem <- mem(sed.lw)
barplot(attr(sed.mem, "values"))
plot(sed.mem)
plot(sed.mem, SpORcoords = sed.coords, nb = sed.nb)
mI <- moran.randtest(sed.mem, sed.lw, 999)
sed.sig <- which(mI$pvalue < 0.05)
plot(sed.mem[,sed.sig], SpORcoords = sed.coords, nb = sed.nb)
sed.mem.df <- as.data.frame(sed.mem[,sed.sig])

sed.env.mod0 <- dbrda(sed.db ~ 1, as.data.frame(sed.env))
sed.env.mod <- dbrda(sed.db ~ ., as.data.frame(sed.env))
sed.env.mod <- ordiR2step(sed.env.mod0, sed.env.mod)
envfit(sed.env.mod, sed.env)
sed.env.mod <- model.matrix(sed.env.mod)
sed.spa.mod0 <- dbrda(sed.db ~ 1, as.data.frame(sed.mem.df))
sed.spa.mod <- dbrda(sed.db ~ ., as.data.frame(sed.mem.df))
sed.spa.mod <- ordiR2step(sed.spa.mod0, sed.spa.mod)
sed.spa.mod <- dbrda(sed.db ~ MEM1 + MEM3, as.data.frame(sed.mem.df))
envfit(sed.spa.mod, sed.mem.df)
sed.spa.mod <- model.matrix(sed.spa.mod)

water.env.mod0 <- dbrda(water.db ~ 1, as.data.frame(water.env))
water.env.mod <- dbrda(water.db ~ ., as.data.frame(water.env))
water.env.mod <- ordiR2step(water.env.mod0, water.env.mod)
envfit(water.env.mod, water.env)
water.env.mod <- model.matrix(water.env.mod)
water.spa.mod0 <- dbrda(water.db ~ 1, as.data.frame(water.mem.df))
water.spa.mod <- dbrda(water.db ~ ., as.data.frame(water.mem.df))
water.spa.mod <- ordiR2step(water.spa.mod0, water.spa.mod)
envfit(water.spa.mod, water.mem.df)
water.spa.mod <- model.matrix(water.spa.mod)

sed.varpart <- varpart(sed.db, sed.env.mod, sed.spa.mod)
water.varpart <- varpart(water.db, water.env.mod, water.spa.mod)

capture.output(
  water.varpart,
  file = "./tables/varpart_water.txt")
capture.output(
  sed.varpart,
  file = "./tables/varpart_sed.txt")


water.env.var <- vegan::capscale(water.db ~ water.env.mod)
water.spa.var <- vegan::capscale(water.db ~ water.spa.mod)
water.env_spa.var <- vegan::capscale(water.db ~ water.env.mod + Condition(water.spa.mod))
water.spa_env.var <- vegan::capscale(water.db ~ water.spa.mod + Condition(water.env.mod))
capture.output(
  permutest(water.env.var, permutations = 999),
  permutest(water.spa.var, permutations = 999),
  permutest(water.env_spa.var, permutations = 999), # NS
  permutest(water.spa_env.var, permutations = 999), # NS
  file = "./tables/varpart_water_permutation_tests.txt")

sed.env.var <- vegan::capscale(sed.db ~ sed.env.mod)
sed.spa.var <- vegan::capscale(sed.db ~ sed.spa.mod)
sed.env_spa.var <- vegan::capscale(sed.db ~ sed.env.mod + Condition(sed.spa.mod))
sed.spa_env.var <- vegan::capscale(sed.db ~ sed.spa.mod + Condition(sed.env.mod))
capture.output(
  permutest(sed.env.var, permutations = 999),
  permutest(sed.spa.var, permutations = 999), # NS
  permutest(sed.env_spa.var, permutations = 999), # NS
  permutest(sed.spa_env.var, permutations = 999), # NS
  file = "./tables/varpart_sed_permutation_tests.txt")
