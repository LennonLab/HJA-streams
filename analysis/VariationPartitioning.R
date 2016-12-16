# setwd("~/GitHub/HJA-streams/")
# source("./analysis/InitialSetup.R")
# source("./analysis/Ordination.R")

# Variation Partitioning


#### dbRDA
habitat <- scale((env[,8] == "sediment") * 1)
hja.coords <- scale((env[,4:5]))
hja.env <- as.data.frame(cbind(habitat, env.mat))
colnames(hja.env)[1] <- "habitat"

trunc.dist <- as.matrix(dist(geo.dists))
dist.pcnm <- pcnm(trunc.dist, dist.ret = T)
rs <- rowSums(OTUsREL) / sum(OTUsREL)
hja.pcnm <- pcnm(trunc.dist, w = rs)

# First 4 eigenvalues contain 98% of variation
(hja.pcnm$values[hja.pcnm$values>0]/sum(hja.pcnm$values[hja.pcnm$values>0]))
hja.pcnm <- scores(hja.pcnm)[,which(hja.pcnm$values>0)]
hja.pcnm <- as.data.frame(hja.pcnm[,1:4])

hja.varpart <- varpart(hja.db, hja.env, hja.pcnm, hja.coords)

capture.output(
  hja.varpart,
  file = "./tables/hja_varpart.txt")

hja.env.mat <- as.matrix(hja.env)
hja.pcnm.mat <- as.matrix(hja.pcnm)

anova.cca(vegan::capscale(hja.db ~ hja.coords), permutations = 999)

hja.env.var <- vegan::capscale(hja.db ~ hja.env.mat)
hja.spa.var <- vegan::capscale(hja.db ~ hja.pcnm.mat + Condition(hja.coords))
hja.env_spa.var <- vegan::capscale(hja.db ~ hja.env.mat + Condition(hja.pcnm.mat) + Condition(hja.coords))
hja.spa_env.var <- vegan::capscale(hja.db ~ hja.pcnm.mat + Condition(hja.env.mat) + Condition(hja.coords))
capture.output(
  permutest(hja.env.var, permutations = 999),
  permutest(hja.spa.var, permutations = 999),
  permutest(hja.env_spa.var, permutations = 999),
  permutest(hja.spa_env.var, permutations = 999),
  file = "./tables/hja_permutation_tests.txt")


# headwaters vs. downstream
headwaters.db <- vegdist(OTUsREL[which(design$order==1),], method = "bray")
downstream.db <- vegdist(OTUsREL[which(design$order!=1),], method = "bray")
headwaters.env <- hja.env[which(design$order==1),] 
downstream.env <- hja.env[which(design$order!=1),]
headwaters.coords <- hja.coords[which(design$order==1),] 
downstream.coords <- hja.coords[which(design$order!=1),]

anova.cca(vegan::capscale(headwaters.db ~ headwaters.coords), permutations = 999)
anova.cca(vegan::capscale(downstream.db ~ downstream.coords), permutations = 999)
# Significant linear trend in downstream, so will condition on coords

headwaters.dist <- as.matrix(dist(geo.dists[which(design$order==1),]))
headwaters.pcnm <- pcnm(headwaters.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$order==1),]) / sum(OTUsREL[which(design$order==1),])
headwaters.pcnm <- pcnm(headwaters.dist, w = rs)
headwaters.pcnm <- scores(headwaters.pcnm)[,which(headwaters.pcnm$values>0)]
headwaters.pcnm <- as.data.frame(headwaters.pcnm)

downstream.dist <- as.matrix(dist(geo.dists[which(design$order!=1),]))
downstream.pcnm <- pcnm(downstream.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$order!=1),]) / sum(OTUsREL[which(design$order!=1),])
downstream.pcnm <- pcnm(downstream.dist, w = rs)
downstream.pcnm <- scores(downstream.pcnm)[,which(downstream.pcnm$values>0)]
downstream.pcnm <- as.data.frame(downstream.pcnm)

headwaters.varpart <- varpart(headwaters.db, headwaters.env, headwaters.pcnm)
downstream.varpart <- varpart(downstream.db, downstream.env, downstream.pcnm, downstream.coords)
capture.output(
  headwaters.varpart,
  file = "./tables/headwaters_varpart.txt")
capture.output(
  downstream.varpart,
  file = "./tables/downstream_varpart.txt")

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
  file = "./tables/headwaters_permutation_tests.txt")

downstream.env.mat <- as.matrix(downstream.env)
downstream.pcnm.mat <- as.matrix(downstream.pcnm)
downstream.env.var <- vegan::capscale(downstream.db ~ downstream.env.mat)
downstream.spa.var <- vegan::capscale(downstream.db ~ downstream.pcnm.mat + Condition(downstream.coords))
downstream.env_spa.var <- vegan::capscale(downstream.db ~ downstream.env.mat + Condition(downstream.pcnm.mat) + Condition(downstream.coords))
downstream.spa_env.var <- vegan::capscale(downstream.db ~ downstream.pcnm.mat + Condition(downstream.env.mat) + Condition(downstream.coords))
capture.output(
  permutest(downstream.env.var, permutations = 999),
  permutest(downstream.spa.var, permutations = 999),
  permutest(downstream.env_spa.var, permutations = 999),
  permutest(downstream.spa_env.var, permutations = 999),
  file = "./tables/downstream_permutation_tests.txt")

# water vs. sediments
water.db <- vegdist(OTUsREL[which(design$habitat=="water"),], method = "bray")
sed.db <- vegdist(OTUsREL[which(design$habitat=="sediment"),], method = "bray")
water.env <- hja.env[which(design$habitat=="water"),] 
sed.env <- hja.env[which(design$habitat=="sediment"),]
water.coords <- hja.coords[which(design$habitat=="water"),]
sed.coords <- hja.coords[which(design$habitat=="sediment"),]

water.dist <- as.matrix(dist(geo.dists[which(design$habitat=="water"),]))
water.pcnm <- pcnm(water.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$habitat=="water"),]) / sum(OTUsREL[which(design$habitat=="water"),])
water.pcnm <- pcnm(water.dist, w = rs)
water.pcnm <- scores(water.pcnm)[,which(water.pcnm$values>0)]
water.pcnm <- as.data.frame(water.pcnm)

sed.dist <- as.matrix(dist(geo.dists[which(design$habitat=="sediment"),]))
sed.pcnm <- pcnm(sed.dist, dist.ret = T)
rs <- rowSums(OTUsREL[which(design$habitat=="sediment"),]) / sum(OTUsREL[which(design$habitat=="sediment"),])
sed.pcnm <- pcnm(sed.dist, w = rs)
sed.pcnm <- scores(sed.pcnm)[,which(sed.pcnm$values>0)]
sed.pcnm <- as.data.frame(sed.pcnm)

anova.cca(vegan::capscale(water.db ~ water.coords), permutations = 999)
anova.cca(vegan::capscale(sed.db ~ sed.coords), permutations = 999)

water.varpart <- varpart(water.db, water.env[,2:7], water.pcnm, water.coords)
sed.varpart <- varpart(sed.db, sed.env[,2:7], sed.pcnm, sed.coords)

capture.output(
  water.varpart,
  file = "./tables/water_varpart.txt")
capture.output(
  sed.varpart,
  file = "./tables/sed_varpart.txt")

water.env.mat <- as.matrix(water.env)
water.pcnm.mat <- as.matrix(water.pcnm)
water.env.var <- vegan::capscale(water.db ~ water.env.mat)
water.spa.var <- vegan::capscale(water.db ~ water.pcnm.mat + Condition(water.coords))
water.env_spa.var <- vegan::capscale(water.db ~ water.env.mat + Condition(water.pcnm.mat) + Condition(water.coords))
water.spa_env.var <- vegan::capscale(water.db ~ water.pcnm.mat + Condition(water.env.mat) + Condition(water.coords))
capture.output(
  permutest(water.env.var, permutations = 999),
  permutest(water.spa.var, permutations = 999),
  permutest(water.env_spa.var, permutations = 999), # NS
  permutest(water.spa_env.var, permutations = 999), # NS
  file = "./tables/water_permutation_tests.txt")

sed.env.mat <- as.matrix(sed.env)
sed.pcnm.mat <- as.matrix(sed.pcnm)
sed.env.var <- vegan::capscale(sed.db ~ sed.env.mat)
sed.spa.var <- vegan::capscale(sed.db ~ sed.pcnm.mat + Condition(sed.coords))
sed.env_spa.var <- vegan::capscale(sed.db ~ sed.env.mat + Condition(sed.pcnm.mat) + Condition(sed.coords))
sed.spa_env.var <- vegan::capscale(sed.db ~ sed.pcnm.mat + Condition(sed.env.mat) + Condition(sed.coords))
capture.output(
  permutest(sed.env.var, permutations = 999),
  permutest(sed.spa.var, permutations = 999), # NS
  permutest(sed.env_spa.var, permutations = 999), # NS
  permutest(sed.spa_env.var, permutations = 999), # NS
  file = "./tables/sed_permutation_tests.txt")

