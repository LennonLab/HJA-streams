# setwd("~/GitHub/HJA-streams/analysis/")
# source("./InitialSetup.R")
# source("./Ordination.R")

## Variation Partitioning

#### dbRDA
habitat <- scale((env[,8] == "sediment") * 1)
hja.env <- as.data.frame(cbind(habitat, env.mat[, 2:6]))
colnames(hja.env)[1] <- "habitat"

trunc.dist <- as.matrix(dist(geo.dists))
dist.pcnm <- pcnm(trunc.dist, dist.ret = T)
rs <- rowSums(OTUsREL) / sum(OTUsREL)
hja.pcnm <- pcnm(trunc.dist, w = rs)
hja.pcnm <- scores(hja.pcnm)[,which(hja.pcnm$values>0)]
hja.pcnm <- as.data.frame(hja.pcnm)

capture.output(
  varpart(hja.db, ~ ., hja.pcnm, data=hja.env),
  file = "../tables/hja_varpart.txt")

hja.env.mat <- as.matrix(hja.env)
hja.pcnm.mat <- as.matrix(hja.pcnm)

hja.env.var <- vegan::capscale(hja.db ~ hja.env.mat)
hja.spa.var <- vegan::capscale(hja.db ~ hja.pcnm.mat)
hja.env_spa.var <- vegan::capscale(hja.db ~ hja.env.mat + Condition(hja.pcnm.mat))
hja.spa_env.var <- vegan::capscale(hja.db ~ hja.pcnm.mat + Condition(hja.env.mat))
capture.output(
  permutest(hja.env.var, permutations = 999),
  permutest(hja.spa.var, permutations = 999),
  permutest(hja.env_spa.var, permutations = 999),
  permutest(hja.spa_env.var, permutations = 999),
  file = "../tables/hja_permutation_tests.txt")

# headwaters vs. downstream
headwaters.db <- vegdist(OTUsREL[which(design$order==1),], method = "bray")
downstream.db <- vegdist(OTUsREL[which(design$order!=1),], method = "bray")
headwaters.env <- hja.env[which(design$order==1),] 
downstream.env <- hja.env[which(design$order!=1),]

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

capture.output(
  varpart(headwaters.db, ~ ., headwaters.pcnm, data=headwaters.env),
  file = "../tables/headwaters_varpart.txt")
capture.output(
  varpart(downstream.db, ~ ., downstream.pcnm, data=downstream.env),
  file = "../tables/downstream_varpart.txt")

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
  file = "../tables/headwaters_permutation_tests.txt")

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
  file = "../tables/downstream_permutation_tests.txt")

# water vs. sediments
water.db <- vegdist(OTUsREL[which(design$habitat=="water"),], method = "bray")
sed.db <- vegdist(OTUsREL[which(design$habitat=="sediment"),], method = "bray")
water.env <- hja.env[which(design$habitat=="water"),] 
sed.env <- hja.env[which(design$habitat=="sediment"),]

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

capture.output(
  varpart(water.db, ~ ., water.pcnm, data=water.env[,2:6]),
  file = "../tables/water_varpart.txt")
capture.output(
  varpart(sed.db, ~ ., sed.pcnm, data=sed.env[,2:6]),
  file = "../tables/sed_varpart.txt")

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
  file = "../tables/water_permutation_tests.txt")

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
  file = "../tables/sed_permutation_tests.txt")
