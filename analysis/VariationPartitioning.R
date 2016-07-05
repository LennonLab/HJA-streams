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

varpart(hja.db, ~ ., hja.pcnm, data=hja.env)

hja.env <- as.matrix(hja.env)
hja.pcnm <- as.matrix(hja.pcnm)

hja.env.var <- vegan::capscale(hja.db ~ hja.env)
hja.spa.var <- vegan::capscale(hja.db ~ hja.pcnm)
hja.env_spa.var <- vegan::capscale(hja.db ~ hja.env + Condition(hja.pcnm))
hja.spa_env.var <- vegan::capscale(hja.db ~ hja.pcnm + Condition(hja.env))
permutest(hja.env.var, permutations = 999)
permutest(hja.spa.var, permutations = 999)
permutest(hja.env_spa.var, permutations = 999)
permutest(hja.spa_env.var, permutations = 999)

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

varpart(headwaters.db, ~ ., headwaters.pcnm, data=headwaters.env)
varpart(downstream.db, ~ ., downstream.pcnm, data=downstream.env)

headwaters.env <- as.matrix(headwaters.env)
headwaters.pcnm <- as.matrix(headwaters.pcnm)
headwaters.env.var <- vegan::capscale(headwaters.db ~ headwaters.env)
headwaters.spa.var <- vegan::capscale(headwaters.db ~ headwaters.pcnm)
headwaters.env_spa.var <- vegan::capscale(headwaters.db ~ headwaters.env + Condition(headwaters.pcnm))
headwaters.spa_env.var <- vegan::capscale(headwaters.db ~ headwaters.pcnm + Condition(headwaters.env))
permutest(headwaters.env.var, permutations = 999)
permutest(headwaters.spa.var, permutations = 999) # NS
permutest(headwaters.env_spa.var, permutations = 999) # NS
permutest(headwaters.spa_env.var, permutations = 999) # NS

downstream.env <- as.matrix(downstream.env)
downstream.pcnm <- as.matrix(downstream.pcnm)
downstream.env.var <- vegan::capscale(downstream.db ~ downstream.env)
downstream.spa.var <- vegan::capscale(downstream.db ~ downstream.pcnm)
downstream.env_spa.var <- vegan::capscale(downstream.db ~ downstream.env + Condition(downstream.pcnm))
downstream.spa_env.var <- vegan::capscale(downstream.db ~ downstream.pcnm + Condition(downstream.env))
permutest(downstream.env.var, permutations = 999)
permutest(downstream.spa.var, permutations = 999)
permutest(downstream.env_spa.var, permutations = 999)
permutest(downstream.spa_env.var, permutations = 999)

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

varpart(water.db, ~ ., water.pcnm, data=water.env[,2:6])
varpart(sed.db, ~ ., sed.pcnm, data=sed.env[,2:6])

water.env <- as.matrix(water.env)
water.pcnm <- as.matrix(water.pcnm)
water.env.var <- vegan::capscale(water.db ~ water.env)
water.spa.var <- vegan::capscale(water.db ~ water.pcnm)
water.env_spa.var <- vegan::capscale(water.db ~ water.env + Condition(water.pcnm))
water.spa_env.var <- vegan::capscale(water.db ~ water.pcnm + Condition(water.env))
permutest(water.env.var, permutations = 999)
permutest(water.spa.var, permutations = 999)
permutest(water.env_spa.var, permutations = 999) # NS
permutest(water.spa_env.var, permutations = 999) # NS

sed.env <- as.matrix(sed.env)
sed.pcnm <- as.matrix(sed.pcnm)
sed.env.var <- vegan::capscale(sed.db ~ sed.env)
sed.spa.var <- vegan::capscale(sed.db ~ sed.pcnm)
sed.env_spa.var <- vegan::capscale(sed.db ~ sed.env + Condition(sed.pcnm))
sed.spa_env.var <- vegan::capscale(sed.db ~ sed.pcnm + Condition(sed.env))
permutest(sed.env.var, permutations = 999)
permutest(sed.spa.var, permutations = 999) # NS
permutest(sed.env_spa.var, permutations = 999) # NS
permutest(sed.spa_env.var, permutations = 999) # NS
