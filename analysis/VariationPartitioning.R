# source("./InitialSetup.R")
# source("./Ordination.R")

### Variation Partitioning

#### Var Part - HJA Catchment
trunc.dist <- as.matrix(dist(geo.dists)) * adj.mat
dist.pcnm <- pcnm(trunc.dist)
ordisurf(geo.dists, scores(dist.pcnm, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(geo.dists, scores(dist.pcnm, choi=2), bubble = 4, main = "PCNM 2")
ordisurf(geo.dists, scores(dist.pcnm, choi=3), bubble = 4, main = "PCNM 3")
ordisplom(dist.pcnm, choices = 1:4)
rs <- rowSums(OTUsREL.log) / sum(OTUsREL.log)
pcnmw <- pcnm(trunc.dist, w = rs)

spatial.var <- vegan::cca(OTUsREL.log ~ scores(pcnmw))
env.var <- vegan::cca(OTUsREL.log ~ env.mat[,c(1:4,6)])
env.spat.var <- vegan::cca(OTUsREL.log ~ scores(pcnmw) +
                             Condition(env.mat[,c(1:4,6)]))
spat.env.var <- vegan::cca(OTUsREL.log ~ env.mat[,c(1:4,6)] +
                             Condition(scores(pcnmw)))



sum(env.var$CCA$eig)
sum(spatial.var$CCA$eig)
sum(env.spat.var$CCA$eig)
sum(spat.env.var$CCA$eig)

hja.ca <- vegan::cca(OTUsREL.log)
sum(hja.ca$CA$eig)

var.part.1 <- sum(env.var$CCA$eig) / sum(hja.ca$CA$eig)
var.part.2 <- sum(spatial.var$CCA$eig) / sum(hja.ca$CA$eig)
var.part.3 <- sum(spat.env.var$CCA$eig) / sum(hja.ca$CA$eig)
var.part.4 <- sum(env.spat.var$CCA$eig) / sum(hja.ca$CA$eig)

print(paste("HJA Pure Env. = ", var.part.1))
print(paste("HJA Pure Space. = ", var.part.2))
print(paste("HJA Spatial + Env = ", (var.part.1 - var.part.3)))
print(paste("HJA Residual Var. = ", (1- sum(var.part.1 + var.part.2))))


#### Var Part - Surface Water
env.mat <- env.mat[,c(1:4,6)]
water.dist.pcnm <- pcnm(dist(xy[which(design$habitat == "water"),]))
water.rs <- rowSums(OTUsREL.log[which(design$habitat == "water"),]) / 
  sum(OTUsREL.log[which(design$habitat == "water"),])
water.pcnmw <- pcnm(dist(xy[which(design$habitat == "water"),]), w = water.rs)

water.spatial.var <- vegan::rda(OTUsREL.log[which(design$habitat == "water"),] ~ scores(water.pcnmw))
water.env.var <- vegan::rda(OTUsREL.log[which(design$habitat == "water"),] ~ env.mat[which(design$habitat == "water"),])
water.env.spat.var <- vegan::rda(OTUsREL.log[which(design$habitat == "water"),] ~
                                   scores(water.pcnmw) + Condition(env.mat[which(design$habitat == "water"),]))
water.spat.env.var <- vegan::rda(OTUsREL.log[which(design$habitat == "water"),] ~ 
                                   env.mat[which(design$habitat == "water"),] + Condition(scores(water.pcnmw)))

sum(water.env.var$CCA$eig)
sum(water.spatial.var$CCA$eig)
sum(water.env.spat.var$CCA$eig)
sum(water.spat.env.var$CCA$eig)

water.ca <- vegan::rda(OTUsREL.log[which(design$habitat == "water"),])
sum(water.ca$CA$eig)

water.var.part.1 <- sum(water.env.var$CCA$eig) / sum(water.ca$CA$eig)
water.var.part.2 <- sum(water.spatial.var$CCA$eig) / sum(water.ca$CA$eig)
water.var.part.3 <- sum(water.spat.env.var$CCA$eig) / sum(water.ca$CA$eig)
water.var.part.4 <- sum(water.env.spat.var$CCA$eig) / sum(water.ca$CA$eig)

water.table <- matrix(c(
  rbind(c("Water", "Pure Env. = ", water.var.part.1)),
  rbind(c("Water", "Pure Space. = ", water.var.part.2)),
  rbind(c("Water", "Spatial + Env = ", (water.var.part.1 - water.var.part.3))),
  rbind(c("Water", "Residual Var. = ", (1- sum(water.var.part.1 + water.var.part.2))))
), nrow = 4, byrow = T)
colnames(water.table) <- c("Habitat", "Variation Partition", "Value")

#### Var Part - Sediments Only
sediment.dist.pcnm <- pcnm(dist(xy[which(design$habitat == "sediment"),]))
sediment.rs <- rowSums(OTUsREL.log[which(design$habitat == "sediment"),]) / 
  sum(OTUsREL.log[which(design$habitat == "sediment"),])
sediment.pcnmw <- pcnm(dist(xy[which(design$habitat == "sediment"),]), w = sediment.rs)

sediment.spatial.var <- vegan::rda(OTUsREL.log[which(design$habitat == "sediment"),] ~ scores(sediment.pcnmw))
sediment.env.var <- vegan::rda(OTUsREL.log[which(design$habitat == "sediment"),] ~ env.mat[which(design$habitat == "sediment"),])
sediment.env.spat.var <- vegan::rda(OTUsREL.log[which(design$habitat == "sediment"),] ~
                                      scores(sediment.pcnmw) + Condition(env.mat[which(design$habitat == "sediment"),]))
sediment.spat.env.var <- vegan::rda(OTUsREL.log[which(design$habitat == "sediment"),] ~ 
                                      env.mat[which(design$habitat == "sediment"),] + Condition(scores(sediment.pcnmw)))

sum(sediment.env.var$CCA$eig)
sum(sediment.spatial.var$CCA$eig)
sum(sediment.env.spat.var$CCA$eig)
sum(sediment.spat.env.var$CCA$eig)

sediment.ca <- vegan::rda(OTUsREL.log[which(design$habitat == "sediment"),])
sum(sediment.ca$CA$eig)

sediment.var.part.1 <- sum(sediment.env.var$CCA$eig) / sum(sediment.ca$CA$eig)
sediment.var.part.2 <- sum(sediment.spatial.var$CCA$eig) / sum(sediment.ca$CA$eig)
sediment.var.part.3 <- sum(sediment.spat.env.var$CCA$eig) / sum(sediment.ca$CA$eig)
sediment.var.part.4 <- sum(sediment.env.spat.var$CCA$eig) / sum(sediment.ca$CA$eig)

print(paste("Sediment Pure Env. = ", sediment.var.part.1))
print(paste("Sediment Pure Space. = ", sediment.var.part.2))
print(paste("Sediment Spatial + Env = ", (sediment.var.part.1 - sediment.var.part.3)))
print(paste("Sediment Residual Var. = ", (1- sum(sediment.var.part.1 + sediment.var.part.2))))

perm.sed.spa <- permutest(sediment.spatial.var, permutations = 999)
perm.sed.env <- permutest(sediment.env.var, permutations = 999)
perm.sed.env_spa <- permutest(sediment.spatial.var, permutations = 999)
perm.sed.spa_env <- permutest(sediment.spatial.var, permutations = 999)

sediment.table <- matrix(c(
  rbind(c("Sediment", "Pure Env. = ", sediment.var.part.1)),
  rbind(c("Sediment", "Pure Space. = ", sediment.var.part.2)),
  rbind(c("Sediment", "Spatial + Env = ", (sediment.var.part.1 - sediment.var.part.3))),
  rbind(c("Sediment", "Residual Var. = ", (1- sum(sediment.var.part.1 + sediment.var.part.2))))
), nrow = 4, byrow = T)
colnames(sediment.table) <- c("Habitat", "Variation Partition", "Value")


habitat <- (env[,8] == "sediment") * 1
env.dummy <- as.matrix(cbind(habitat, env[,c(10, 11, 12, 13, 15)]))

d.spatial.var <- spatial.var
d.env.var <- vegan::cca(OTUsREL.log ~ env.dummy)
d.env.spat.var <- vegan::cca(OTUsREL.log ~ scores(pcnmw) +
                               Condition(env.dummy))
d.spat.env.var <- vegan::cca(OTUsREL.log ~ env.dummy +
                               Condition(scores(pcnmw)))

sum(d.env.var$CCA$eig)
sum(d.spatial.var$CCA$eig)
sum(d.env.spat.var$CCA$eig)
sum(d.spat.env.var$CCA$eig)

d.var.part.1 <- sum(d.env.var$CCA$eig) / sum(hja.ca$CA$eig)
d.var.part.2 <- sum(d.spatial.var$CCA$eig) / sum(hja.ca$CA$eig)
d.var.part.3 <- sum(d.spat.env.var$CCA$eig) / sum(hja.ca$CA$eig)
d.var.part.4 <- sum(d.env.spat.var$CCA$eig) / sum(hja.ca$CA$eig)
