# bNTI analysis

# read sediment null dists 
sed.nulldists <- readRDS("./data/mntds-sed-null-dist.rda")

# read sed bmntd vals
sed.mntds <- readRDS("./data/mntds-sed.rda")


# calculate bNTIs
obs.mntds <- as.matrix(sed.mntds)
site.compares <- expand.grid(site1 = 1:ncol(obs.mntds), site2 = 1:ncol(obs.mntds))
bNTI <- matrix(NA, nrow = nrow(obs.mntds), ncol = ncol(obs.mntds))
for(row.i in 1:nrow(site.compares)){
  site1 <- site.compares[row.i,1]
  site2 <- site.compares[row.i,2]
  pairwise.null <- sed.nulldists[site1,site2,]
  pairwise.mntd <- obs.mntds[site1,site2]
  null.mean <- mean(pairwise.null, na.rm = TRUE)
  null.sd <- sd(pairwise.null, na.rm = TRUE)
  val <- (pairwise.mntd - null.mean) / null.sd
  bNTI[site1, site2] <- val
}

bNTI.sed.dist <- as.dist(bNTI)

hist(bNTI.sed.dist)

# read sediment null dists 
water.nulldists <- readRDS("./data/mntds-water-null-dist.rda")

# read sed bmntd vals
water.mntds <- readRDS("./data/mntds-water.rda")


# calculate bNTIs
obs.mntds <- as.matrix(water.mntds)
site.compares <- expand.grid(site1 = 1:ncol(obs.mntds), site2 = 1:ncol(obs.mntds))
bNTI <- matrix(NA, nrow = nrow(obs.mntds), ncol = ncol(obs.mntds))
for(row.i in 1:nrow(site.compares)){
  site1 <- site.compares[row.i,1]
  site2 <- site.compares[row.i,2]
  pairwise.null <- sed.nulldists[site1,site2,]
  pairwise.mntd <- obs.mntds[site1,site2]
  null.mean <- mean(pairwise.null, na.rm = TRUE)
  null.sd <- sd(pairwise.null, na.rm = TRUE)
  val <- (pairwise.mntd - null.mean) / null.sd
  bNTI[site1, site2] <- val
}

bNTI.water.dist <- as.dist(bNTI)


