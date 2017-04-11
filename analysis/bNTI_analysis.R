library(colorspace)
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
colnames(bNTI) <- rownames(design[which(design$habitat == "sediment"),])
rownames(bNTI) <- rownames(design[which(design$habitat == "sediment"),])

bNTI.sed.dist <- as.dist(bNTI)
hist(bNTI.sed.dist, breaks = 20)
sum(bNTI.sed.dist < 2 & bNTI.sed.dist > -2) / length(bNTI.sed.dist) # undom
sum(bNTI.sed.dist > 2) / length(bNTI.sed.dist) # variable selection
sum(bNTI.sed.dist < -2) / length(bNTI.sed.dist) # homogeneous selection


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
  pairwise.null <- water.nulldists[site1,site2,]
  pairwise.mntd <- obs.mntds[site1,site2]
  null.mean <- mean(pairwise.null, na.rm = TRUE)
  null.sd <- sd(pairwise.null, na.rm = TRUE)
  val <- (pairwise.mntd - null.mean) / null.sd
  bNTI[site1, site2] <- val
}
colnames(bNTI) <- rownames(design[which(design$habitat == "water"),])
rownames(bNTI) <- rownames(design[which(design$habitat == "water"),])

bNTI.water.dist <- as.dist(bNTI)
hist(bNTI.water.dist, breaks = 20)
sum(bNTI.water.dist < 2 & bNTI.water.dist > -2) / length(bNTI.water.dist) # undom
sum(bNTI.water.dist > 2) / length(bNTI.water.dist) # variable selection
sum(bNTI.water.dist < -2) / length(bNTI.water.dist) # homogeneous selection




water.bnti.dist.ls <- liste(bNTI.water.dist)
sed.bnti.dist.ls <- liste(bNTI.sed.dist)
water.rc.dist.ls
sed.rc.dist.ls

water.assembly <- as.data.frame(cbind(water.bnti.dist.ls, water.rc.dist.ls))
names(water.assembly)[c(3,4)] <- c("bNTI", "RC.bray")
sed.assembly <- as.data.frame(cbind(sed.bnti.dist.ls, sed.rc.dist.ls))
names(sed.assembly)[c(3,4)] <- c("bNTI", "RC.bray")

assembly.colors <- rainbow_hcl(4)



sed.mechanism <- vector(length = nrow(sed.assembly))
for(row.i in 1:nrow(sed.assembly)){
  if(sed.assembly[row.i,3] >= -2 && sed.assembly[row.i,3] <= 2){
    if(sed.assembly[row.i,4] < -0.95){
      sed.mechanism[row.i] <- "mass effects"
    }
    if(sed.assembly[row.i,4] > 0.95){
      sed.mechanism[row.i] <- "dispersal limitation"
    }
    if(sed.assembly[row.i,4] <= 0.95 && sed.assembly[row.i,4] >= -0.95){
      sed.mechanism[row.i] <- "undominated"
    }
  }
  if(sed.assembly[row.i,3] < -2){
    sed.mechanism[row.i] <- "homogenizing selection"
  }
  if(sed.assembly[row.i,3] > 2){
    sed.mechanism[row.i] <- "variable selection"
  }
  
}
sed.assembly$habitat <- "sediment"
sed.assembly$mechanism <- sed.mechanism

water.mechanism <- vector(length = nrow(water.assembly))
for(row.i in 1:nrow(water.assembly)){
  if(water.assembly[row.i,3] >= -2 && water.assembly[row.i,3] <= 2){
    if(water.assembly[row.i,4] < -0.95){
      water.mechanism[row.i] <- "mass effects"
    }
    if(water.assembly[row.i,4] > 0.95){
      water.mechanism[row.i] <- "dispersal limitation"
    }
    if(water.assembly[row.i,4] <= 0.95 && water.assembly[row.i,4] >= -0.95){
      water.mechanism[row.i] <- "undominated"
    }
  }
  if(water.assembly[row.i,3] < -2){
    water.mechanism[row.i] <- "homogenizing selection"
  }
  if(water.assembly[row.i,3] > 2){
    water.mechanism[row.i] <- "variable selection"
  }
  
}
water.assembly$mechanism <- water.mechanism
water.assembly$habitat <- "water"

community.assembly <- rbind(water.assembly, sed.assembly)

community.assembly.plot <- ggplot(data = community.assembly, aes(x = bNTI, y = RC.bray, col = mechanism)) +
  facet_grid(~habitat) +
  geom_point(show.legend = T) + 
  labs(x = "bNTI", y = "Raup-Crick_bray-curtis", 
       title = "Community Assembly")
ggsave("figures/comm_assembly.pdf", width = 8, height = 8, units = "in")

