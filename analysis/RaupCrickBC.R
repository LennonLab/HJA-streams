source("analysis/InitialSetup.R")
source("analysis/DistanceCalcs.R")

regional.abunds <- t(as.matrix(colSums(OTUs)))
regional.relabunds <- decostand(regional.abunds, method = "total")
occupancy.probs <- t(as.matrix(colSums(decostand(OTUs, method = "pa")) / nrow(OTUs)))
site.abunds <- rowSums(OTUs)
site.rich <- specnumber(OTUs)
a <- regional.relabunds * occupancy.probs

# Create a null community based on Stegen et al. 2015

spec.vec <- 1:ncol(OTUs)
RCbc.nulls <- array(NA, c(56, 56, 999))

# stochastic community assembly nulls
for(i in 1:999){
  
  null.comm <- OTUs * 0
  # for first simulation:
  for(row.i in 1:nrow(null.comm)){
    print(paste("run :", i, " -> ", row.i, " : ", site.abunds[row.i], " inds"))
    while(rowSums(null.comm)[row.i] < site.abunds[row.i]){
      
      
      # choose a species based on its occupancy
      local.specs <- sample(x = spec.vec, size = site.rich[row.i], 
                            prob = as.vector(occupancy.probs), replace = FALSE)
      
      local.probs <- decostand(t(as.matrix(regional.abunds[,local.specs])), method = "total")
      
      local.inds <- sample(x = local.specs, size = site.abunds[row.i],
                           prob = as.vector(local.probs), replace = TRUE)
      
      local.abunds <- rle(sort(local.inds))
      
      # add an individual to the local community
      null.comm[row.i, local.abunds$values] <- local.abunds$lengths
    }
  }
  null.bc <- as.matrix(vegdist(decostand(null.comm, method = "total"), method = "bray"))
  RCbc.nulls[,,i] <- null.bc
}
saveRDS(RCbc.nulls, file = "data/null_models/RCbc.null.rda")

obs.bc <- as.matrix(vegdist(OTUsREL, method = "bray"))
site.compares <- expand.grid(site1 = 1:56, site2 = 1:56)
RC.bray <- matrix(NA, nrow = 56, ncol = 56)
for(row.i in 1:nrow(site.compares)){
  site1 <- site.compares[row.i,1]
  site2 <- site.compares[row.i,2]
  pairwise.null <- RCbc.nulls[site1,site2,]
  pairwise.bray <- obs.bc[site1,site2]
  num.greater <- sum(pairwise.null > pairwise.bray)
  num.ties <- sum(pairwise.null == pairwise.bray)
  val <- (((1 * num.greater) + (0.5 * num.ties))/1000 - 0.5) * 2
  RC.bray[site1, site2] <- val
}
RC.bray.dist <- as.dist(RC.bray)

