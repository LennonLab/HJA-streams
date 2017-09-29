#source("./analysis/InitialSetup.R")


ddr.summary <- data.frame()

response.matrix <- "otus"
update.plots <- FALSE

# Water DDR
water.dists <- list()  
water.dists$otus <- vegdist(OTUsREL[which(design$habitat == "water"),], method = "bray")
water.dists$geo <- earth.dist(xy[which(design$habitat == "water"),])*1000
water.dists$den <- as.dist(as.matrix(den.dists)[which(design$habitat == "water"), which(design$habitat == "water")])
water.dists$env <- dist(env.mat[which(design$habitat == "water"),])
water.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$habitat == "water"), which(design$habitat == "water")])
water.dists$phylo <- as.dist(hja.unifrac[which(design$habitat == "water"), which(design$habitat == "water")])
water.beta.part <- beta.pair(OTUs.PA[which(design$habitat == "water"),], index.family = "sorensen")
water.dists$nest <- water.beta.part$beta.sne
water.dists$turn <- water.beta.part$beta.sim
water.dists$sor <- water.beta.part$beta.sor
water.dists$bNTI <- normalize.matrix(bNTI.water.dist)
water.lms <- DDR(dists = water.dists, comm = response.matrix)
ddr.summary <- fill.table(water.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-water-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(water.lms)
  dev.off()
}

# Sediment Distance Decay
sediment.dists <- list()  
sediment.dists$otus <- vegdist(OTUsREL[which(design$habitat == "sediment"),], method = "bray")
sediment.dists$geo <- earth.dist(xy[which(design$habitat == "sediment"),])*1000
sediment.dists$den <- as.dist(as.matrix(den.dists)[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sediment.dists$env <- dist(env.mat[which(design$habitat == "sediment"),])
sediment.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sediment.dists$phylo <- as.dist(hja.unifrac[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sediment.beta.part <- beta.pair(OTUs.PA[which(design$habitat == "sediment"),], index.family = "sorensen")
sediment.dists$nest <- sediment.beta.part$beta.sne
sediment.dists$turn <- sediment.beta.part$beta.sim
sediment.dists$sor <- sediment.beta.part$beta.sor
sediment.dists$bNTI <- normalize.matrix(bNTI.sed.dist)
sediment.lms <- DDR(dists = sediment.dists, comm = response.matrix)
ddr.summary <- fill.table(sediment.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-sediment-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(sediment.lms)
  dev.off()
}

# Headwaters DDRs
headwater.dists <- list()  
headwater.dists$otus <- vegdist(OTUsREL[which(design$order == 1),], method = "bray")
headwater.dists$geo <- earth.dist(xy[which(design$order == 1),])*1000
headwater.dists$den <- as.dist(as.matrix(den.dists)[which(design$order == 1), which(design$order == 1)])
headwater.dists$env <- dist(env.mat[which(design$order == 1),])
headwater.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order == 1), which(design$order == 1)])
headwater.dists$phylo <- as.dist(hja.unifrac[which(design$order == 1), which(design$order == 1)])
headwater.beta.part <- beta.pair(OTUs.PA[which(design$order == 1),], index.family = "sorensen")
headwater.dists$nest <- headwater.beta.part$beta.sne
headwater.dists$turn <- headwater.beta.part$beta.sim
headwater.dists$sor <- headwater.beta.part$beta.sor
headwater.lms <- DDR(dists = headwater.dists, comm = response.matrix)
ddr.summary <- fill.table(headwater.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-headwater-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(headwater.lms)
  dev.off()
}

# headwater seds
headwater.seds.dists <- list()  
headwater.seds.dists$otus <- vegdist(OTUsREL[which(design$order == 1 & design$habitat == "sediment"),], method = "bray")
headwater.seds.dists$geo <- earth.dist(xy[which(design$order == 1 & design$habitat == "sediment"),])*1000
headwater.seds.dists$den <- as.dist(as.matrix(den.dists)[which(design$order == 1 & design$habitat == "sediment"),which(design$order == 1 & design$habitat == "sediment")])
headwater.seds.dists$env <- dist(env.mat[which(design$order == 1 & design$habitat == "sediment"),])
headwater.seds.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order == 1 & design$habitat == "sediment"),which(design$order == 1 & design$habitat == "sediment")])
headwater.seds.dists$phylo <- as.dist(hja.unifrac[which(design$order == 1 & design$habitat == "sediment"),which(design$order == 1 & design$habitat == "sediment")])
headwater.seds.beta.part <- beta.pair(OTUs.PA[which(design$order == 1 & design$habitat == "sediment"),], index.family = "sorensen")
headwater.seds.dists$nest <- headwater.seds.beta.part$beta.sne
headwater.seds.dists$turn <- headwater.seds.beta.part$beta.sim
headwater.seds.dists$sor <- headwater.seds.beta.part$beta.sor
headwater.seds.lms <- DDR(dists = headwater.seds.dists, comm = response.matrix)
ddr.summary <- fill.table(headwater.seds.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-headwater-seds-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(headwater.seds.lms)
  dev.off()
}

# headwater water
headwater.water.dists <- list()  
headwater.water.dists$otus <- vegdist(OTUsREL[which(design$order == 1 & design$habitat == "water"),], method = "bray")
headwater.water.dists$geo <- earth.dist(xy[which(design$order == 1 & design$habitat == "water"),])*1000
headwater.water.dists$den <- as.dist(as.matrix(den.dists)[which(design$order == 1 & design$habitat == "water"),which(design$order == 1 & design$habitat == "water")])
headwater.water.dists$env <- dist(env.mat[which(design$order == 1 & design$habitat == "water"),])
headwater.water.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order == 1 & design$habitat == "water"),which(design$order == 1 & design$habitat == "water")])
headwater.water.dists$phylo <- as.dist(hja.unifrac[which(design$order == 1 & design$habitat == "water"),which(design$order == 1 & design$habitat == "water")])
headwater.water.beta.part <- beta.pair(OTUs.PA[which(design$order == 1 & design$habitat == "water"),], index.family = "sorensen")
headwater.water.dists$nest <- headwater.water.beta.part$beta.sne
headwater.water.dists$turn <- headwater.water.beta.part$beta.sim
headwater.water.dists$sor <- headwater.water.beta.part$beta.sor
headwater.water.lms <- DDR(dists = headwater.water.dists, comm = response.matrix)
ddr.summary <- fill.table(headwater.water.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-headwater-water-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(headwater.water.lms)
  dev.off()
}

# Downstream DDRs
downstream.dists <- list()  
downstream.dists$otus <- vegdist(OTUsREL[which(design$order != 1),], method = "bray")
downstream.dists$geo <- earth.dist(xy[which(design$order != 1),])*1000
downstream.dists$den <- as.dist(as.matrix(den.dists)[which(design$order != 1), which(design$order != 1)])
downstream.dists$env <- dist(env.mat[which(design$order != 1),])
downstream.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order != 1), which(design$order != 1)])
downstream.dists$phylo <- as.dist(hja.unifrac[which(design$order != 1), which(design$order != 1)])
downstream.beta.part <- beta.pair(OTUs.PA[which(design$order != 1),], index.family = "sorensen")
downstream.dists$nest <- downstream.beta.part$beta.sne
downstream.dists$turn <- downstream.beta.part$beta.sim
downstream.dists$sor <- downstream.beta.part$beta.sor
downstream.lms <- DDR(dists = downstream.dists, comm = response.matrix)
ddr.summary <- fill.table(downstream.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-downstream-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(downstream.lms)
  dev.off()
}


# Downstream seds
downstream.seds.dists <- list()  
downstream.seds.dists$otus <- vegdist(OTUsREL[which(design$order != 1 & design$habitat == "sediment"),], method = "bray")
downstream.seds.dists$geo <- earth.dist(xy[which(design$order != 1 & design$habitat == "sediment"),])*1000
downstream.seds.dists$den <- as.dist(as.matrix(den.dists)[which(design$order != 1 & design$habitat == "sediment"),which(design$order != 1 & design$habitat == "sediment")])
downstream.seds.dists$env <- dist(env.mat[which(design$order != 1 & design$habitat == "sediment"),])
downstream.seds.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order != 1 & design$habitat == "sediment"),which(design$order != 1 & design$habitat == "sediment")])
downstream.seds.dists$phylo <- as.dist(hja.unifrac[which(design$order != 1 & design$habitat == "sediment"),which(design$order != 1 & design$habitat == "sediment")])
downstream.seds.beta.part <- beta.pair(OTUs.PA[which(design$order != 1 & design$habitat == "sediment"),], index.family = "sorensen")
downstream.seds.dists$nest <- downstream.seds.beta.part$beta.sne
downstream.seds.dists$turn <- downstream.seds.beta.part$beta.sim
downstream.seds.dists$sor <- downstream.seds.beta.part$beta.sor
downstream.seds.lms <- DDR(dists = downstream.seds.dists, comm = response.matrix)
ddr.summary <- fill.table(downstream.seds.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-downstream-seds-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(downstream.seds.lms)
  dev.off()
}

# Downstream water
downstream.water.dists <- list()  
downstream.water.dists$otus <- vegdist(OTUsREL[which(design$order != 1 & design$habitat == "water"),], method = "bray")
downstream.water.dists$geo <- earth.dist(xy[which(design$order != 1 & design$habitat == "water"),])*1000
downstream.water.dists$den <- as.dist(as.matrix(den.dists)[which(design$order != 1 & design$habitat == "water"),which(design$order != 1 & design$habitat == "water")])
downstream.water.dists$env <- dist(env.mat[which(design$order != 1 & design$habitat == "water"),])
downstream.water.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order != 1 & design$habitat == "water"),which(design$order != 1 & design$habitat == "water")])
downstream.water.dists$phylo <- as.dist(hja.unifrac[which(design$order != 1 & design$habitat == "water"),which(design$order != 1 & design$habitat == "water")])
downstream.water.beta.part <- beta.pair(OTUs.PA[which(design$order != 1 & design$habitat == "water"),], index.family = "sorensen")
downstream.water.dists$nest <- downstream.water.beta.part$beta.sne
downstream.water.dists$turn <- downstream.water.beta.part$beta.sim
downstream.water.dists$sor <- downstream.water.beta.part$beta.sor
downstream.water.lms <- DDR(dists = downstream.water.dists, comm = response.matrix)
ddr.summary <- fill.table(downstream.water.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-downstream-water-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(downstream.water.lms)
  dev.off()
}

### Catchment Scale DDRs
hja.dists <- list()  
hja.dists$otus <- vegdist(OTUsREL, method = "bray")
hja.dists$geo <- earth.dist(xy)*1000
hja.dists$den <- as.dist(as.matrix(den.dists))
hja.dists$env <- dist(env.mat)
hja.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist))
hja.dists$phylo <- as.dist(hja.unifrac)
hja.beta.part <- beta.pair(OTUs.PA, index.family = "sorensen")
hja.dists$turn <- hja.beta.part$beta.sim
hja.dists$nest <- hja.beta.part$beta.sne
hja.dists$sor <- hja.beta.part$beta.sor
hja.lms <- DDR(dists = hja.dists, comm = response.matrix)
ddr.summary <- fill.table(hja.lms, ddr.summary = ddr.summary)

if(update.plots){
  pdf(file = paste("figures/DDR-hja-",response.matrix,".pdf", sep = ""), bg = "white")
  par(oma = c(1,1,1,1))
  plot.DDRs(hja.lms)
  dev.off()
}

rownames(ddr.summary) <- NULL
write.table(ddr.summary, file = paste("tables/DDR-summary-table-",response.matrix,".txt", sep = ""), sep = ",", row.names = F)
