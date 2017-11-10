#source("./analysis/InitialSetup.R")
require('betapart')

ddr.summary <- data.frame()

response.matrix <- "turn"
update.plots <- TRUE

# Water DDR
water.dists <- list()  
water.dists$otus <- 1-vegdist(OTUsREL[which(design$habitat == "water"),], method = "bray")
water.dists$geo <- earth.dist(xy[which(design$habitat == "water"),])*1000
water.dists$den <- as.dist(as.matrix(den.dists)[which(design$habitat == "water"), which(design$habitat == "water")])
water.dists$env <- dist(env.mat[which(design$habitat == "water"),])
water.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$habitat == "water"), which(design$habitat == "water")])
water.dists$phylo <- 1-as.dist(hja.unifrac[which(design$habitat == "water"), which(design$habitat == "water")])
water.beta.part <- beta.pair(OTUs.PA[which(design$habitat == "water"),], index.family = "sorensen")
water.dists$nest <- 1-water.beta.part$beta.sne
water.dists$turn <- 1-water.beta.part$beta.sim
water.dists$sor <- 1-water.beta.part$beta.sor
water.dists$bNTI <- as.dist(normalize.matrix(as.matrix(bNTI.dist)[which(design$habitat == "water"), which(design$habitat == "water")]))
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
sediment.dists$otus <- 1-vegdist(OTUsREL[which(design$habitat == "sediment"),], method = "bray")
sediment.dists$geo <- earth.dist(xy[which(design$habitat == "sediment"),])*1000
sediment.dists$den <- as.dist(as.matrix(den.dists)[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sediment.dists$env <- dist(env.mat[which(design$habitat == "sediment"),])
sediment.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sediment.dists$phylo <- 1-as.dist(hja.unifrac[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
sediment.beta.part <- beta.pair(OTUs.PA[which(design$habitat == "sediment"),], index.family = "sorensen")
sediment.dists$nest <- 1-sediment.beta.part$beta.sne
sediment.dists$turn <- 1-sediment.beta.part$beta.sim
sediment.dists$sor <- 1-sediment.beta.part$beta.sor
sediment.dists$bNTI <- as.dist(normalize.matrix(as.matrix(bNTI.dist)[which(design$habitat == "sediment"), which(design$habitat == "sediment")]))
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
headwater.dists$otus <- 1-vegdist(OTUsREL[which(design$order == 1),], method = "bray")
headwater.dists$geo <- earth.dist(xy[which(design$order == 1),])*1000
headwater.dists$den <- as.dist(as.matrix(den.dists)[which(design$order == 1), which(design$order == 1)])
headwater.dists$env <- dist(env.mat[which(design$order == 1),])
headwater.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order == 1), which(design$order == 1)])
headwater.dists$phylo <- 1-as.dist(hja.unifrac[which(design$order == 1), which(design$order == 1)])
headwater.beta.part <- beta.pair(OTUs.PA[which(design$order == 1),], index.family = "sorensen")
headwater.dists$nest <- 1-headwater.beta.part$beta.sne
headwater.dists$turn <- 1-headwater.beta.part$beta.sim
headwater.dists$sor <- 1-headwater.beta.part$beta.sor
headwater.dists$bNTI <- as.dist(normalize.matrix(as.matrix(bNTI.dist)[which(design$order == 1), which(design$order == 1)]))
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
headwater.seds.dists$otus <- 1-vegdist(OTUsREL[which(design$order == 1 & design$habitat == "sediment"),], method = "bray")
headwater.seds.dists$geo <- earth.dist(xy[which(design$order == 1 & design$habitat == "sediment"),])*1000
headwater.seds.dists$den <- as.dist(as.matrix(den.dists)[which(design$order == 1 & design$habitat == "sediment"),which(design$order == 1 & design$habitat == "sediment")])
headwater.seds.dists$env <- dist(env.mat[which(design$order == 1 & design$habitat == "sediment"),])
headwater.seds.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order == 1 & design$habitat == "sediment"),which(design$order == 1 & design$habitat == "sediment")])
headwater.seds.dists$phylo <- 1-as.dist(hja.unifrac[which(design$order == 1 & design$habitat == "sediment"),which(design$order == 1 & design$habitat == "sediment")])
headwater.seds.beta.part <- beta.pair(OTUs.PA[which(design$order == 1 & design$habitat == "sediment"),], index.family = "sorensen")
headwater.seds.dists$nest <- 1-headwater.seds.beta.part$beta.sne
headwater.seds.dists$turn <- 1-headwater.seds.beta.part$beta.sim
headwater.seds.dists$sor <- 1-headwater.seds.beta.part$beta.sor
headwater.seds.dists$bNTI <- as.dist(normalize.matrix(as.matrix(bNTI.dist)[which(design$order == 1 & design$habitat == "sediment"), which(design$order == 1 & design$habitat == "sediment")]))
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
headwater.water.dists$otus <- 1-vegdist(OTUsREL[which(design$order == 1 & design$habitat == "water"),], method = "bray")
headwater.water.dists$geo <- earth.dist(xy[which(design$order == 1 & design$habitat == "water"),])*1000
headwater.water.dists$den <- as.dist(as.matrix(den.dists)[which(design$order == 1 & design$habitat == "water"),which(design$order == 1 & design$habitat == "water")])
headwater.water.dists$env <- dist(env.mat[which(design$order == 1 & design$habitat == "water"),])
headwater.water.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order == 1 & design$habitat == "water"),which(design$order == 1 & design$habitat == "water")])
headwater.water.dists$phylo <- 1-as.dist(hja.unifrac[which(design$order == 1 & design$habitat == "water"),which(design$order == 1 & design$habitat == "water")])
headwater.water.beta.part <- beta.pair(OTUs.PA[which(design$order == 1 & design$habitat == "water"),], index.family = "sorensen")
headwater.water.dists$nest <- 1-headwater.water.beta.part$beta.sne
headwater.water.dists$turn <- 1-headwater.water.beta.part$beta.sim
headwater.water.dists$sor <- 1-headwater.water.beta.part$beta.sor
headwater.water.dists$bNTI <- as.dist(normalize.matrix(as.matrix(bNTI.dist)[which(design$order == 1 & design$habitat == "water"), which(design$order == 1 & design$habitat == "water")]))
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
downstream.dists$otus <- 1-vegdist(OTUsREL[which(design$order != 1),], method = "bray")
downstream.dists$geo <- earth.dist(xy[which(design$order != 1),])*1000
downstream.dists$den <- as.dist(as.matrix(den.dists)[which(design$order != 1), which(design$order != 1)])
downstream.dists$env <- dist(env.mat[which(design$order != 1),])
downstream.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order != 1), which(design$order != 1)])
downstream.dists$phylo <- 1-as.dist(hja.unifrac[which(design$order != 1), which(design$order != 1)])
downstream.beta.part <- beta.pair(OTUs.PA[which(design$order != 1),], index.family = "sorensen")
downstream.dists$nest <- 1-downstream.beta.part$beta.sne
downstream.dists$turn <- 1-downstream.beta.part$beta.sim
downstream.dists$sor <- 1-downstream.beta.part$beta.sor
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
downstream.seds.dists$otus <- 1-vegdist(OTUsREL[which(design$order != 1 & design$habitat == "sediment"),], method = "bray")
downstream.seds.dists$geo <- earth.dist(xy[which(design$order != 1 & design$habitat == "sediment"),])*1000
downstream.seds.dists$den <- as.dist(as.matrix(den.dists)[which(design$order != 1 & design$habitat == "sediment"),which(design$order != 1 & design$habitat == "sediment")])
downstream.seds.dists$env <- dist(env.mat[which(design$order != 1 & design$habitat == "sediment"),])
downstream.seds.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order != 1 & design$habitat == "sediment"),which(design$order != 1 & design$habitat == "sediment")])
downstream.seds.dists$phylo <- 1-as.dist(hja.unifrac[which(design$order != 1 & design$habitat == "sediment"),which(design$order != 1 & design$habitat == "sediment")])
downstream.seds.beta.part <- beta.pair(OTUs.PA[which(design$order != 1 & design$habitat == "sediment"),], index.family = "sorensen")
downstream.seds.dists$nest <- 1-downstream.seds.beta.part$beta.sne
downstream.seds.dists$turn <- 1-downstream.seds.beta.part$beta.sim
downstream.seds.dists$sor <- 1-downstream.seds.beta.part$beta.sor
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
downstream.water.dists$otus <- 1-vegdist(OTUsREL[which(design$order != 1 & design$habitat == "water"),], method = "bray")
downstream.water.dists$geo <- earth.dist(xy[which(design$order != 1 & design$habitat == "water"),])*1000
downstream.water.dists$den <- as.dist(as.matrix(den.dists)[which(design$order != 1 & design$habitat == "water"),which(design$order != 1 & design$habitat == "water")])
downstream.water.dists$env <- dist(env.mat[which(design$order != 1 & design$habitat == "water"),])
downstream.water.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist)[which(design$order != 1 & design$habitat == "water"),which(design$order != 1 & design$habitat == "water")])
downstream.water.dists$phylo <- 1- as.dist(hja.unifrac[which(design$order != 1 & design$habitat == "water"),which(design$order != 1 & design$habitat == "water")])
downstream.water.beta.part <- beta.pair(OTUs.PA[which(design$order != 1 & design$habitat == "water"),], index.family = "sorensen")
downstream.water.dists$nest <- 1-downstream.water.beta.part$beta.sne
downstream.water.dists$turn <- 1-downstream.water.beta.part$beta.sim
downstream.water.dists$sor <- 1-downstream.water.beta.part$beta.sor
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
hja.dists$otus <- 1-vegdist(OTUsREL, method = "bray")
hja.dists$geo <- earth.dist(xy)*1000
hja.dists$den <- as.dist(as.matrix(den.dists))
hja.dists$env <- dist(env.mat)
hja.dists$rc.bray <- as.dist(as.matrix(RC.bray.dist))
hja.dists$phylo <- 1-as.dist(hja.unifrac)
hja.beta.part <- beta.pair(OTUs.PA, index.family = "sorensen")
hja.dists$turn <- 1-hja.beta.part$beta.sim
hja.dists$nest <- 1-hja.beta.part$beta.sne
hja.dists$sor <- 1-hja.beta.part$beta.sor
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


# multiple regression
summary(lm(water.dists$otus ~ water.dists$den + water.dists$env))
summary(lm(sediment.dists$otus ~ sediment.dists$den + sediment.dists$env))
summary(lm(headwater.dists$otus ~ headwater.dists$den + headwater.dists$env))
summary(lm(headwater.water.dists$otus ~ headwater.water.dists$den + headwater.water.dists$env))
summary(lm(headwater.seds.dists$otus ~ headwater.seds.dists$den + headwater.seds.dists$env))
summary(lm(downstream.dists$otus ~ downstream.dists$den + downstream.dists$env))
summary(lm(downstream.seds.dists$otus ~ downstream.seds.dists$den + downstream.seds.dists$env))
summary(lm(downstream.water.dists$otus ~ downstream.water.dists$den + downstream.water.dists$env))
