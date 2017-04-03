source("analysis/InitialSetup.R")
source("analysis/DistanceCalcs.R")
source("analysis/Ordination.R")
require(picante)
require(png)
require(grid)

#------------------------------------------------#
# Read TREE
hja.tree <- read.tree(file = "./data/hja_streams.rename.tree")

#-------------------------------------------------#
# UniFac Distances 

# unifrac.pcoa <- cmdscale(hja.unifrac.dist, eig=TRUE)
# unifvar1 <- round(unifrac.pcoa$eig[1] / sum(unifrac.pcoa$eig),3) * 100
# unifvar2 <- round(unifrac.pcoa$eig[2] / sum(unifrac.pcoa$eig),3) * 100
# unifvar3 <- round(unifrac.pcoa$eig[3] / sum(unifrac.pcoa$eig),3) * 100
# 
# png(filename = "./figures/HJA_PCoA_UniFrac.png",
#     width = 1200, height = 1200, res = 96*2)
# adonis(hja.unifrac.dist ~ design$habitat * design$watershed + design$order, permutations = 999)
# par(mar = c(5, 5, 3, 2) + 0.1)
# plot(unifrac.pcoa$points[ ,1], unifrac.pcoa$points[ ,2], ylim = c(-0.17, 0.2), xlim = c(-.17, .2),
#      xlab = paste("PCoA 1 (", unifvar1, "%)", sep = ""),
#      ylab = paste("PCoA 2 (", unifvar2, "%)", sep = ""), 
#      pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
# axis(side = 1, labels = T, at = c(-0.4,-.2,0,.2,.4), lwd.ticks = 2, cex.axis = 1.2, las = 1)
# axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
# axis(side = 3, labels = F, at = c(-0.4,-.2,0,.2,.4), lwd.ticks = 2, cex.axis = 1.2, las = 1)
# axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
# abline(h = 0, v = 0, lty = 3)
# box(lwd = 2)
# points(unifrac.pcoa$points[which(design$habitat == "sediment"),1], 
#        unifrac.pcoa$points[which(design$habitat == "sediment"),2],
#        pch=21, cex=2, bg="grey")
# points(unifrac.pcoa$points[which(design$habitat == "water"),1], 
#        unifrac.pcoa$points[which(design$habitat == "water"),2],
#        pch=24, cex=2, bg="white")
# legend("topright", c("Water", "Sediment"),
#        pt.bg = c("white", "grey"), pch = c(24,21), cex = 1.5, bty = "n")
# ordiellipse(unifrac.pcoa, design$habitat, conf = 0.95)
# dev.off()
# graphics.off()
# img <- readPNG("./figures/HJA_PCoA_UniFrac.png")
# grid.raster(img)

# #----------------------------------------------------#
# OTUs.water <- OTUs[which(design$habitat == "water"),]
# OTUs.water <- OTUs.water[,which(colSums(OTUs.water) < 2)]
# OTUs.water <- decostand(OTUs.water, method = 'total')
# OTUs.sed <- OTUs[which(design$habitat == "sediment"),]
# OTUs.sed <- OTUs.sed[,which(colSums(OTUs.sed) < 2)]

#----------------------------------------------------#
# OTUs.water <- OTUs[which(design$habitat == "water"),]
# OTUs.water <- OTUs.water[,-which(colSums(OTUs.water) < 5)]
# OTUs.water <- decostand(OTUs.water, method = 'total')
# OTUs.sed <- OTUs[which(design$habitat == "sediment"),]
# OTUs.sed <- OTUs.sed[,-which(colSums(OTUs.sed) < 4)]
# OTUs.sed <- decostand(OTUs.sed, method = 'total')
# matched.phylo.water <- match.phylo.comm(hja.tree, OTUs.water)
# matched.phylo.sed <- match.phylo.comm(hja.tree, OTUs.sed)
# saveRDS(matched.phylo.water, file = "data/matched.phylo.water.rda")
# saveRDS(matched.phylo.sed, file = "data/matched.phylo.sed.rda")
matched.phylo.water <- readRDS(file = "data/matched.phylo.water.rda")
matched.phylo.sed <- readRDS(file = "data/matched.phylo.sed.rda")


# hja.cor <- comm.phylo.cor(samp = matched.phylo$comm, phylo = matched.phylo$phy,
#                           metric = "cij", null.model = "sample.taxa.labels", runs = 999)
# hja.cor$obs.corr.p
# hja.cor$obs.rand.p
# hja.cor$obs.corr
# hja.cor$obs.rank
# 

# hja.pd <- pd(samp = matched.phylo$com, tree = matched.phylo$phy, include.root = F)
# hja.pd$PDavg <- hja.pd$PD / hja.pd$SR
# write.table(hja.pd, file = "./data/hja.pd.txt", sep="\t")
# hja.pd <- read.table("./data/hja.pd.txt", sep = "\t")
# hja.pd.df <- data.frame(cbind(design, hja.pd, env.mat))
# t.test(hja.pd.df$PD[which(design$habitat == "water")],
#        hja.pd.df$PD[which(design$habitat == "sediment")])
# hja.pdsr.lm <- lm(log10(PD)~log10(SR), data=hja.pd)
# summary(hja.pdsr.lm)
# plot(log10(PD)~log10(SR), data=hja.pd)
# abline(hja.pdsr.lm)
# 
# source("analysis/DDRs.R")
# uf1 <- (as.dist(hja.unifrac[which(design$order == 1), which(design$order == 1)]))
# uf2 <- (as.dist(hja.unifrac[which(design$order == 2), which(design$order == 2)]))
# uf3 <- (as.dist(hja.unifrac[which(design$order == 3), which(design$order == 3)]))
# uf4 <- (as.dist(hja.unifrac[which(design$order == 4), which(design$order == 4)]))
# uf5 <- (as.dist(hja.unifrac[which(design$order == 5), which(design$order == 5)]))
# 
# 
# require(pez)


### subset tree and species matrix
# n <- 100
# top.taxa = matrix(nrow = nrow(OTUsREL), ncol = n+1)
# rownames(top.taxa) <- rownames(OTUsREL)
# for(i in 1:nrow(OTUsREL)){
#   top.taxa[i,1:n] <- names(sort(tail(sort(OTUsREL[i,]),n), decreasing = T))
#   top.taxa[i,n+1] <- sum(OTUsREL[i,top.taxa[i,1:n]])
# }
# 
# OTUsREL.top <- OTUsREL[,which(colnames(OTUsREL) %in% top.taxa)]
# 
# hja.phylo <- match.phylo.comm(phy = hja.tree, comm = OTUsREL.top)
# hja.mntd.ses.tip <- ses.mntd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "taxa.labels")
# hja.mpd.ses.tip <- ses.mpd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "taxa.labels")
# hja.pd.ses.tip <- ses.pd(samp = hja.phylo$comm, tree = hja.phylo$phy, null.model = "taxa.labels")
# hja.mntd.ses.ind <- ses.mntd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "independentswap")
# hja.mpd.ses.ind <- ses.mpd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "independentswap")
# hja.pd.ses.ind <- ses.pd(samp = hja.phylo$comm, tree = hja.phylo$phy, null.model = "independentswap")
# hja.mntd.ses.trial <- ses.mntd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "trialswap")
# hja.mpd.ses.trial <- ses.mpd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "trialswap")
# hja.pd.ses.trial <- ses.pd(samp = hja.phylo$comm, tree = hja.phylo$phy, null.model = "trialswap")
# 
# hja.phylo.div <- NULL
# hja.phylo.div <- list(
# mntd = list(trial=hja.mntd.ses.trial,tip=hja.mntd.ses.tip,ind=hja.mntd.ses.ind),
# pd = list(trial=hja.pd.ses.trial,tip=hja.pd.ses.tip,ind=hja.pd.ses.ind),
# mpd = list(trial=hja.mpd.ses.trial,tip=hja.mpd.ses.tip,ind=hja.mpd.ses.ind))
# 
# 
# write.table(hja.mntd.ses.trial, file = "./data/hja-mntd-ses-trial.txt", sep="\t")
# write.table(hja.mpd.ses.trial, file = "./data/hja-mpd-ses-trial.txt", sep="\t")
# write.table(hja.pd.ses.trial, file = "./data/hja-pd-ses-trial.txt", sep="\t")
# write.table(hja.mntd.ses.tip, file = "./data/hja-mntd-ses-tip.txt", sep="\t")
# write.table(hja.mpd.ses.tip, file = "./data/hja-mpd-ses-tip.txt", sep="\t")
# write.table(hja.pd.ses.tip, file = "./data/hja-pd-ses-tip.txt", sep="\t")
# write.table(hja.mntd.ses.ind, file = "./data/hja-mntd-ses-ind.txt", sep="\t")
# write.table(hja.mpd.ses.ind, file = "./data/hja-mpd-ses-ind.txt", sep="\t")
# write.table(hja.pd.ses.ind, file = "./data/hja-pd-ses-ind.txt", sep="\t")  
#   
# hja.mntd.ses.trial <- read.table(file = "./data/hja-mntd-ses-trial.txt", sep="\t")
# hja.mpd.ses.trial <- read.table(file = "./data/hja-mpd-ses-trial.txt", sep="\t")
# hja.pd.ses.trial <- read.table(file = "./data/hja-pd-ses-trial.txt", sep="\t")
# hja.mntd.ses.tip <- read.table(file = "./data/hja-mntd-ses-tip.txt", sep="\t")
# hja.mpd.ses.tip <- read.table(file = "./data/hja-mpd-ses-tip.txt", sep="\t")
# hja.pd.ses.tip <- read.table(file = "./data/hja-pd-ses-tip.txt", sep="\t")
# hja.mntd.ses.ind <- read.table(file = "./data/hja-mntd-ses-ind.txt", sep="\t")
# hja.mpd.ses.ind <- read.table(file = "./data/hja-mpd-ses-ind.txt", sep="\t")
# hja.pd.ses.ind <- read.table(file = "./data/hja-pd-ses-ind.txt", sep="\t")  

# hja.phylo.div <- NULL
# hja.phylo.div <- list(
# mntd = list(trial=hja.mntd.ses.trial,tip=hja.mntd.ses.tip,ind=hja.mntd.ses.ind),
# pd = list(trial=hja.pd.ses.trial,tip=hja.pd.ses.tip,ind=hja.pd.ses.ind),
# mpd = list(trial=hja.mpd.ses.trial,tip=hja.mpd.ses.tip,ind=hja.mpd.ses.ind))


# plot(hja.phylo.div$mntd$trial$mntd.obs.z[which(design$habitat!="water")] 
#      ~ design$order[which(design$habitat!="water")])
# abline(h=0)
# plot(hja.phylo.div$mntd$trial$mntd.obs.z[which(design$habitat=="water")] 
#      ~ design$order[which(design$habitat=="water")])
# abline(h=0)
# 
# plot(hja.phylo.div$mpd$trial$ntaxa, hja.phylo.div$mpd$trial$mpd.obs.z)
# abline(h=0)
# summary(lm(hja.phylo.div$mntd$trial$mntd.obs.z ~ hja.phylo.div$mntd$trial$ntaxa))
# 

### bNTI 
water.phy <- matched.phylo.water$phy
water.comm <- matched.phylo.water$comm
sed.phy <- matched.phylo.sed$phy
sed.comm <- matched.phylo.sed$comm

# Calc. obs. mntds
mntds.water <- (comdistnt(water.comm,
                        cophenetic(water.phy),
                        abundance.weighted=T))
mntds.sed <- (comdistnt(sed.comm,
                        cophenetic(sed.phy),
                        abundance.weighted=T))
saveRDS(mntds.water, file = "data/mntds-water.rda")
saveRDS(mntds.sed, file = "data/mntds-sed.rda")
# mntds.water <- readRDS(file = "data/mntds-water.rda")
# mntds.sed <- readRDS(file = "data/mntds-sed.rda")

# Create null comms
mntd.null.water <- NULL
write(paste("i",",","bMNTDvalue.water", sep = ""), file = "./analysis/logs/mntd.water.null.log")
for(i in 1:999){
  temp.mntd <- liste(
    comdistnt(water.comm,
              cophenetic(tipShuffle(water.phy)),
              abundance.weighted=T)
  )[,3]
  mntd.null.water[i] <- mean(temp.mntd)
  write(paste(i,",",mntd.null.water[i], sep = ""), file = "./analysis/logs/mntd.water.null.log", append = T)
}
saveRDS(mntd.null.water, file = "data/mntds-water-null-dist.rda")
#mntd.null.water <- readRDS(file = "data/mntds-water-null-dist.rda")
#hist(mntd.null.water)

mntd.null.sed <- NULL
write(paste("i",",","bMNTDvalue.sed", sep = ""), file = "./analysis/logs/mntd.sed.null.log")
for(i in 1:999){
  temp.mntd <- liste(
    comdistnt(sed.comm,
              cophenetic(tipShuffle(sed.phy)),
              abundance.weighted=T)
  )[,3]
  mntd.null.sed[i] <- mean(temp.mntd)
  write(paste(i,",",mntd.null.sed[i], sep = ""), file = "./analysis/logs/mntd.sed.null.log", append = T)
}
saveRDS(mntd.null.sed, file = "data/mntds-sed-null-dist.rda")
#mntd.null.sed <- readRDS(file = "data/mntds-sed-null-dist.rda")
#hist(mntd.null.sed)

# Calculate bNTIs
# mntds.water.rank <- NULL
# mntds.water.pval <- NULL
# bNTI.water <- NULL
# mntds.sed.rank <- NULL
# mntds.sed.pval <- NULL
# bNTI.sed <- NULL
# i <- 1
# 
# for(each in liste(mntds.water)[,3]){
#   # mntds.rank[i] <- rank(c(each, mntd.null))[1]
#   # pval <- mntds.rank[i] / (length(mntd.null) + 1)
#   # if (pval > 0.5) mntds.pval[i] <- (1-pval)*2
#   # if (pval <= 0.5) mntds.pval[i] <- pval*2
#   # i <- i+1
#   bNTI.water[i] <- (each - mean(mntd.null.water)) / sd(mntd.null.water)
#   i <- i+1
# }
# 
# i <- 1
# for(each in liste(mntds.sed)[,3]){
#   # mntds.rank[i] <- rank(c(each, mntd.null))[1]
#   # pval <- mntds.rank[i] / (length(mntd.null) + 1)
#   # if (pval > 0.5) mntds.pval[i] <- (1-pval)*2
#   # if (pval <= 0.5) mntds.pval[i] <- pval*2
#   # i <- i+1
#   bNTI.sed[i] <- (each - mean(mntd.null.sed)) / sd(mntd.null.sed)
#   i <- i+1
# }
# 
# write.csv()
# length(bNTI.sed) - (sum(bNTI.sed < -2) + sum(bNTI.sed > 2))
# 
# 
# water.den.dists <- liste(as.dist(
#   den.dists[which(design$habitat == "water"), 
#             which(design$habitat == "water")]))[,3]
# sed.den.dists <- liste(as.dist(
#   den.dists[which(design$habitat == "sediment"), 
#             which(design$habitat == "sediment")]))[,3]
# summary(lm(bNTI.sed ~ sed.env.dist.ls * sed.den.dists))
# summary(lm(bNTI.water ~ water.env.dist.ls * water.den.dists))
# plot(water.den.dists, bNTI.water)
# plot(sed.env.dist.ls, bNTI.sed)
# plot(water.env.dist.ls, bNTI.water)
# 
# sed.mntd <- ses.mntd(samp = sed.comm, dis = cophenetic(sed.phy), 
#                      abundance.weighted = T, iterations = 999,
#                       null.model = "taxa.labels")
# water.mntd <- ses.mntd(samp = water.comm, dis = cophenetic(water.phy), 
#                      abundance.weighted = T, iterations = 999,
#                      null.model = "taxa.labels")
