#source("analysis/InitialSetup.R")
#source("analysis/DistanceCalcs.R")
#source("analysis/Ordination.R")
require(picante)
require(png)
require(grid)

hja.tree <- read.tree(file = "./data/hja_streams.rename.tree")
# hja.unifrac.raw <- read.delim(file = "./data/hja_streams.tree1.weighted.phylip.dist", header = F, skip = 1, row.names = 1)
# colnames(hja.unifrac.raw) <- as.vector(lapply(strsplit(rownames(hja.unifrac.raw)," "), function(x) x[1]))
# rownames(hja.unifrac.raw) <- colnames(hja.unifrac.raw)
# hja.unifrac <- hja.unifrac.raw[which(rownames(hja.unifrac.raw) %in% rownames(OTUs)),
#                                which(rownames(hja.unifrac.raw) %in% rownames(OTUs))]
# hja.unifrac.dist <- as.dist(hja.unifrac)

unifrac.pcoa <- cmdscale(hja.unifrac.dist, eig=TRUE)
unifvar1 <- round(unifrac.pcoa$eig[1] / sum(unifrac.pcoa$eig),3) * 100
unifvar2 <- round(unifrac.pcoa$eig[2] / sum(unifrac.pcoa$eig),3) * 100
unifvar3 <- round(unifrac.pcoa$eig[3] / sum(unifrac.pcoa$eig),3) * 100

png(filename = "./figures/HJA_PCoA_UniFrac.png",
    width = 1200, height = 1200, res = 96*2)
adonis(hja.unifrac.dist ~ design$habitat * design$watershed + design$order, permutations = 999)
par(mar = c(5, 5, 3, 2) + 0.1)
plot(unifrac.pcoa$points[ ,1], unifrac.pcoa$points[ ,2], ylim = c(-0.17, 0.2), xlim = c(-.17, .2),
     xlab = paste("PCoA 1 (", unifvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", unifvar2, "%)", sep = ""), 
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, at = c(-0.4,-.2,0,.2,.4), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, at = c(-0.4,-.2,0,.2,.4), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(unifrac.pcoa$points[which(design$habitat == "sediment"),1], 
       unifrac.pcoa$points[which(design$habitat == "sediment"),2],
       pch=21, cex=2, bg="grey")
points(unifrac.pcoa$points[which(design$habitat == "water"),1], 
       unifrac.pcoa$points[which(design$habitat == "water"),2],
       pch=24, cex=2, bg="white")
legend("topright", c("Water", "Sediment"),
       pt.bg = c("white", "grey"), pch = c(24,21), cex = 1.5, bty = "n")
ordiellipse(unifrac.pcoa, design$habitat, conf = 0.95)
dev.off()
graphics.off()
img <- readPNG("./figures/HJA_PCoA_UniFrac.png")
grid.raster(img)

matched.phylo <- match.phylo.comm(hja.tree, OTUs)

which(rownames(hja.unifrac.raw) %in% rownames(OTUs))



hja.phylostruc <- phylostruct()


hja.cor <- comm.phylo.cor(samp = matched.phylo$comm, phylo = matched.phylo$phy,
                          metric = "cij", null.model = "sample.taxa.labels", runs = 999)
hja.cor$obs.corr.p
hja.cor$obs.rand.p
hja.cor$obs.corr
hja.cor$obs.rank


hja.pd <- pd(samp = matched.phylo$com, tree = matched.phylo$phy, include.root = F)
hja.pd$PDavg <- hja.pd$PD / hja.pd$SR
write.table(hja.pd, file = "./data/hja.pd.txt", sep="\t")
read.table("./data/hja.pd.txt", sep = "\t")
hja.pd.df <- data.frame(cbind(design, hja.pd, env.mat))
t.test(hja.pd.df$PD[which(design$habitat == "water")],
       hja.pd.df$PD[which(design$habitat == "sediment")])
hja.pdsr.lm <- lm(log10(PD)~log10(SR), data=hja.pd)
summary(hja.pdsr.lm)
plot(log10(PD)~log10(SR), data=hja.pd)
abline(hja.pdsr.lm)

source("analysis/DDRs.R")
uf1 <- (as.dist(hja.unifrac[which(design$order == 1), which(design$order == 1)]))
uf2 <- (as.dist(hja.unifrac[which(design$order == 2), which(design$order == 2)]))
uf3 <- (as.dist(hja.unifrac[which(design$order == 3), which(design$order == 3)]))
uf4 <- (as.dist(hja.unifrac[which(design$order == 4), which(design$order == 4)]))
uf5 <- (as.dist(hja.unifrac[which(design$order == 5), which(design$order == 5)]))


require(pez)


### subset tree and species matrix
n <- 100
top.taxa = matrix(nrow = nrow(OTUsREL), ncol = n+1)
rownames(top.taxa) <- rownames(OTUsREL)
for(i in 1:nrow(OTUsREL)){
  top.taxa[i,1:n] <- names(sort(tail(sort(OTUsREL[i,]),n), decreasing = T))
  top.taxa[i,n+1] <- sum(OTUsREL[i,top.taxa[i,1:n]])
}

OTUsREL.top <- OTUsREL[,which(colnames(OTUsREL) %in% top.taxa)]

hja.phylo <- match.phylo.comm(phy = hja.tree, comm = OTUsREL.top)
hja.mntd.ses.tip <- ses.mntd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "taxa.labels")
hja.mpd.ses.tip <- ses.mpd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "taxa.labels")
hja.pd.ses.tip <- ses.pd(samp = hja.phylo$comm, tree = hja.phylo$phy, null.model = "taxa.labels")
hja.mntd.ses.ind <- ses.mntd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "independentswap")
hja.mpd.ses.ind <- ses.mpd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "independentswap")
hja.pd.ses.ind <- ses.pd(samp = hja.phylo$comm, tree = hja.phylo$phy, null.model = "independentswap")
hja.mntd.ses.trial <- ses.mntd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "trialswap")
hja.mpd.ses.trial <- ses.mpd(samp = hja.phylo$comm, dis = cophenetic(hja.phylo$phy), null.model = "trialswap")
hja.pd.ses.trial <- ses.pd(samp = hja.phylo$comm, tree = hja.phylo$phy, null.model = "trialswap")

hja.phylo.div <- NULL
hja.phylo.div <- list(
mntd = list(trial=hja.mntd.ses.trial,tip=hja.mntd.ses.tip,ind=hja.mntd.ses.ind),
pd = list(trial=hja.pd.ses.trial,tip=hja.pd.ses.tip,ind=hja.pd.ses.ind),
mpd = list(trial=hja.mpd.ses.trial,tip=hja.mpd.ses.tip,ind=hja.mpd.ses.ind))


write.table(hja.mntd.ses.trial, file = "./data/hja-mntd-ses-trial.txt", sep="\t")
write.table(hja.mpd.ses.trial, file = "./data/hja-mpd-ses-trial.txt", sep="\t")
write.table(hja.pd.ses.trial, file = "./data/hja-pd-ses-trial.txt", sep="\t")
write.table(hja.mntd.ses.tip, file = "./data/hja-mntd-ses-tip.txt", sep="\t")
write.table(hja.mpd.ses.tip, file = "./data/hja-mpd-ses-tip.txt", sep="\t")
write.table(hja.pd.ses.tip, file = "./data/hja-pd-ses-tip.txt", sep="\t")
write.table(hja.mntd.ses.ind, file = "./data/hja-mntd-ses-ind.txt", sep="\t")
write.table(hja.mpd.ses.ind, file = "./data/hja-mpd-ses-ind.txt", sep="\t")
write.table(hja.pd.ses.ind, file = "./data/hja-pd-ses-ind.txt", sep="\t")  
  
read.table(file = "./data/hja-mntd-ses-trial.txt", sep="\t")
read.table(file = "./data/hja-mpd-ses-trial.txt", sep="\t")
read.table(file = "./data/hja-pd-ses-trial.txt", sep="\t")
read.table(file = "./data/hja-mntd-ses-tip.txt", sep="\t")
read.table(file = "./data/hja-mpd-ses-tip.txt", sep="\t")
read.table(file = "./data/hja-pd-ses-tip.txt", sep="\t")
read.table(file = "./data/hja-mntd-ses-ind.txt", sep="\t")
read.table(file = "./data/hja-mpd-ses-ind.txt", sep="\t")
read.table(file = "./data/hja-pd-ses-ind.txt", sep="\t")  

plot(hja.phylo.div$mntd$trial$mntd.obs.z[which(design$habitat!="water")] 
     ~ design$order[which(design$habitat!="water")])
abline(h=0)
plot(hja.phylo.div$mntd$trial$mntd.obs.z[which(design$habitat=="water")] 
     ~ design$order[which(design$habitat=="water")])
abline(h=0)

plot(hja.phylo.div$mpd$trial$ntaxa, hja.phylo.div$mpd$trial$mpd.obs.z)
abline(h=0)
summary(lm(hja.phylo.div$mntd$trial$mntd.obs.z ~ hja.phylo.div$mntd$trial$ntaxa))
