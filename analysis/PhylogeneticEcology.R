source("analysis/InitialSetup.R")
source("analysis/DistanceCalcs.R")
source("analysis/Ordination.R")
require(picante)
require(png)
require(grid)

hja.tree <- read.tree(file = "./data/hja_streams.rename.tree")
hja.unifrac.raw <- read.delim(file = "./data/hja_streams.tree1.weighted.phylip.dist", header = F, skip = 1, row.names = 1)
colnames(hja.unifrac.raw) <- as.vector(lapply(strsplit(rownames(hja.unifrac.raw)," "), function(x) x[1]))
rownames(hja.unifrac.raw) <- colnames(hja.unifrac.raw)
hja.unifrac <- hja.unifrac.raw[which(rownames(hja.unifrac.raw) %in% rownames(OTUs)), 
                               which(rownames(hja.unifrac.raw) %in% rownames(OTUs))]
hja.unifrac.dist <- as.dist(hja.unifrac)

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


source("analysis/DDRs.R")
uf1 <- (as.dist(hja.unifrac[which(design$order == 1), which(design$order == 1)]))
uf2 <- (as.dist(hja.unifrac[which(design$order == 2), which(design$order == 2)]))
uf3 <- (as.dist(hja.unifrac[which(design$order == 3), which(design$order == 3)]))
uf4 <- (as.dist(hja.unifrac[which(design$order == 4), which(design$order == 4)]))
uf5 <- (as.dist(hja.unifrac[which(design$order == 5), which(design$order == 5)]))
