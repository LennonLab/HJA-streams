#### Sed-water comparisons within sites

paired.div <- alpha.div
paired.OTUsREL.log <- OTUsREL.log
remove.sites <- 0
for(site in unique(paired.div$site)){
  if(sum(paired.div$site == site) < 2){
    paired.div <- paired.div[-which(paired.div$site == site),]
    remove.sites <- c(remove.sites, which(design$site == site))
  }
}
paired.OTUsREL.log <- paired.OTUsREL.log[-remove.sites,]

paired.div.difs <- (paired.div$S.rare[which(paired.div$habitat == "water")] -
                      paired.div$S.rare[which(paired.div$habitat == "sediment")])
mean(paired.div.difs)
se(paired.div.difs)
t.test(paired.div$S.rare[which(paired.div$habitat == "water")],
       paired.div$S.rare[which(paired.div$habitat == "sediment")], paired = TRUE)

#### Figure 2: Water vs. Sed. diversity
png(filename = "./figures/Figure2.png",
    width = 1200, height = 1200, res = 96*2)
par(c(5,5,3,2) + 0.3)
barplot(height = sort(paired.div.difs, decreasing = T), names.arg = unique(paired.div$site),
        horiz = TRUE, yaxt = "n", xaxt = "n", xlim = c(-1500,1500),
        col = "black", border = NA)
axis(side = 1, labels = T, lwd = 2, lwd.ticks = 2)
axis(side = 3, labels = F, lwd = 2, tck = 0)
abline(v = 0, lty = 3)
mtext(side = 1, line = 3, cex = 1.5, expression(paste("Difference in ",alpha,"-diversity")))
mtext(side = 3, line = 0.5, at = c(-1500), adj = 0, "Higher in\nsediment")
mtext(side = 3, line = 0.5, at = c(1500), adj = 1, "Higher in\nwater")
#mtext(side = 2, line = 2.5, cex = 1.2, "Site")
dev.off()
graphics.off()

img <- readPNG("./figures/Figure2.png")
grid.raster(img)


### Organismal analysis
top.taxa = matrix(nrow = nrow(OTUsREL), ncol = 5)
rownames(top.taxa) <- rownames(OTUsREL)
for(i in 1:nrow(OTUsREL)){
  top.taxa[i,] <- names(sort(tail(sort(OTUsREL[i,]),5), decreasing = T))
}


OTU.tax[OTU.tax$OTU %in% unique(as.vector(top.taxa[,5])),]

OTU.tax[OTU.tax$OTU %in% unique(as.vector(top.taxa[which(design$habitat == "water")])),]
OTU.tax[OTU.tax$OTU %in% unique(as.vector(top.taxa[which(design$habitat == "sediment")])),]

taxa.summary <- as.data.frame(cbind(rowSums(sediment.only > 1), rowSums(water.only > 1), 
                                    rowSums(water.sed > 1), rowSums(in.sediment > 1),
                                    rowSums(in.water > 1), rowSums(OTUs.PA)))
colnames(taxa.summary) <- c("sediments.only", "water.only", "shared", "in.sediment", "in.water", "total")

taxa.summary$water.only[which(design$habitat == "water")] / taxa.summary$total[which(design$habitat == "water")]
taxa.summary$sediments.only[which(design$habitat == "sediment")] / taxa.summary$total[which(design$habitat == "sediment")]
taxa.summary$shared / taxa.summary$total

# Unique to stream order
order.1.OTUs <- OTUs[, which(colSums(OTUs[which(design$order != "1"),]) == 0)]
order.2.OTUs <- OTUs[, which(colSums(OTUs[which(design$order != "2"),]) == 0)]
order.3.OTUs <- OTUs[, which(colSums(OTUs[which(design$order != "3"),]) == 0)]
order.4.OTUs <- OTUs[, which(colSums(OTUs[which(design$order != "4"),]) == 0)]
order.5.OTUs <- OTUs[, which(colSums(OTUs[which(design$order != "5"),]) == 0)]

unique(OTU.tax[OTU.tax$OTU %in% unique(as.vector(colnames(order.1.OTUs))),]$Order)
unique(OTU.tax[OTU.tax$OTU %in% unique(as.vector(colnames(order.2.OTUs))),]$Order)
unique(OTU.tax[OTU.tax$OTU %in% unique(as.vector(colnames(order.3.OTUs))),]$Order)
unique(OTU.tax[OTU.tax$OTU %in% unique(as.vector(colnames(order.4.OTUs))),]$Order)
unique(OTU.tax[OTU.tax$OTU %in% unique(as.vector(colnames(order.5.OTUs))),]$Order)

# site.by.phylum <- matrix(nrow = length(rownames(OTUsREL)), ncol = length(unique(OTU.tax[OTU.tax$OTU %in% colnames(OTUsREL),]$Phylum)))
# rownames(site.by.phylum) <- rownames(OTUsREL)
# colnames(site.by.phylum) <- unique(phyla.present)
# 
# # Taxonomic Analysis
# for(site in 1:length(rownames(OTUsREL))){
#   comm <- OTUsREL[site,]
#   taxa.present <- rownames(as.matrix(comm))
#   phyla.present <- (OTU.tax[OTU.tax$OTU %in% taxa.present,]$Phylum)
# 
#   
#   for(phylum in unique(phyla.present)){
#     site.by.phylum[site,phylum] <- sum(phyla.present == phylum) / length(phyla.present)
#   }
#   rownames(OTUsREL)[site]
# }
### Common Bacteria
sediment.common <- OTUsREL[, which(colSums(OTUs.PA[which(design$habitat == "sediment"),]) > 0.4*length(rownames(OTUs)))]
water.common <- OTUsREL[, which(colSums(OTUs.PA[which(design$habitat == "water"),]) > 0.4*length(rownames(OTUs)))]
intersect(colnames(sediment.common),colnames(water.common))

OTU.tax[OTU.tax$OTU %in% intersect(colnames(sediment.common),colnames(water.common)),]

### Mantel Test
# Water Distances
water.lat <- as.numeric(env$latitude[which(design$habitat == "water")])
water.lon <- as.numeric(env$longitude[which(design$habitat == "water")])

water.struc.dist <- as.matrix(1 - vegdist(OTUsREL.log[which(design$habitat == "water"),]))
water.coord.dist <- as.matrix(dist(as.matrix(water.lat, water.lon)))

mantel(water.struc.dist, water.coord.dist, method = "pearson", permutations = 999)

# Sediment Distances
sed.lat <- as.numeric(env$latitude[which(design$habitat == "sediment")])
sed.lon <- as.numeric(env$longitude[which(design$habitat == "sediment")])

sed.struc.dist <- as.matrix(1 - vegdist(OTUsREL.log[which(design$habitat == "sediment"),]))
sed.coord.dist <- as.matrix(dist(as.matrix(sed.lat, sed.lon)))

mantel(sed.struc.dist, sed.coord.dist, method = "pearson", permutations = 999)

#### Figure S1: Surface-water / Sediment differences are robust across scales
png(filename = "./figures/FigureS1.png",
    width = 2400, height = 1200, res = 96*2)
par(mfrow = c(1,2))

# Lookout Creek PCoA
par(mar = c(5, 5, 3, 2) + 0.1)
plot(lc.pcoa$points[ ,1], lc.pcoa$points[ ,2], ylim = c(-0.4, 0.35),
     xlab = paste("PCoA 1 (", lc.var1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", lc.var2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
mtext("Lookout Creek Watershed", side = 3, line = 2, cex = 1.5)
points(lc.pcoa$points[which(lc.design$habitat == "water"),1], 
       lc.pcoa$points[which(lc.design$habitat == "water"),2],
       pch=21, cex=2, bg="skyblue")
points(lc.pcoa$points[which(lc.design$habitat == "sediment"),1], 
       lc.pcoa$points[which(lc.design$habitat == "sediment"),2],
       pch=24, cex=2, bg="wheat")
legend("topleft", c("Sediment", "Water"),
       pt.bg = c("wheat", "skyblue"), pch = c(21,24), cex = 1.5, bty = "n")


# WS01 PCoA
par(mar = c(5, 5, 3, 2) + 0.1)
plot(w1.pcoa$points[ ,1], w1.pcoa$points[ ,2], ylim = c(-0.4, 0.35),
     xlab = paste("PCoA 1 (", w1.var1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", w1.var2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
mtext("Watershed 01", side = 3, line = 2, cex = 1.5)
points(w1.pcoa$points[which(w1.design$habitat == "water"),1], 
       w1.pcoa$points[which(w1.design$habitat == "water"),2],
       pch=21, cex=2, bg="skyblue")
points(w1.pcoa$points[which(w1.design$habitat == "sediment"),1], 
       w1.pcoa$points[which(w1.design$habitat == "sediment"),2],
       pch=24, cex=2, bg="wheat")

dev.off()
graphics.off()

img <- readPNG("./figures/FigureS1.png")
grid.raster(img)


### CCA across watersheds and habitats


# LC Watershed
lc.cca <- vegan::cca(OTUsREL.log[which(design$watershed == "LC"),] ~ env.mat[which(design$watershed == "LC"),]) 
#anova(lc.cca, by = "axis")
lc.cca.fit <- envfit(lc.cca, env.mat[which(design$watershed == "LC"),], perm = 999)
lc.cca.explainvar1 <- round(lc.cca$CCA$eig[1] /
                              sum(c(lc.cca$CCA$eig, lc.cca$CA$eig)), 3) * 100
lc.cca.explainvar2 <- round(lc.cca$CCA$eig[2] /
                              sum(c(lc.cca$CCA$eig, lc.cca$CA$eig)), 3) * 100

# WS01
w1.cca <- vegan::cca(OTUsREL.log[which(design$watershed == "WS01"),] ~ env.mat[which(design$watershed == "WS01"),]) 
#anova(w1.cca, by = "axis")
w1.cca.fit <- envfit(w1.cca, env.mat[which(design$watershed == "WS01"),], perm = 999)
w1.cca.explainvar1 <- round(w1.cca$CCA$eig[1] /
                              sum(c(w1.cca$CCA$eig, w1.cca$CA$eig)), 3) * 100
w1.cca.explainvar2 <- round(w1.cca$CCA$eig[2] /
                              sum(c(w1.cca$CCA$eig, w1.cca$CA$eig)), 3) * 100

# Water Only
water.cca <- vegan::cca(OTUsREL.log[which(design$habitat == "water"),] ~ env.mat[which(design$habitat == "water"),]) 
#anova(water.cca, by = "axis")
water.cca.fit <- envfit(water.cca, env.mat[which(design$habitat == "water"),], perm = 999)
water.cca.explainvar1 <- round(water.cca$CCA$eig[1] /
                                 sum(c(water.cca$CCA$eig, water.cca$CA$eig)), 3) * 100
water.cca.explainvar2 <- round(water.cca$CCA$eig[2] /
                                 sum(c(water.cca$CCA$eig, water.cca$CA$eig)), 3) * 100

# Sediments Only
sediment.cca <- vegan::cca(OTUsREL.log[which(design$habitat == "sediment"),] ~ env.mat[which(design$habitat == "sediment"),]) 
#anova(sediment.cca, by = "axis")
sediment.cca.fit <- envfit(sediment.cca, env.mat[which(design$habitat == "sediment"),], perm = 999)
sediment.cca.explainvar1 <- round(sediment.cca$CCA$eig[1] /
                                    sum(c(sediment.cca$CCA$eig, sediment.cca$CA$eig)), 3) * 100
sediment.cca.explainvar2 <- round(sediment.cca$CCA$eig[2] /
                                    sum(c(sediment.cca$CCA$eig, sediment.cca$CA$eig)), 3) * 100


#### Figure S2: Environmental variables explain little variation in structure. Elevation seems important.
png(filename = "./figures/FigureS2.png",
    width = 2400, height = 1800, res=96*2)
par(mfrow=c(2,2))

# Lookout Plot
par(mar = c(5, 5, 4, 4) + 0.1)
plot(scores(lc.cca, display = "wa"), xlim = c(-2, 3), ylim = c(-2.5, 2),
     xlab = paste("CCA 1 (", lc.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2 (", lc.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE,
     main="Lookout Creek")
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
abline(h=0, v=0, lty=3)
box(lwd=2)
points(scores(lc.cca, display = "wa"),
       pch=19, cex=2, bg="gray", col="gray")
text(scores(lc.cca, display = "wa"), cex=0.5,
     labels = row.names(scores(lc.cca, display = "wa")))

lc.vectors <- scores(lc.cca, display = "bp")
row.names(lc.vectors) <- c("elev", "temp", "cond", "pH", "TN", "TP")
arrows(0, 0, lc.vectors[, 1] * 2, lc.vectors[, 2]*2,
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(lc.vectors[, 1] * 2, lc.vectors[, 2] * 2, pos=3,
     labels = row.names(lc.vectors), col = "red")

# WS01 Plot
par(mar = c(5, 5, 4, 4) + 0.1)
plot(scores(w1.cca, display = "wa"), xlim = c(-1,3), ylim = c(-3, 2),
     xlab = paste("CCA 1 (", w1.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2 (", w1.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE,
     main = "Watershed 1")
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
abline(h=0, v=0, lty=3)
box(lwd=2)
points(scores(w1.cca, display = "wa"),
       pch=19, cex=2, bg="gray", col="gray")
text(scores(w1.cca, display = "wa"), cex = 0.5,
     labels = row.names(scores(w1.cca, display = "wa")))

w1.vectors <- scores(w1.cca, display = "bp")
row.names(w1.vectors) <- c("elev", "temp", "cond", "pH", "TN", "TP")
arrows(0, 0, w1.vectors[, 1] * 2, w1.vectors[, 2]*2,
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(w1.vectors[, 1] * 2, w1.vectors[, 2] * 2, pos=3,
     labels = row.names(w1.vectors), col = "red")


#Water PLot
par(mar = c(5, 5, 4, 4) + 0.1)
plot(scores(water.cca, display = "wa"), xlim = c(-2, 2), ylim = c(-3, 2.5),
     xlab = paste("CCA 1 (", water.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2 (", water.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE,
     main = "Water")
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
abline(h=0, v=0, lty=3)
box(lwd=2)
points(scores(water.cca, display = "wa"),
       pch=19, cex=2, bg="gray", col="gray")
text(scores(water.cca, display = "wa"), cex = 0.5,
     labels = row.names(scores(water.cca, display = "wa")))

water.vectors <- scores(water.cca, display = "bp")
row.names(water.vectors) <- c("elev", "temp", "cond", "pH", "TN", "TP")
arrows(0, 0, water.vectors[, 1] * 2, water.vectors[, 2]*2,
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(water.vectors[, 1] * 2, water.vectors[, 2] * 2, pos=3,
     labels = row.names(water.vectors), col = "red")

# Sed plot

par(mar = c(5, 5, 4, 4) + 0.1)
plot(scores(sediment.cca, display = "wa"), xlim = c(-2, 3), ylim = c(-1.25, 1.25),
     xlab = paste("CCA 1 (", sediment.cca.explainvar1, "%)", sep = ""), 
     ylab = paste("CCA 2 (", sediment.cca.explainvar2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE,
     main = "Sediment")
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
abline(h=0, v=0, lty=3)
box(lwd=2)
points(scores(sediment.cca, display = "wa"),
       pch=19, cex=2, bg="gray", col="gray")
text(scores(sediment.cca, display = "wa"), cex = 0.5,
     labels = row.names(scores(sediment.cca, display = "wa")))

sediment.vectors <- scores(sediment.cca, display = "bp")
row.names(sediment.vectors) <- c("elev", "temp", "cond", "pH", "TN", "TP")
arrows(0, 0, sediment.vectors[, 1] * 2, sediment.vectors[, 2]*2,
       lwd = 2, lty = 1, length = 0.2, col = "red")
text(sediment.vectors[, 1] * 2, sediment.vectors[, 2] * 2, pos=3,
     labels = row.names(sediment.vectors), col = "red")

dev.off()
graphics.off()

img <- readPNG("./figures/FigureS2.png")
grid.raster(img)





#### Figure 6: Surface Water vs. Sediment Communities

png(filename = "./figures/Figure6.png",
    width = 2400, height = 1200, res = 96*2)
par(mfrow = c(1,2))

# Sediments Only PCoA
par(mar = c(5, 5, 3, 2) + 0.1)
plot(sediment.pcoa$points[ ,1], sediment.pcoa$points[ ,2], ylim = c(-0.4, 0.4), xlim = c(-.4,.4),
     xlab = paste("PCoA 1 (", s.var1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", s.var2, "%)", sep = ""), 
     pch = 19, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
mtext(side = 3, line = 2, "Sediments", cex = 1.5)
points(sediment.pcoa$points[which(sed.design$watershed == "LC"),1], 
       sediment.pcoa$points[which(sed.design$watershed == "LC"),2],
       pch=21, cex=2, bg="slategray2")
points(sediment.pcoa$points[which(sed.design$watershed == "WS01"),1], 
       sediment.pcoa$points[which(sed.design$watershed == "WS01"),2],
       pch=24, cex=2, bg="yellowgreen")
legend("topleft", c("Lookout Creek", "Watershed 01"),
       pt.bg = c("slategray2", "yellowgreen"), pch = c(21,24), cex = 1.5, bty = "n")

# Water Only PCoA
par(mar = c(5, 5, 3, 2) + 0.1)
plot(water.pcoa$points[ ,1], water.pcoa$points[ ,2], ylim = c(-0.4, 0.4), xlim = c(-.5,.4),
     xlab = paste("PCoA 1 (", w.var1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", w.var2, "%)", sep = ""), 
     pch = 16, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
mtext(side = 3, line = 2, "Water", cex = 1.5)
points(water.pcoa$points[which(water.design$watershed == "LC"),1], 
       water.pcoa$points[which(water.design$watershed == "LC"),2],
       pch=21, cex=2, bg="slategray2")
points(water.pcoa$points[which(water.design$watershed == "WS01"),1], 
       water.pcoa$points[which(water.design$watershed == "WS01"),2],
       pch=24, cex=2, bg="yellowgreen")
dev.off()
graphics.off()

img <- readPNG("./figures/Figure6.png")
grid.raster(img)


# Co-occurrence analysis
require(cooccur)
hja.cooccur <- cooccur(mat = as.data.frame(OTUs.PA),
                       type = "site_spp",
                       spp_names = T,
                       thresh = T)


# Rank-Abundance Curve

RAC <- function(x = ""){
  x = as.vector(x)
  x.ab = x[x > 0]
  x.ab.ranked = x.ab[order(x.ab, decreasing = TRUE)]
  return(x.ab.ranked)
}

sed.sites <- OTUsREL.log[which(design$habitat == "sediment"),]
water.sites <- OTUsREL.log[which(design$habitat == "water"),]
rac <- RAC(sed.sites[1,])
ranks <- as.vector(seq(1, length(rac)))
par(mar=c(5.1, 5.1, 4.1, 2.1))                        
plot(ranks, log(rac), type = 'l', axes=F, col = "brown",           
     xlab = "Rank in abundance", ylab = "Abundance",
     las = 1, cex.lab = 1.4, cex.axis = 1.25,
     ylim = c(0, log(20)), xlim = c(0, max(alpha.div$S.obs)))
for(site in 2:nrow(sed.sites)){
  rac <- RAC(sed.sites[site,])
  ranks <- as.vector(seq(1, length(rac)))
  points(ranks, log(rac), type = 'l', col = "brown")
}
for(site in 1:nrow(water.sites)){
  rac <- RAC(water.sites[site,])
  ranks <- as.vector(seq(1, length(rac)))
  points(ranks, log(rac), type = 'l', col = "blue")
}

box()                                                 
axis(side = 1, labels=T, cex.axis = 1.25)             
axis(side = 2, las = 1, cex.axis = 1.25,              
     labels=c(1, 2, 5, 10, 20), at=log(c(1, 2, 5, 10, 20)))



# Elements of Metacommunity Structure (Leibold and Mikkelson 2002)

# water.null <- nullmodel(OTUs.PA[which(design$habitat == "water"),], method = "r0_ind")
# sed.null <- nullmodel(OTUs.PA[which(design$habitat == "sediment"),], method = "r0_ind")
# 
# OTUs.common <- OTUsREL.log * 0
# for(site in 1:nrows(OTUsREL.log)){
#   this.site <- OTUsREL.log[site,]
#   for otu in sort(this.site){
#     rel.abund 
#   }
# }



# ### GIS and Imagery
# DEM <- raster("~/GitHub/HJA-streams/imagery/bare_earth_DEM_1m/gi01001.e00")
# DEM
# image(DEM, col=terrain.colors(1000))
# streams <- readOGR(dsn = "./imagery/stream_network_2008/", layer = "lidar_stream")
# streams@data
# plot(streams, add = TRUE)
# summary(streams)

# heatmaps and clustering
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", 
                                 "#7FFF7F", "yellow", "#FF7F00", "red", 
                                 "#7F0000"))
order <- c(design$site[design$habitat == "water"], design$site[design$habitat == "sediment"])
levelplot(as.matrix(hja.db)[order,order], aspect = "iso", col.regions = jet.colors,
          xlab = "site", ylab = "site", scales = list(cex = 0.5),
          main = "Bray-Curtis Distance")

hja.ward <- hclust(hja.db, method = "ward.D2")
par(mar = c(1,5,2,2) + 0.1)
plot(hja.ward)

require(gplots)
heatmap.2(as.matrix(hja.db), distfun = function(x) vegdist(x, method = "bray"),
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          col = jet.colors(100), trace = "none", density.info = "none")
