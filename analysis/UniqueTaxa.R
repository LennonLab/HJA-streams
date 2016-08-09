# source("./analysis/InitialSetup.R")

sediment.only <- OTUs[, which(colSums(OTUs[which(design$habitat == "water"),]) == 0)]
water.only <- OTUs[, which(colSums(OTUs[which(design$habitat == "sediment"),]) == 0)]
water.sed <- OTUs[, which(colSums(OTUs[which(design$habitat == "sediment"),]) > 1 & 
                            colSums(OTUs[which(design$habitat == "water"),]) > 1)]
in.water <- OTUs[, which(colSums(OTUs[which(design$habitat == "water"),]) > 1)]
in.water.pa <- decostand(in.water, method = "pa")
in.sediment <- OTUs[, which(colSums(OTUs[which(design$habitat == "sediment"),]) > 1)]
in.sediment.pa <- decostand(in.sediment, method = "pa")

water.frac <- as.vector(rowSums(in.water.pa[which(design$habitat == "water"), 
                                            which(colSums(in.water.pa[which(design$habitat == "sediment"),]) > 0)]) 
                        / rowSums(in.water.pa[which(design$habitat == "water"),]))

se(rowSums(in.water[which(design$habitat == "water"), 
                    which(colSums(in.water[which(design$habitat == "sediment"),]) > 1)]) 
   / rowSums(in.water[which(design$habitat == "water"),]))

sed.frac <- as.vector(rowSums(in.sediment.pa[which(design$habitat == "sediment"), 
                                             which(colSums(in.sediment.pa[which(design$habitat == "water"),]) > 1)])
                      / rowSums(in.sediment.pa[which(design$habitat == "sediment"),]))

se(rowSums(in.sediment[which(design$habitat == "sediment"), 
                       which(colSums(in.sediment[which(design$habitat == "water"),]) > 1)]) 
   / rowSums(in.sediment[which(design$habitat == "sediment"),]))

unique.fracs <- as.data.frame(rbind(cbind(1-water.frac, 1),cbind(1-sed.frac, 0)))
colnames(unique.fracs) <- c("unique_frac", "habitat")

capture.output(summary(lm(unique_frac ~ habitat, data = unique.fracs)),
               file = "./tables/unique_compare.txt")

### Figure 2: Proportion unique taxa
png(filename = "./figures/UniqueTaxa.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 2, 2) + 0.1)
boxplot(unique_frac ~ habitat, data = unique.fracs, names = c("Sediments", "Water"), col = c("grey", "white"),
        at = c(1,3), xlim = c(0,4), yaxt = "n", cex.axis = 1.2)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 1, labels = F, at = c(1,3), lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd=2)
mtext("Proportion Unique Taxa", side = 2, line = 3.5, cex = 1.5)
mtext("Habitat Type", side = 1, line = 3, cex = 1.5)
text(x = 1, y = 0.25, "*", cex = 2)
dev.off()
graphics.off()

img <- readPNG("./figures/UniqueTaxa.png")
grid.raster(img)