# script for making figures for hja 2015 project

source("analysis/InitialSetup.R")
source("analysis/Ordination.R")
source("analysis/VariationPartitioning.R")
source("analysis/Map.R")

# FIGURE 1. Site Map of HJA
pdf("figures/Figure1_hja-map.pdf", height = 6, width = 6)
hja.map
dev.off()

# FIGURE 2. Constrained Ordination
pdf("figures/Figure2_hja-dbrda.pdf", height = 6, width = 6, bg = "white")
par(oma = c(.1,.1,.1,.1), pty = "s")
plot(scores(hja.dbrda.env, display = "sites"),
     xlab = paste("dbRDA 1 (", explain.var(hja.dbrda.env, axis = 1), "%)", sep = ""),
     ylab = paste("dbRDA 2 (", explain.var(hja.dbrda.env, axis = 2), "%)", sep = ""),
     pch = 21, cex = 2.0, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = F,
     xlim = c(-1, 2), ylim = c(-2.1, 1.1), asp = 1)
add.axes()
abline(h = 0, v = 0, lty = 3)
points(scores(hja.dbrda.env, display = "sites")[which(design$habitat == "sediment"),],
       pch = 21, bg = "wheat", cex = 2)
points(scores(hja.dbrda.env, display = "sites")[which(design$habitat == "water"),],
       pch = 24, bg = "skyblue", cex = 2)
ordiellipse(hja.dbrda.env, design$habitat, conf = 0.95, kind = "se", lwd=1,
            draw = "polygon", col=c("wheat", "skyblue"), border="black", alpha=63)
vectors <- scores(hja.dbrda.env, display = "bp")[-3,]
arrows(0, 0, 1.5*vectors[,1], 1.5*vectors[,2],
       lwd = 2, lty = 1, length = 0.1)
text(1.5*vectors[,1], 1.5*vectors[,2], pos = c(1, 1), 
     labels = rownames(vectors), offset = 1)
legend("topright", legend = c("Planktonic", "Sediment"), pch = c(24, 21), 
       bty = "n", pt.bg = c("skyblue", "wheat"))
dev.off()

# FIGURE 3. 




# FIGURE 4.
#source("analysis/RaupCrickBC.R")
pdf("figures/Figure4_RaupCrickBC.pdf", width = 6, height = 6)
RC.brayplot <- ggplot(data = rc.dat) +
  geom_density(aes(x = rc.bray, y = ..scaled.., fill = Habitat), alpha = .5) +
  labs(x = expression(paste("RC"["bray"])), y = "Frequency") + 
  scale_fill_manual(values = c("skyblue", "wheat"))+
  theme_cowplot() + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
RC.brayplot
dev.off()

# FIGURE 5
#source("analysis/bNTI_analysis.R")
pdf("figures/Figure5_CommAssembly.pdf", width = 10, height = 6)
assembly.counts
dev.off()
