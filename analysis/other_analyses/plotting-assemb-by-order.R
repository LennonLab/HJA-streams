sed.design <- design[which(design$habitat == "sediment"),]
water.design <- design[which(design$habitat == "water"),]

bnti.hw.sed <- as.dist(
  as.matrix(bNTI.sed.dist)[which(sed.design$order == 1), which(sed.design$order == 1)])
bnti.ds.sed <- as.dist(
  as.matrix(bNTI.sed.dist)[which(sed.design$order != 1), which(sed.design$order != 1)])
bnti.hw.water <- as.dist(
  as.matrix(bNTI.water.dist)[which(water.design$order == 1), which(water.design$order == 1)])
bnti.ds.water <- as.dist(
  as.matrix(bNTI.water.dist)[which(water.design$order != 1), which(water.design$order != 1)])

bntihwsed.df <- cbind.data.frame(liste(bnti.hw.sed, entry = "bNTI"), order = "headwater", habitat = "sediment")
bntihwwat.df <- cbind.data.frame(liste(bnti.hw.water, entry = "bNTI"), order = "headwater", habitat = "plankton")
bntidssed.df <- cbind.data.frame(liste(bnti.ds.sed, entry = "bNTI"), order = "downstream", habitat = "sediment")
bntidswat.df <- cbind.data.frame(liste(bnti.ds.water, entry = "bNTI"), order = "downstream", habitat = "plankton")

bnti.df <- rbind.data.frame(bntihwsed.df, bntihwwat.df, bntidssed.df, bntidswat.df)

rcbray.hw.sed <- as.dist(
  as.matrix(RC.bray.dist)[which(design$order == 1 & design$habitat == "sediment"), 
                          which(design$order == 1 & design$habitat == "sediment")])

rcbray.ds.sed <- as.dist(
  as.matrix(RC.bray.dist)[which(design$order != 1 & design$habitat == "sediment"), 
                           which(design$order != 1 & design$habitat == "sediment")])
rcbray.hw.water <- as.dist(
  as.matrix(RC.bray.dist)[which(design$order == 1 & design$habitat == "water"),
                               which(design$order == 1 & design$habitat == "water")])
rcbray.ds.water <- as.dist(
  as.matrix(RC.bray.dist)[which(design$order != 1 & design$habitat == "water"),
                               which(design$order != 1 & design$habitat == "water")])

rcbhwsed.df <- cbind.data.frame(liste(rcbray.hw.sed, entry = "RCbray"), order = "headwater", habitat = "sediment")
rcbhwwat.df <- cbind.data.frame(liste(rcbray.hw.water, entry = "RCbray"), order = "headwater", habitat = "plankton")
rcbdssed.df <- cbind.data.frame(liste(rcbray.ds.sed, entry = "RCbray"), order = "downstream", habitat = "sediment")
rcbdswat.df <- cbind.data.frame(liste(rcbray.ds.water, entry = "RCbray"), order = "downstream", habitat = "plankton")

rcb.df <- rbind.data.frame(rcbhwsed.df, rcbhwwat.df, rcbdssed.df, rcbdswat.df)

comm.assemb.df <- full_join(rcb.df, bnti.df)
assemb.mech <- vector(length = nrow(comm.assemb.df))
for(row.i in 1:nrow(comm.assemb.df)){
  if(comm.assemb.df[row.i,"bNTI"] >= -2 && comm.assemb.df[row.i,"RCbray"] <= -0.95){
    assemb.mech[row.i] <- str_wrap("Mass Effects", width = 12)
  }
  if(comm.assemb.df[row.i,"bNTI"] <= 2 & comm.assemb.df[row.i,"RCbray"] >= 0.95){
    assemb.mech[row.i] <- str_wrap("Dispersal Limitation", width = 12)
  }
  if(comm.assemb.df[row.i,"bNTI"] < 2 & comm.assemb.df[row.i,"bNTI"] > -2 & 
     comm.assemb.df[row.i,"RCbray"] < 0.95 & comm.assemb.df[row.i,"RCbray"] > -0.95){
    assemb.mech[row.i] <- str_wrap("Undominated", width = 12)
  }
  if(comm.assemb.df[row.i,"RCbray"] <= .95 & comm.assemb.df[row.i,"bNTI"] <= -2){
    assemb.mech[row.i] <- str_wrap("Selection (Convergent)", width = 12)
  }
  if(comm.assemb.df[row.i,"RCbray"] >= -.95 & comm.assemb.df[row.i,"bNTI"] >= 2){
    assemb.mech[row.i] <- str_wrap("Selection (Divergent)", width = 12)
  }
  
}
comm.assemb.df$mechanism <- assemb.mech

assembly.table <- comm.assemb.df %>% count(mechanism, habitat, order)

assembly.table <- rbind(
  rbind(assembly.table, c("Selection\n(Convergent)", "sediment", "headwater", 0)),
  c("Selection\n(Convergent)", "sediment", "downstream", 0))
assembly.counts <- ggplot(assembly.table, aes(x = mechanism, y = n, fill = habitat)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("wheat", "skyblue")) +
  theme_minimal() + 
  labs(y = "Count", x = "Community Assembly Mechanism") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)) +
  facet_grid(facets = order ~ ., scales = "free_y")
assembly.counts

pdf("~/Desktop/assemb.pdf", height = 7, width = 9, bg = "white")
ggplot(comm.assemb.df, aes(x = mechanism, fill = habitat)) + 
  geom_bar() + 
  scale_fill_manual(values = c("wheat", "skyblue")) +
  theme_minimal() + 
  labs(y = "Count", x = "\nCommunity Assembly Mechanism") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14), 
        text = element_text(size = 16)) +
  facet_grid(facets = order ~ ., scales = "free_y")
dev.off()

ggplot(comm.assemb.df, aes(y = RCbray, x = bNTI, col = habitat)) + 
  geom_point() + 
  facet_grid(~ order)

water.scores <- data.frame(water.pcoa$pcoa$points, order = water.design$order)
sed.scores <- data.frame(sed.pcoa$pcoa$points, order = sed.design$order)

require(cowplot)
ggplot(data = water.scores, aes(x = PCoA1, y = PCoA2, col = order)) +
  geom_point(size = 2)+
  scale_color_gradientn(colors = viridis::viridis(5))

ggplot(data = sed.scores, aes(x = PCoA1, y = PCoA2, col = order)) +
  geom_point(size = 2) + 
  scale_color_gradientn(colors = viridis::inferno(5))

water.betadisp <- betadisper(water.dists$otus, group = water.design$order, type = "centroid")
sed.betadisp <- betadisper(sediment.dists$otus, group = sed.design$order, type = "centroid")
permutest(water.betadisp)
pdf(file = "figures/betadisper_water.pdf", bg = "white")
par(oma = c(1,1,1,1), pty = "s")
plot(water.betadisp, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", sub = "")
add.axes(xlab = "PCoA 1", ylab = "PCoA 2", main = "Plankton")
dev.off()

permutest(sed.betadisp)
pdf(file = "figures/betadisper_sediment.pdf", bg = "white")
par(oma = c(1,1,1,1), pty = "s")
plot(sed.betadisp, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", sub = "")
add.axes(xlab = "PCoA 1", ylab = "PCoA 2", main = "Sediments")
dev.off()



hja.betadisp <- betadisper(hja.dists$otus, group = design$flow, type = "centroid")
permutest(hja.betadisp)
 