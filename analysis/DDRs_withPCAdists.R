
# Principal Components Analysis on HJA Environment
hja.pca <- princomp(env.mat)
summary(hja.pca)
plot(hja.pca, type = "l")
biplot(hja.pca)
pc1 <- hja.pca$scores[,1]
pc2 <- hja.pca$scores[,2]
pc3 <- hja.pca$scores[,3]
pc4 <- hja.pca$scores[,4]
pc5 <- hja.pca$scores[,5]
cor(pc1, env.mat)
cor(pc2, env.mat)
cor(pc3, env.mat)
cor(pc4, env.mat)
cor(pc5, env.mat)

?dist
pc.dists <- sqrt((dist(pc1))^2 + (dist(pc2))^2 + (dist(pc3)^2))
pc.dists <- liste(pc.dists)[,3]
plot(liste(hja.db)[,3] ~ pc.dists)


water.pc1 <- pc1[which(design$habitat == "water")]
water.pc2 <- pc2[which(design$habitat == "water")]
water.pc3 <- pc3[which(design$habitat == "water")]
water.pc.dists <- sqrt((dist(water.pc1))^2 + (dist(water.pc2))^2 + (dist(water.pc3)^2))
water.pc.dists <- liste(water.pc.dists)[,3]
plot(water.struc.dist.ls ~ water.pc.dists)
summary(lm(water.struc.dist.ls ~ water.pc.dists))

sediment.pc1 <- pc1[which(design$habitat == "sediment")]
sediment.pc2 <- pc2[which(design$habitat == "sediment")]
sediment.pc3 <- pc3[which(design$habitat == "sediment")]
sediment.pc.dists <- sqrt((dist(sediment.pc1))^2 + (dist(sediment.pc2))^2 + (dist(sediment.pc3)^2))
sediment.pc.dists <- liste(sediment.pc.dists)[,3]
plot(sed.struc.dist.ls ~ sediment.pc.dists)
summary(lm(sed.struc.dist.ls ~ sediment.pc.dists))

headwater.pc1 <- pc1[which(design$order == 1)]
headwater.pc2 <- pc2[which(design$order == 1)]
headwater.pc3 <- pc2[which(design$order == 1)]
headwater.pc.dists <- sqrt((dist(headwater.pc1))^2 + (dist(headwater.pc2))^2 + (dist(headwater.pc3)^2))
headwater.pc.dists <- liste(headwater.pc.dists)[,3]
plot(headwater.dists$comm.struc ~ headwater.pc.dists)
summary(lm(headwater.dists$comm.struc ~ headwater.pc.dists))

downstream.pc1 <- pc1[which(design$order != 1)]
downstream.pc2 <- pc2[which(design$order != 1)]
downstream.pc3 <- pc3[which(design$order != 1)]
downstream.pc.dists <- sqrt((dist(downstream.pc1))^2 + (dist(downstream.pc2))^2 + (dist(downstream.pc3)^2))
downstream.pc.dists <- liste(downstream.pc.dists)[,3]
plot(mainstem.dists$comm.struc ~ downstream.pc.dists)
summary(lm(mainstem.dists$comm.struc ~ downstream.pc.dists))

