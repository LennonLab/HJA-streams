# AEM 

design.sed <- subset(design, habitat == "sediment")
design.water <- subset(design, habitat == "water")


E.sed.connection <- matrix(NA, nrow = nrow(design.sed), 
                           ncol = nrow(design.sed))
rownames(E.sed.connection) <- rownames(design.sed)
colnames(E.sed.connection) <- rownames(design.sed)

#write.csv(E.sed.connection, file = "data/AEMsed.csv", row.names = T, col.names = T)

E.sed.connection <- read.csv("data/AEMsed.csv", header = T, row.names = 1)
E.sed.connection <- as.matrix(E.sed.connection)

sed <- OTUsREL[which(design$habitat == "sediment"),]

cmdscale(dist(E.sed.connection), eig = T)
F.eigen

sed.aem <- aem(binary.mat = E.sed.connection)

water.dist <- vegdist(water)



sitebyedge <- read.csv(file = "data/sites-by-edges.csv", header = T, row.names = 1)


# Construct side-by-edge matrix for sediments
sed.E <- matrix(0, nrow = 24, ncol = 31)
rownames(sed.E) <- rownames(design[which(design$habitat == "sediment"),])
sed.E[1,c(1:17)] <- 1
sed.E[2,c(2:17)] <- 1
sed.E[3,c(5,6,7,8,9,10,12,13,14,15)] <- 1
sed.E[4,17] <- 1
sed.E[5,13] <- 1
sed.E[6,12] <- 1
sed.E[7,16] <- 1
sed.E[8,c(8,9,10,14,15)] <- 1
sed.E[9,11] <- 1
sed.E[10,c(9,10,14,15)] <- 1
sed.E[11,15] <- 1
sed.E[12,14] <- 1
sed.E[13,c(18:31)] <- 1
sed.E[14,c(19:31)] <- 1
sed.E[15,c(22:31)] <- 1
sed.E[16,c(23:31)] <- 1
sed.E[17,25] <- 1
sed.E[18,c(26:31)] <- 1
sed.E[19,c(27:31)] <- 1
sed.E[20,c(29,30)] <- 1
sed.E[21,30] <- 1
sed.E[22,31] <- 1
sed.E[23,c(20:31)] <- 1
sed.E[24,c(21:31)] <- 1

sed.aem <- aem(binary.mat = sed.E)
plot(sed.aem$values/sum(sed.aem$values))
sed.db <- vegdist(OTUsREL[which(design$habitat == "sediment"),])
mod0 <- dbrda(sed.db ~ 1, data = as.data.frame(sed.aem$vectors))
mod1 <- dbrda(sed.db ~ ., data = as.data.frame(sed.aem$vectors))
envfit(mod1, as.data.frame(sed.aem$vectors))
anova.cca(mod1, by = "axis")
sed.aem.mod <- ordiR2step(mod0, mod1)

sed.aem.2 <- cmdscale(dist(sed.E), eig = T, k = 23)
plot(sed.aem.2$eig/sum(sed.aem.2$eig))
sed.aem.eigvec <- as.data.frame(sed.aem.2$points)
mod1 <- dbrda(sed.db ~ ., data = sed.aem.eigvec)
envfit(mod1, sed.aem.eigvec)
sed.eigs <- subset(sed.aem.eigvec, select = c(V2, V3, V7, V11, V12, V13, V19, V21, V23))
sed.spatial0 <- dbrda(sed.db ~ 1, sed.eigs)
sed.spatial1 <- dbrda(sed.db ~ ., sed.eigs[,c("V2", "V7")])
sed.spatial <- ordiR2step(sed.spatial0, sed.spatial1)
plot(sed.spatial)
envfit(sed.spatial1, sed.eigs)
