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

