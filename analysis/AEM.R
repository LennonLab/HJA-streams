# AEM 

site.by.edges.lc <- as.matrix(read.csv("./data/sites-by-edges_LC.csv", header=T))
row.names(site.by.edges.lc) <- site.by.edges.lc[,1]
site.by.edges.lc <- site.by.edges.lc[,-1]
site.by.edges.lc <- (site.by.edges.lc == 1) * 1

site.by.edges.w1 <- as.matrix(read.csv("./data/sites-by-edges_WS1.csv", header=T))
row.names(site.by.edges.w1) <- site.by.edges.w1[,1]
site.by.edges.w1 <- site.by.edges.w1[,-1]
site.by.edges.w1 <- (site.by.edges.w1 == 1) * 1

w1.dists <- c(1,
              1,
              255,
              10,
              4,
              6.46,
              201.47,
              196.598,
              3,
              4,
              100.3,
              116.85,
              802.3,
              1065,
              40.39,
              1)


w1.F.mat <- princomp(site.by.edges.w1)
summary(w1.F.mat)
plot(w1.F.mat)

w1.env <- princomp(env.mat[which(design$watershed == "WS01"),])
summary(w1.env)
plot(w1.env)

