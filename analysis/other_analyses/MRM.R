require(ecodist)

den.dist <- as.dist(den.dists)

mant.corr <- ecodist::mgram(hja.db, den.dist)
plot(mant.corr)

pmg <- ecodist::pmgram(hja.db, den.dist, vegdist(env.mat, method = "gower"),
                       resids = T)
pmg$resids
plot(pmg)
plot(mantel(pmg$resids))


hja.mrm <- ecodist::MRM(formula = hja.db ~ den.dist, vegdist(env.mat, method = "gower"))


# Role of env after controlling for space
den.resids <- residuals(ecodist::pmgram(hja.db, den.dist, resids = T))
env.resids <- residuals(ecodist::pmgram(vegdist(env.mat, method = "gower"), den.dist, resids = T))
ecodist::mantel(den.resids ~ env.resids)


# Role of space after controlling for env
env.resids <- residuals(ecodist::pmgram(hja.db, vegdist(env.mat, method = "gower"), resids = T))
den.resids <- residuals(ecodist::pmgram(den.dist, vegdist(env.mat, method = "gower"), resids = T))
ecodist::mantel(env.resids ~ den.resids)
