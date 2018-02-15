source("analysis/InitialSetup.R")
source("analysis/DistanceCalcs.R")
source("analysis/DDRs.R")

require(betapart)
require(vegan)

OTUs.PA <- decostand(OTUs, method = "pa")
hja.beta.pair <- beta.pair(OTUs.PA[which(design$habitat == "water"),], index.family = "jac")
nest <- as.matrix(hja.beta.pair$beta.jne)[1,-1]
turn <- as.matrix(hja.beta.pair$beta.jtu)[1,-1]
usdis <- design$upstreamdist[which(design$habitat == "water")][-1]
beta <- as.matrix(hja.beta.pair$beta.jac)[1,-1]
plot(nest~usdis)
plot(turn~usdis)
plot(beta~usdis)
water.env <- env.mat[which(design$habitat == "water"),]
water.env <- as.matrix(dist(water.env))[1,-1]
plot(nest~water.env)
plot(turn~water.env)
plot(beta~water.env)

hja.beta.pair <- beta.pair(OTUs.PA[which(design$habitat != "water"),], index.family = "jac")
nest <- as.matrix(hja.beta.pair$beta.jne)[1,-1]
turn <- as.matrix(hja.beta.pair$beta.jtu)[1,-1]
beta <- as.matrix(hja.beta.pair$beta.jac)[1,-1]
usdis <- design$upstreamdist[which(design$habitat != "water")][-1]
plot(nest~usdis)
plot(turn~usdis)
plot(beta~usdis)
water.env <- env.mat[which(design$habitat != "water"),]
water.env <- as.matrix(dist(water.env))[1,-1]
plot(nest~water.env)
plot(turn~water.env)
plot(beta~water.env)

hja.beta.pair <- beta.pair(OTUs.PA[which(design$order == 1),], index.family = "jac")
nest <- as.matrix(hja.beta.pair$beta.jne)[1,-1]
turn <- as.matrix(hja.beta.pair$beta.jtu)[1,-1]
beta <- as.matrix(hja.beta.pair$beta.jac)[1,-1]
usdis <- design$upstreamdist[which(design$order == 1)][-1]
plot(nest~usdis)
plot(turn~usdis)
plot(beta~usdis)
water.env <- env.mat[which(design$order == 1),]
water.env <- as.matrix(dist(water.env))[1,-1]
plot(nest~water.env)
plot(turn~water.env)
plot(beta~water.env)

hja.beta.pair <- beta.pair(OTUs.PA[which(design$order != 1),], index.family = "jac")
nest <- as.matrix(hja.beta.pair$beta.jne)[1,-1]
turn <- as.matrix(hja.beta.pair$beta.jtu)[1,-1]
beta <- as.matrix(hja.beta.pair$beta.jac)[1,-1]
usdis <- design$upstreamdist[which(design$order != 1)][-1]
plot(nest~usdis)
plot(turn~usdis)
plot(beta~usdis)
water.env <- env.mat[which(design$order != 1),]
water.env <- as.matrix(dist(water.env))[1,-1]
plot(nest~water.env)
plot(turn~water.env)
plot(beta~water.env)

hja.beta.pair$beta.sim
hja.beta.pair$beta.sne
hja.beta.pair$beta.sor

# Distance Decay with Env
plot(hja.beta.pair$beta.sor ~ hja.dists$env)
plot(hja.beta.pair$beta.sim ~ hja.dists$env) # turnover
plot(hja.beta.pair$beta.sne ~ hja.dists$env) # nestedness

water.betapart <- beta.pair(decostand(water, method = "pa"))
sed.betapart <- beta.pair(decostand(sediment, method = "pa"))

plot(water.betapart$beta.sor ~ water.dists$env)
plot(water.betapart$beta.sim ~ water.dists$env) # turnover
plot(water.betapart$beta.sne ~ water.dists$env) # nestedness

plot(sed.betapart$beta.sor ~ sed.dists$env)
plot(sed.betapart$beta.sim ~ sed.dists$env) # turnover
plot(sed.betapart$beta.sne ~ sed.dists$env) # nestedness


downstream.betapart <- beta.pair(decostand(OTUs.PA[which(design$order>1),], method = "pa"))
headwater.betapart <- beta.pair(decostand(OTUs.PA[which(design$order==1),], method = "pa"))

plot(headwater.betapart$beta.sor ~ headwater.dists$env)
plot(headwater.betapart$beta.sim ~ headwater.dists$env) # turnover
plot(headwater.betapart$beta.sne ~ headwater.dists$env) # nestedness

plot(downstream.betapart$beta.sor ~ downstream.dists$env)
plot(downstream.betapart$beta.sim ~ downstream.dists$env) # turnover
plot(downstream.betapart$beta.sne ~ downstream.dists$env) # nestedness


downstream.sed.betapart <- beta.pair(decostand(OTUs.PA[which(design$order>1 & design$habitat == "sediment"),], method = "pa"))
plot(downstream.sed.betapart$beta.sor ~ downstream.sed.dists$env)
plot(downstream.sed.betapart$beta.sim ~ downstream.sed.dists$env) # turnover
plot(downstream.sed.betapart$beta.sne ~ downstream.sed.dists$env) # nestedness

downstream.water.betapart <- beta.pair(decostand(OTUs.PA[which(design$order>1 & design$habitat == "water"),], method = "pa"))
plot(downstream.water.betapart$beta.sor ~ downstream.water.dists$env)
plot(downstream.water.betapart$beta.sim ~ downstream.water.dists$env) # turnover
plot(downstream.water.betapart$beta.sne ~ downstream.water.dists$env) # nestedness


# Distance Decay with Space
plot(hja.beta.pair$beta.sor ~ hja.dists$den)
plot(hja.beta.pair$beta.sim ~ hja.dists$den) # turnover
plot(hja.beta.pair$beta.sne ~ hja.dists$den) # nestedness

plot(water.betapart$beta.sor ~ water.dists$den)
plot(water.betapart$beta.sim ~ water.dists$den) # turnover
plot(water.betapart$beta.sne ~ water.dists$den) # nestedness

plot(sed.betapart$beta.sor ~ sed.dists$den)
plot(sed.betapart$beta.sim ~ sed.dists$den) # turnover
plot(sed.betapart$beta.sne ~ sed.dists$den) # nestedness

plot(headwater.betapart$beta.sor ~ headwater.dists$den)
plot(headwater.betapart$beta.sim ~ headwater.dists$den) # turnover
plot(headwater.betapart$beta.sne ~ headwater.dists$den) # nestedness

plot(downstream.betapart$beta.sor ~ downstream.dists$den)
plot(downstream.betapart$beta.sim ~ downstream.dists$den) # turnover
plot(downstream.betapart$beta.sne ~ downstream.dists$den) # nestedness

plot(downstream.sed.betapart$beta.sor ~ downstream.sed.dists$den)
plot(downstream.sed.betapart$beta.sim ~ downstream.sed.dists$den) # turnover
plot(downstream.sed.betapart$beta.sne ~ downstream.sed.dists$den) # nestedness

plot(downstream.water.betapart$beta.sor ~ downstream.water.dists$den)
plot(downstream.water.betapart$beta.sim ~ downstream.water.dists$den) # turnover
plot(downstream.water.betapart$beta.sne ~ downstream.water.dists$den) # nestedness



dt.hja.env <- residuals(lm(hja.dists$env ~ hja.dists$den))

# Distance Decay with Env
plot(hja.beta.pair$beta.sor ~ dt.hja.env)
plot(hja.beta.pair$beta.sim ~ dt.hja.env) # turnover
plot(hja.beta.pair$beta.sne ~ dt.hja.env) # nestedness

plot(water.betapart$beta.sor ~ residuals(lm(water.dists$env ~ water.dists$den)))
plot(water.betapart$beta.sim ~ residuals(lm(water.dists$env ~ water.dists$den))) # turnover
plot(water.betapart$beta.sne ~ residuals(lm(water.dists$env ~ water.dists$den))) # nestedness

plot(sed.betapart$beta.sor ~ residuals(lm(sed.dists$env ~ sed.dists$den)))
plot(sed.betapart$beta.sim ~ residuals(lm(sed.dists$env ~ sed.dists$den))) # turnover
plot(sed.betapart$beta.sne ~ residuals(lm(sed.dists$env ~ sed.dists$den))) # nestedness


downstream.betapart <- beta.pair(decostand(OTUs.PA[which(design$order>1),], method = "pa"))
headwater.betapart <- beta.pair(decostand(OTUs.PA[which(design$order==1),], method = "pa"))

plot(headwater.betapart$beta.sor ~ residuals(lm(headwater.dists$env ~ headwater.dists$den)))
plot(headwater.betapart$beta.sim ~ residuals(lm(headwater.dists$env ~ headwater.dists$den))) # turnover
plot(headwater.betapart$beta.sne ~ residuals(lm(headwater.dists$env ~ headwater.dists$den))) # nestedness

plot(downstream.betapart$beta.sor ~ residuals(lm(downstream.dists$env ~ downstream.dists$den)))
plot(downstream.betapart$beta.sim ~ residuals(lm(downstream.dists$env ~ downstream.dists$den))) # turnover
plot(downstream.betapart$beta.sne ~ residuals(lm(downstream.dists$env ~ downstream.dists$den))) # nestedness


downstream.sed.betapart <- beta.pair(decostand(OTUs.PA[which(design$order>1 & design$habitat == "sediment"),], method = "pa"))
plot(downstream.sed.betapart$beta.sor ~ downstream.sed.dists$env)
plot(downstream.sed.betapart$beta.sim ~ downstream.sed.dists$env) # turnover
plot(downstream.sed.betapart$beta.sne ~ downstream.sed.dists$env) # nestedness

downstream.water.betapart <- beta.pair(decostand(OTUs.PA[which(design$order>1 & design$habitat == "water"),], method = "pa"))
plot(downstream.water.betapart$beta.sor ~ downstream.water.dists$env)
plot(downstream.water.betapart$beta.sim ~ downstream.water.dists$env) # turnover
plot(downstream.water.betapart$beta.sne ~ downstream.water.dists$env) # nestedness

