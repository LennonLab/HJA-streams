require(betapart)
require(vegan)
OTUs.PA <- decostand(OTUs, method = "pa")
hja.beta.pair <- beta.pair(OTUs.PA)
ordiplot(cmdscale(hja.beta.pair$beta.sne))


hja.beta.samp.water <- beta.sample(OTUs.PA[which(design$habitat=="water"),],sites=2,samples = 999)
hja.beta.samp.water$mean.values

hja.beta.samp.sed <- beta.sample(OTUs.PA[which(design$habitat!="water"),],sites=2,samples = 999)
hja.beta.samp.sed$mean.values




#### EMS Framework
require(metacom)
EMS.water <- Metacommunity(comm = OTUs.PA[which(design$habitat=="water"),], method = "tswap", verbose = TRUE)
EMS.sed <- Metacommunity(comm = OTUs.PA[which(design$habitat=="sediment"),], method = "tswap", verbose = TRUE)