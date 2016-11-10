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
EMS.water <- Coherence(comm = OTUs.PA[which(design$habitat=="water"),
                                          which(colnames(OTUsREL) %in% top.taxa)], 
                       method = "tswap", sims = 100, verbose = T)
EMS.sed <- Metacommunity(comm = OTUs.PA[which(design$habitat=="sediment"),
                                        which(colnames(OTUsREL) %in% top.taxa)], method = "tswap", verbose = TRUE)