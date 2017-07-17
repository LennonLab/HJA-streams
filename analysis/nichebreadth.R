# Calculate niche breadth
source("analysis/InitialSetup.R")

library(tidyverse)

OTUsPA <- decostand(OTUs, method = "pa")
OTUsREL <- decostand(OTUs, method = "total")


nb.calc <- function(species = ""){
  Pij2 = species^2
  B = 1/sum(Pij2)
  return(B)
}

x <- rep(0.1, 10)
y <- c(rep(0,9), 1)


OTUsRELbySpec <- decostand(OTUs, method = "total", MARGIN = 2)
OTUsRELbySpec <- OTUsRELbySpec[,colSums(OTUsRELbySpec) > 0]

nb <- apply(X = OTUsRELbySpec, FUN = nb.calc, MARGIN = 2)
hist(nb)
site1 <- OTUsREL[5,]


water <- OTUsREL[which(env$habitat == "water"),which(colSums(OTUsREL)>0)]
sediment <- OTUsREL[which(env$habitat == "sediment"),which(colSums(OTUsREL)>0)]

water.df <- cbind.data.frame(nb, habitat = "water", t(water))
sed.df <- cbind.data.frame(nb, habitat = "sediment", t(sediment))

water.df.long <- water.df %>% gather(site, relabund, -nb, -habitat)
sed.df.long <- sed.df %>% gather(site, relabund, -nb, -habitat)

total.df.long <- rbind.data.frame(water.df.long, sed.df.long)

# plot(nb, waterabund, ylab = "Median Relative Abundance", xlab = "Niche Breadth")
# plot(nb, sedabund, ylab = "Median Relative Abundance", xlab = "Niche Breadth")
# 
# df <- rbind.data.frame(
#   cbind.data.frame(nb = nb, abund = waterabund, habitat = "water"),
#   cbind.data.frame(nb = nb, abund = sedabund, habitat = "sediment"))

nbplot <- ggplot(data = total.df.long, mapping = aes(y = relabund, x = nb)) +
  geom_point() +
  facet_grid(~habitat)+
  labs(y = "Relative Abundance", x = "Niche Breadth")+
  theme_bw()
nbplot

ls2 <- subset(water.df.long, site == "LC_02_W")
nbplot <- ggplot(data = ls2, mapping = aes(y = relabund, x = nb)) +
  geom_smooth() + 
  labs(y = "Relative Abundance", x = "Niche Breadth")
nbplot
