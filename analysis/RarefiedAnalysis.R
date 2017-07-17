# Initial setup
rm(list=ls())
#setwd("~/GitHub/HJA-streams/")
opar <- par()


# Check for and install required packages
package.list <- c('vegan', 'png', 'simba', 'grid', 
                  'vegetarian', 'pander', 'SoDA', 'fossil',
                  'ggplot2')
# 'sp', 'vegetarian', 
# 'SoDA', 'geoR',

for (package in package.list) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# Load packages and other tools
source("./analysis/MothurTools.R")

se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


## Import Shared, Design, and Environment Files

# Define Inputs
# Design = general design file for experiment
# shared = OTU table from mothur with sequence similarity clustering
# Taxonomy = Taxonomic information for each OTU
design <- "./data/design.txt"
shared <- "./data/hja_streams.final.shared"
taxon  <- "./data/hja_streams.final.0.03.taxonomy"
env    <- "./data/hja_env.csv"

# Import Design
design.total <- read.delim(design, header=T, row.names=1)

# Import Shared Files
OTUs <- read.otu(shared = shared, cutoff = "0.03") # 97% Similarity

# Import Taxonomy
OTU.tax <- read.tax(taxonomy = taxon, format = "rdp")

# Import Env
env.total <- read.csv(env, header=T)

### Data Transformations
# Remove OTUs with less than two occurances across all sites
OTUs <- OTUs[, which(colSums(OTUs) >= 2)]

# Sequencing an Good's Coverage
# Sequencing Coverage
coverage <- rowSums(OTUs)

# Good's Coverage
goods <- function(x = ""){
  1 - (rowSums(x == 1) / rowSums(x))
}
goods.c <- goods(OTUs)

# Remove Low Coverage Samples
#lows <- which(coverage < 7000)
OTUs <- OTUs[-which(goods.c < 0.95), ]
design <- design.total[-which(goods.c < 0.95), ]
env <- env.total[-which(goods.c < 0.95), ]
summary(goods(OTUs))



ComLoop <- vector(mode = "list", length =  10)
for (i in seq_along(ComLoop)) {
  ComLoop[[i]] <- rrarefy(OTUs, sample = min(rowSums(OTUs)))
}

lapply(X = ComLoop, FUN = diversity, index = "shannon")
OTUs.rarefied.hel <- lapply(X = ComLoop, FUN = decostand, method = "hel")
OTUs.bray.d <- lapply(X = OTUs.rarefied.hel, FUN = vegdist, method = "bray")
OTUs.bray.mean <- lapply(OTUs.bray.d, mean)


OTUs.log <- decostand(OTUs, method = "log")
OTUs.hel <- decostand(OTUs, method = "hel")
OTUs.rel <- decostand(OTUs, method = "total")





hja.db <- vegdist(OTUs.hel, method = "bray")
hja.d.sorensen <- vegdist(OTUs.hel, method = "bray", binary = T)
hja.d.jaccard <- vegdist(OTUs.hel, method = "jaccard")
hja.pcoa <- cmdscale(hja.db, eig=TRUE)
hja.pcoa.sorensen <- cmdscale(hja.d.sorensen, eig = TRUE)
hja.pcoa.jaccard <- cmdscale(hja.d.jaccard, eig = TRUE)
var1 <- round(hja.pcoa$eig[1] / sum(hja.pcoa$eig),3) * 100
var2 <- round(hja.pcoa$eig[2] / sum(hja.pcoa$eig),3) * 100
var3 <- round(hja.pcoa$eig[3] / sum(hja.pcoa$eig),3) * 100
var1.sor <- round(hja.pcoa.sorensen$eig[1] / sum(hja.pcoa.sorensen$eig),3) * 100
var2.sor <- round(hja.pcoa.sorensen$eig[2] / sum(hja.pcoa.sorensen$eig),3) * 100
var3.sor <- round(hja.pcoa.sorensen$eig[3] / sum(hja.pcoa.sorensen$eig),3) * 100
var1.jac <- round(hja.pcoa.jaccard$eig[1] / sum(hja.pcoa.jaccard$eig),3) * 100
var2.jac <- round(hja.pcoa.jaccard$eig[2] / sum(hja.pcoa.jaccard$eig),3) * 100
var3.jac <- round(hja.pcoa.jaccard$eig[3] / sum(hja.pcoa.jaccard$eig),3) * 100

den.d <- as.dist(den.dists)
plot(hja.db ~ den.d)

env.d <- cluster::daisy(x = env.mat, metric = "gower")


plot(hja.db ~ den.d)
