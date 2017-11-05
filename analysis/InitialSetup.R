# Initial setup
rm(list=ls())
#setwd("~/GitHub/HJA-streams/")
opar <- par()


# Check for and install required packages
package.list <- c('vegan', 'png', 'grid', 'simba',
                  'vegetarian', 'pander', 'SoDA', 'fossil',
                  'tidyverse', 'cluster', 'adespatial', 'spdep', 
                  'stringr', 'picante')
# 'sp', 'vegetarian', 
# 'SoDA', 'geoR',

for (package in package.list) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


source("./analysis/MothurTools.R")
source("analysis/HJA-Functions.R")

## Import Shared, Design, and Environment Files

# Define Inputs
# Design = general design file for experiment
# shared = OTU table from mothur with sequence similarity clustering
# Taxonomy = Taxonomic information for each OTU

# # Import Design
# design.total <- read.delim("./data/design.txt", header=T, row.names=1)
# 
# # Import Shared Files
# OTUs <- read.otu(shared = "./data/hja_streams.final.shared", cutoff = "0.03") # 97% Similarity
# 
# # Import Taxonomy
# OTU.tax <- read.tax(taxonomy = "./data/hja_streams.final.0.03.taxonomy", format = "rdp")
# 
# # Import Env
# env.total <- read.csv("./data/hja_env.csv", header=T)
# 
# ### Data Transformations
# 
# # Sequencing an Good's Coverage
# # Sequencing Coverage
# coverage <- rowSums(OTUs)
# 
# # Good's Coverage
# goods <- function(x = ""){
#   1 - (rowSums(x == 1) / rowSums(x))
# }
# goods.c <- goods(OTUs)
# 
# # Remove Low Coverage Samples
# cutoff <- 10000
# lows <- which(coverage < cutoff)
# OTUs <- OTUs[-which(coverage < cutoff), ]
# design <- design.total[-which(coverage < cutoff), ]
# env <- env.total[-which(coverage < cutoff), ]
# 
# # Remove OTUs with less than two occurances across all sites
# OTUs <- OTUs[, which(colSums(OTUs) >= 2)]
# 
# 
# # Write and read data files
# saveRDS(OTUs, file = "./data/SiteBySpecies.rda")
# saveRDS(env, file = "./data/SiteByEnv.rda")
# saveRDS(OTU.tax, file = "./data/Taxonomy.rda")
# saveRDS(design, file = "./data/SiteDesign.rda")

OTUs <- readRDS(file = "./data/SiteBySpecies.rda")
env <- readRDS(file = "./data/SiteByEnv.rda")
OTU.tax <- readRDS(file = "./data/Taxonomy.rda")
design <- readRDS(file = "./data/SiteDesign.rda")
hja.unifrac.dist <- readRDS(file = "data/UnifracDists.rda")
den.dists <- as.dist(readRDS(file = "data/DendriticDists.rda"))
design$upstreamdist <- as.matrix(den.dists)[1,]


# Remove orthogonal vectors
env.mat <- as.matrix(env[10:20])
env.mat["83",9] <- 150 # These were the highest samples, overflow
#env.mat[52,5] <- 150
for(i in 1:nrow(env.mat)){
  if(env.mat[i, 9] < 0){
    env.mat[i, 9] <- 0.0001 # these were below detection of the machine
  }
  if(env.mat[i, 10] < 0){
    env.mat[i, 10] <- 0.0001
  }
}
habitat.dummy <- simba::mad(as.factor(env$habitat))
env.mat <- cbind(habitat.dummy, env.mat)
env.mat <- scale(env.mat)
env.mat <- env.mat[,c("sediment", "elevation", "temperature", "conductivity",
           "ph", "TN", "TP", "DOC")]

# Rarefy communities
#OTUs <- rrarefy(OTUs, sample = min(rowSums(OTUs)))
#OTUs <- OTUs[,-which(colSums(OTUs) == 0)]
#saveRDS(OTUs, file = "data/SiteBySpeciesRarefied.rda")
OTUs <- readRDS("data/SiteBySpeciesRarefied.rda")

# Transformations and Standardizations
OTUsREL <- decostand(OTUs, method = "total")
OTUs.PA <- decostand(OTUs, method = "pa")
OTUsREL.log <- decostand(OTUs, method = "log")
OTUsREL.hel <- decostand(OTUs, method = "hellinger")

# Read in Distances
# Geo distance Matrix
xy <- cbind(jitter(env$longitude, amount = .0001),
            jitter(env$latitude, amount = .0001))
#geo.dists <- geoXY(env$latitude, env$longitude)
#xy <- project(xy, "+proj=utm +zone=10 +ellps=WGS84")
#dist.mat <- as.matrix(dist(xy, method = "euclidean"))
dist.mat <- fossil::earth.dist(xy) * 1000

# # Read in phylodist
hja.unifrac.raw <- read.delim(file = "./data/hja_streams.tree1.weighted.phylip.dist", header = F, skip = 1, row.names = 1)
colnames(hja.unifrac.raw) <- as.vector(lapply(strsplit(rownames(hja.unifrac.raw)," "), function(x) x[1]))
rownames(hja.unifrac.raw) <- colnames(hja.unifrac.raw)
hja.unifrac <- hja.unifrac.raw[which(rownames(hja.unifrac.raw) %in% rownames(OTUs)),
                               which(rownames(hja.unifrac.raw) %in% rownames(OTUs))]
hja.unifrac.dist <- as.dist(hja.unifrac)

RC.bray.dist <- readRDS(file = "./data/RCbraydist.csv")
bNTI.water.dist <- readRDS(file = "data/bNTIwater.rda")
bNTI.sed.dist <- readRDS(file = "data/bNTIsed.rda")
