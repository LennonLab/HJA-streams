# Initial setup
rm(list=ls())
#setwd("~/GitHub/HJA-streams/")
opar <- par()


# Check for and install required packages
package.list <- c('vegetarian', 'vegan', 'png', 'sp', 'rgdal',
                  'SoDA', 'grid', 'simba', 'geoR', 'raster')
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
OTUs <- read.otu(shared = shared, cutoff = "0.03")         # 97% Similarity

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
  1 - (sum(x == 1) / rowSums(x))
}
goods.c <- goods(OTUs)

# Remove Low Coverage Samples
lows <- which(coverage < 7000)
OTUs <- OTUs[-which(coverage < 7000), ]
design <- design.total[-which(coverage < 7000), ]
env <- env.total[-which(coverage < 7000), ]

# Remove orthogonal vectors
env <- env[c(1:11, 13, 16, 18, 19)]
env.mat <- as.matrix(env[10:15])
env.mat[51,5] <- 150
env.mat[52,5] <- 150
for(i in 1:nrow(env.mat)){
  if(env.mat[i, 5] < 0){
    env.mat[i, 5] <- 0.0001
  }
  if(env.mat[i, 6] < 0){
    env.mat[i, 6] <- 0.0001
  }
}
env.mat <- scale(env.mat)
env.pca <- princomp(env.mat, scores = T)

# Distance Matrix
xy <- cbind(env$longitude, env$latitude)
geo.dists <- geoXY(env$latitude, env$longitude)
xy <- project(xy, "+proj=utm +zone=10 +ellps=WGS84")
dist.mat <- as.matrix(dist(xy, method = "euclidean"))

# Make Relative Abundence Matrices
OTUsREL <- decostand(OTUs, method = "total")

# Transform Relative Abundances
OTUsREL.log <- decostand(OTUs, method = "log")
