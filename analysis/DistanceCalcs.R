# source("./analysis/InitialSetup.R")

# Calculate Dendritic Distances
require('igraph')
adj.mat <- as.matrix(read.csv("./data/undirected-matrix.csv", header=T))
row.names(adj.mat) <- adj.mat[,1]
adj.mat <- adj.mat[,-1]
adj.mat <- (adj.mat == 1) * 1

stream.network <- graph_from_adjacency_matrix(adjmatrix = adj.mat)

png(filename = "./figures/stream-network.png", 
    height = 4800, width = 4800, res = 2*96)
plot.igraph(stream.network)
dev.off()

# Create Dist Matrix
xy.total <- cbind(env.total$longitude, env.total$latitude)
xy.total <- project(xy.total, "+proj=utm +zone=10 +ellps=WGS84")
euc.dist.mat <- as.matrix(dist(xy.total, method = "euclidean"))
rownames(euc.dist.mat) <- env.total$sample
colnames(euc.dist.mat) <- env.total$sample

# Function to find paths along the network.
dend.dist <- function(graph = "", dist.mat = ""){
  
  dend.dist.mat <- matrix(data = NA, nrow = length(env.total$sample), 
                          ncol = length(env.total$sample))
  rownames(dend.dist.mat) <- env.total$sample
  colnames(dend.dist.mat) <- env.total$sample
  
  for(row in rownames(dend.dist.mat)){
    for(col in colnames(dend.dist.mat)){
      path.dist <- 0
      path <- shortest_paths(graph = graph, from = row, to = col)$vpath[[1]]
      
      if(length(path) > 1){
        for(i in 1:(length(path) - 1)){
          path.dist <- path.dist + dist.mat[path[[i]]$name, path[[i+1]]$name]
        }
      }
      
      dend.dist.mat[row,col] <- path.dist
    }
  }
  
  return(dend.dist.mat)
}

dend.dist.mat <- dend.dist(graph = stream.network, dist.mat = euc.dist.mat)

# Make distance lists
xy <- cbind(env$longitude, env$latitude)
xy <- project(xy, "+proj=utm +zone=10 +ellps=WGS84")

dist.mat <- as.matrix(dist(xy, method = "euclidean"))
dist.mat[!lower.tri(dist.mat)] <- NA 
den.dists <- dend.dist.mat[which(rownames(dend.dist.mat) %in% rownames(design)),
                       which(rownames(dend.dist.mat) %in% rownames(design))]
den.dists[!lower.tri(den.dists)] <- NA 
env.2 <- env.mat[,2:6]+2.685105 # Make all values positive
