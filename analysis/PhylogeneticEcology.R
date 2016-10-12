require(picante)

hja.tree <- read.tree(file = "./data/hja_streams.tree")

matched.phylo <- match.phylo.comm(hja.tree, OTUs)

