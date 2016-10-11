require(picante)

hja.tree <- read.tree(file = "./data/hja_streams.tree")

matched.phylo <- match.phylo.comm(hja.tree, OTUs)

cophenetic(matched.phylo$phy)

pd(samp = matched.phylo$comm, tree = matched.phylo$phy, include.root = F)
