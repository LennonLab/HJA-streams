# source("analysis/InitialSetup.R")
# source("analysis/RaupCrickBC.R")
# source("analysis/DDRs.R")
require("cowplot")
require("progress")

# Calculate abundance-weighted Raup-Crick dissimilarities 
regional.abunds <- t(as.matrix(colSums(OTUs)))
regional.relabunds <- decostand(regional.abunds, method = "total")
occupancy.probs <- t(as.matrix(colSums(decostand(OTUs, method = "pa")) / nrow(OTUs)))
site.abunds <- rowSums(OTUs)
site.rich <- specnumber(OTUs)
a <- regional.relabunds * occupancy.probs

# Create a null community based on Stegen et al. 2015
r <- nrow(OTUs)
c <- ncol(OTUs)
spec.vec <- 1:ncol(OTUs)
RCbc.nulls <- array(NA, c(r, c, 999))

# stochastic community assembly nulls
# for(i in 1:999){
#   if(i == 1) pb <- progress_bar$new(total = 999, force = T)
#   pb$update(ratio = i/999)
#   
#   null.comm <- OTUs * 0
#   # for first simulation:
#   for(row.i in 1:nrow(null.comm)){
#     #print(paste("run :", i, " -> ", row.i, " : ", site.abunds[row.i], " inds"))
#     
#     while(rowSums(null.comm)[row.i] < site.abunds[row.i]){
#       
#       
#       # choose a species based on its occupancy
#       local.specs <- sample(x = spec.vec, size = site.rich[row.i],
#                             prob = as.vector(occupancy.probs), replace = FALSE)
#       
#       local.probs <- decostand(t(as.matrix(regional.abunds[,local.specs])), method = "total")
#       
#       local.inds <- sample(x = local.specs, size = site.abunds[row.i],
#                            prob = as.vector(local.probs), replace = TRUE)
#       
#       local.abunds <- rle(sort(local.inds))
#       
#       # add an individual to the local community
#       null.comm[row.i, local.abunds$values] <- local.abunds$lengths
#     }
#   }
#   null.bc <- as.matrix(vegdist(decostand(null.comm, method = "total"), method = "bray"))
#   RCbc.nulls[,,i] <- null.bc
# }
# saveRDS(RCbc.nulls, file = "data/null_models/RCbc.null.rda")
RCbc.nulls <- readRDS(file = "data/null_models/RCbc.null.rda")
obs.bc <- as.matrix(vegdist(OTUsREL, method = "bray"))
site.compares <- expand.grid(site1 = 1:r, site2 = 1:r)
RC.bray <- matrix(NA, nrow = r, ncol = r)

for(row.i in 1:nrow(site.compares)){
  site1 <- site.compares[row.i,1]
  site2 <- site.compares[row.i,2]
  pairwise.null <- RCbc.nulls[site1,site2,]
  pairwise.bray <- obs.bc[site1,site2]
  num.greater <- sum(pairwise.null > pairwise.bray)
  num.ties <- sum(pairwise.null == pairwise.bray)
  val <- (((1 * num.greater) + (0.5 * num.ties))/1000 - 0.5) * 2
  RC.bray[site1, site2] <- val
}
rownames(RC.bray) <- rownames(design)
colnames(RC.bray) <- rownames(design)
RC.bray.dist <- as.dist(RC.bray)


# write.csv(RC.bray, "data/RCbray.csv")
# saveRDS(RC.bray.dist, "data/RCbraydist.rda")
rc.water <- as.dist(RC.bray[which(design$habitat == "water"), which(design$habitat == "water")])
rc.sed <- as.dist(RC.bray[which(design$habitat == "sediment"), which(design$habitat == "sediment")])
hist(rc.water, breaks = 30)
hist(rc.sed, breaks = 30)

rc.dat <- rbind.data.frame(
  cbind(RC.bray = liste(rc.water)[3], Habitat = rep("Planktonic")),
  cbind(RC.bray = liste(rc.sed)[3], Habitat = rep("Sediment")))
names(rc.dat)[1] <- "rc.bray"


# bNTI analysis
matched.phylo <- readRDS("data/matched.phylo.rda")
hja.comm <- matched.phylo$comm
hja.phy <- matched.phylo$phy

# mntd.hja <- comdistnt(hja.comm, cophenetic(hja.phy), abundance.weighted = T)
# saveRDS(mntd.hja, file = "data/mntds.rda")

# # Create null comms
# mntd.null <- array(NA, c(50, 50, 999))
# for(i in 1:999){
#   if(i == 1) pb <- progress_bar$new(total = 999, force = T)
#   pb$update(ratio = i/999)
#   #print(paste("creating null community ", i, " of 999"))
#   temp.mntd <- comdistnt(hja.comm,
#                          cophenetic(tipShuffle(hja.phy)),
#                          abundance.weighted=T)
#   mntd.null[,,i] <- as.matrix(temp.mntd)
#   if(i %% 50 == 0 | i == 999) saveRDS(mntd.null, file = "data/mntds-null-dist.rda")
# }

# read null dists 
mntd.hja <- readRDS(file = "data/mntds.rda")
mntds.null <- readRDS(file = "data/mntds-null-dist.rda")

obs.mntds <- as.matrix(mntd.hja)
site.compares <- expand.grid(site1 = 1:ncol(obs.mntds), site2 = 1:ncol(obs.mntds))
bNTI <- matrix(NA, nrow = nrow(obs.mntds), ncol = ncol(obs.mntds))
for(row.i in 1:nrow(site.compares)){
  site1 <- site.compares[row.i,1]
  site2 <- site.compares[row.i,2]
  pairwise.null <- mntds.null[site1,site2,]
  pairwise.mntd <- obs.mntds[site1,site2]
  null.mean <- mean(pairwise.null, na.rm = TRUE)
  null.sd <- sd(pairwise.null, na.rm = TRUE)
  val <- (pairwise.mntd - null.mean) / null.sd
  bNTI[site1, site2] <- val
}
colnames(bNTI) <- rownames(hja.comm)
rownames(bNTI) <- rownames(hja.comm)

bNTI.dist <- as.dist(bNTI)
hist(bNTI.dist, breaks = 30)
sum(bNTI.dist < 2 & bNTI.dist > -2) / length(bNTI.dist) # undom
sum(bNTI.dist > 2) / length(bNTI.dist) # variable selection
sum(bNTI.dist < -2) / length(bNTI.dist) # homogeneous selection


hja.bnti.dist.ls <- liste(bNTI.dist, entry = "bNTI")
hja.rcbray.dist.ls <- liste(RC.bray.dist, entry = "RC.bray")
hja.assembly <- as.data.frame(cbind(hja.bnti.dist.ls, hja.rcbray.dist.ls))
hja.assembly <- full_join(hja.bnti.dist.ls, hja.rcbray.dist.ls)

hja.mechanism <- vector(length = nrow(hja.assembly))
for(row.i in 1:nrow(hja.assembly)){
  if(hja.assembly[row.i,"bNTI"] < -2){
    hja.mechanism[row.i] <- str_wrap("Selection (Convergent)", width = 12)
  }
  if(hja.assembly[row.i,"bNTI"] > 2){
    hja.mechanism[row.i] <- str_wrap("Selection (Divergent)", width = 12)
  }
  if(hja.assembly[row.i,"bNTI"] <= 2 & hja.assembly[row.i,"bNTI"] >= -2){
    
    if (hja.assembly[row.i,"RC.bray"] <= 0.95 & hja.assembly[row.i,"RC.bray"] >= -0.95){
    hja.mechanism[row.i] <- str_wrap("Undominated", width = 12)
    }
    if(hja.assembly[row.i,"RC.bray"] < -.95){
    hja.mechanism[row.i] <- str_wrap("Mass Effects", width = 12)
    }
    if(hja.assembly[row.i,"RC.bray"] > .95){
    hja.mechanism[row.i] <- str_wrap("Dispersal Limitation", width = 12)
    }
  }
}
hja.assembly$Mechanism <- hja.mechanism



# identify comparisons

hja.assembly$comparison.habitat <- NA
hja.assembly[str_detect(hja.assembly$NBX, "_W") & str_detect(hja.assembly$NBY, "_W"),]$comparison.habitat <- "water_water" 
hja.assembly[str_detect(hja.assembly$NBX, "_W") & str_detect(hja.assembly$NBY, "_S"),]$comparison.habitat <- "water_sed" 
hja.assembly[str_detect(hja.assembly$NBX, "_S") & str_detect(hja.assembly$NBY, "_W"),]$comparison.habitat <- "water_sed" 
hja.assembly[str_detect(hja.assembly$NBX, "_S") & str_detect(hja.assembly$NBY, "_S"),]$comparison.habitat <- "sed_sed" 

hja.assembly$comparison.location <- NA
hja.assembly[design[hja.assembly$NBX,]$order < 2 & design[hja.assembly$NBY,]$order < 2,]$comparison.location <- "headwater_headwater"
hja.assembly[design[hja.assembly$NBX,]$order >= 2 & design[hja.assembly$NBY,]$order >= 2,]$comparison.location <- "downstream_downstream"
hja.assembly[design[hja.assembly$NBX,]$order < 2 & design[hja.assembly$NBY,]$order >= 2,]$comparison.location <- "headwater_downstream"
hja.assembly[design[hja.assembly$NBX,]$order >= 2 & design[hja.assembly$NBY,]$order < 2,]$comparison.location <- "headwater_downstream"

head(hja.assembly)

community.assembly <- hja.assembly
ggplot(data = community.assembly, aes(x = bNTI, y = RC.bray, col = Mechanism)) +
  geom_point(show.legend = T) + 
  geom_hline(yintercept = c(0.95, -0.95), lty = "dashed", col = "lightgrey")+
  geom_vline(xintercept = c(-2, +2), lty = "dashed", col = "lightgrey")+
  theme_cowplot()+
  labs(x = expression(paste(beta,"NTI")), y = expression(paste("Raup-Crick"[BC])))+
  labs(title = "Community Assembly Mechanisms")
community.assembly.plot
ggsave("figures/comm_assembly.pdf", width = 8, height = 6, units = "in")

assembly.table <- community.assembly %>% group_by(comparison.habitat) %>% count(Mechanism)
ggplot(assembly.table, aes(x = factor(Mechanism), y = n, fill = comparison.habitat)) + 
  geom_bar(stat = "identity") +
  theme_cowplot() + 
  labs(y = "Count", x = "Community Assembly Mechanism") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))
assembly.counts

community.assembly %>% group_by(comparison.habitat, comparison.location) %>% 
  count(Mechanism) %>% complete(comparison.habitat, comparison.location)
  
  mutate(tot = sum(n)) %>% mutate(proportion = n/tot) %>%
  ggplot(aes(x = factor(Mechanism), y = proportion, fill = comparison.habitat, alpha = comparison.location)) + 
  geom_bar(stat = "identity", position = 'dodge') +
  theme_classic() + 
  labs(y = "Count", x = "Community Assembly Mechanism") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))


assembly.table <- community.assembly %>% group_by(comparison.habitat, comparison.location) %>% count(Mechanism)
ggplot(assembly.table, aes(x = factor(Mechanism), y = n, fill = comparison.habitat, alpha = comparison.location)) + 
  geom_bar(stat = "identity") +
  theme_cowplot() + 
  labs(y = "Count", x = "Community Assembly Mechanism") +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14))

