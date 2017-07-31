# source("analysis/InitialSetup.R")
# source("analysis/DistanceCalcs.R")

# regional.abunds <- t(as.matrix(colSums(OTUs)))
# regional.relabunds <- decostand(regional.abunds, method = "total")
# occupancy.probs <- t(as.matrix(colSums(decostand(OTUs, method = "pa")) / nrow(OTUs)))
# site.abunds <- rowSums(OTUs)
# site.rich <- specnumber(OTUs)
# a <- regional.relabunds * occupancy.probs
# 
# # Create a null community based on Stegen et al. 2015
# 
# spec.vec <- 1:ncol(OTUs)
# RCbc.nulls <- array(NA, c(56, 56, 999))
# 
# # stochastic community assembly nulls
# for(i in 1:999){
#   
#   null.comm <- OTUs * 0
#   # for first simulation:
#   for(row.i in 1:nrow(null.comm)){
#     print(paste("run :", i, " -> ", row.i, " : ", site.abunds[row.i], " inds"))
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
site.compares <- expand.grid(site1 = 1:56, site2 = 1:56)
RC.bray <- matrix(NA, nrow = 56, ncol = 56)

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
hist(rc.water)
hist(rc.sed)

rc.dat <- rbind.data.frame(
  cbind(RC.bray = liste(rc.water)[3], Habitat = rep("Planktonic")),
  cbind(RC.bray = liste(rc.sed)[3], Habitat = rep("Sediment")))
names(rc.dat)[1] <- "rc.bray"

require(cowplot)
# 
# rc.bp <- ggplot(data = rc.dat, aes(x = habitat, y = rc.bray)) + 
#   geom_boxplot() + 
#   geom_hline(aes(yintercept = 0.95), lty = "dashed") +
#   geom_hline(aes(yintercept = -0.95), lty = "dashed") +
#   geom_jitter(aes(colour = abs(rc.bray) > 0.95)) + 
#   scale_color_manual(values = c('grey', 'red'))
# ggsave(filename = "figures/RC_bray_boxplot.pdf", plot = rc.bp, width = 8, height = 8, units = "in")

# RC.brayplot <- ggplot(data = rc.dat) +
#   geom_density(aes(x = rc.bray, y = ..scaled.., fill = Habitat), alpha = .5) +
#   labs(x = expression(paste("RC"["bray"])), y = "Frequency") + 
#   scale_fill_manual(values = c("skyblue", "wheat"))+
#   theme_cowplot() + 
#   theme(axis.title = element_text(size = 16), axis.text = element_text(size = 12))
# RC.brayplot
