# source("./analysis/InitialSetup.R")
# 
# require("vegan")
# require("vegetarian")

#### Diversity Analysis

### Richness

# Observed Richness
S.obs <- rowSums((OTUs > 0) * 1)

# Simpson's Evenness
SimpE <- function(x = ""){
  x <- as.data.frame(x)
  D <- vegan::diversity(x, "inv")
  S <- sum((x > 0) * 1) 
  E <- (D)/S 
  return(E)
}
simpsE <- round(apply(OTUsREL, 1, SimpE), 3)
shan <- vegan::diversity(OTUsREL, index = "shannon")
N1 <- exp(shan)
simpsD <- vegan::diversity(OTUsREL, index = "invsimpson")
E1 <- N1/S.obs

# Rarefaction
min.N <- min(rowSums(OTUs))
S.rarefy <- rarefy(x = OTUs, sample = (min.N), se = TRUE)
#rarecurve(x = OTUs, step = 20, col = "blue", cex = 0.6, las = 1)
rare.d <- t(S.rarefy)
colnames(rare.d) <- c("S.rare", "se.rare")
alpha.div <- cbind(design, S.obs, simpsE, shan, N1, simpsD, rare.d)


num.coms <- 10
rarecoms <- vector("list", num.coms)
for(each in 1:num.coms){
  rarecoms[[each]] <- rrarefy(OTUs, sample = min.N)
  #rarecoms[[each]] <- OTUs[1:10,1:10]
}

dists <- lapply(X=rarecoms, FUN = vegdist, method = "bray")


# Seperate data based on water and sediment samples
alpha.water <- alpha.div[alpha.div$habitat == "water",]
alpha.sed <- alpha.div[alpha.div$habitat == "sediment", ]

# Differences in diversity (richness and evenness) between habitats
t.test(alpha.water$S.rare, alpha.sed$S.rare)
t.test(alpha.water$N1, alpha.sed$N1)
boxplot(alpha.water$N1, alpha.sed$N1)
capture.output(summary(lm(alpha.div$S.rare ~ design$habitat == "water")),
               file = "./tables/richness_compare_model.txt")

### Diversity Partitioning
s.size <- as.vector(c(rowSums(OTUsREL)))
hja.alpha <- d(OTUsREL, lev = "alpha", q = 1, boot = T, 
               boot.arg = list(s.sizes = s.size))
hja.beta <- d(OTUsREL, lev = "beta", q = 1, boot = T, 
              boot.arg = list(s.sizes = s.size))
hja.gamma <- d(OTUsREL, lev = "gamma", q = 1, boot = T, 
               boot.arg = list(s.sizes = s.size))

water.alpha <- d(OTUsREL[which(design$habitat == "water"),], 
                 lev = "alpha", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "water")]))
sediment.alpha <- d(OTUsREL[which(design$habitat == "sediment"),], 
                    lev = "alpha", q = 1, boot = T, 
                    boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment")]))

water.beta <- d(OTUsREL[which(design$habitat == "water"),], 
                lev = "beta", q = 1, boot = T, 
                boot.arg = list(s.sizes = s.size[which(design$habitat == "water")]))
sediment.beta <- d(OTUsREL[which(design$habitat == "sediment"),], 
                   lev = "beta", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment")]))

water.alpha
sediment.alpha
water.alpha$D.Value / sediment.alpha$D.Value
water.beta$D.Value / sediment.beta$D.Value

# Calculate Water Alpha Diversity
water.alpha.1 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 1),], 
                   lev = "alpha", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 1)]))
water.alpha.2 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 2),], 
                   lev = "alpha", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 2)]))
water.alpha.3 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 3),], 
                   lev = "alpha", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 3)]))
water.alpha.4 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 4),], 
                   lev = "alpha", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 4)]))
water.alpha.5 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 5),], 
                   lev = "alpha", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 5)]))
# Calculate Water Beta Diversity
water.beta.1 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 1),], 
                  lev = "beta", q = 1, boot = T, 
                  boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                           design$order == 1)]))
water.beta.2 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 2),], 
                  lev = "beta", q = 1, boot = T, 
                  boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                           design$order == 2)]))
water.beta.3 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 3),], 
                  lev = "beta", q = 1, boot = T, 
                  boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                           design$order == 3)]))
water.beta.4 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 4),], 
                  lev = "beta", q = 1, boot = T, 
                  boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                           design$order == 4)]))
water.beta.5 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 5),], 
                  lev = "beta", q = 1, boot = T, 
                  boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                           design$order == 5)]))

# Calculate Water Gamma Diversity
water.gamma.1 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 1),], 
                   lev = "gamma", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 1)]))
water.gamma.2 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 2),], 
                   lev = "gamma", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 2)]))
water.gamma.3 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 3),], 
                   lev = "gamma", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 3)]))
water.gamma.4 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 4),], 
                   lev = "gamma", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 4)]))
water.gamma.5 <- d(OTUsREL.log[which(design$habitat == "water" & design$order == 5),], 
                   lev = "gamma", q = 1, boot = T, 
                   boot.arg = list(s.sizes = s.size[which(design$habitat == "water" & 
                                                            design$order == 5)]))

# Calculate sediment alpha diversity
sed.alpha.1 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 1),], 
                 lev = "alpha", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 1)]))
sed.alpha.2 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 2),], 
                 lev = "alpha", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 2)]))
sed.alpha.3 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 3),], 
                 lev = "alpha", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 3)]))
sed.alpha.4 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 4),], 
                 lev = "alpha", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 4)]))
sed.alpha.5 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 5),], 
                 lev = "alpha", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 5)]))

# Calculate sediment beta diversity
sed.beta.1 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 1),], 
                lev = "beta", q = 1, boot = T, 
                boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                         design$order >= 1)]))
sed.beta.2 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 2),], 
                lev = "beta", q = 1, boot = T, 
                boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                         design$order >= 2)]))
sed.beta.3 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 3),], 
                lev = "beta", q = 1, boot = T, 
                boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                         design$order >= 3)]))
sed.beta.4 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 4),], 
                lev = "beta", q = 1, boot = T, 
                boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                         design$order >= 4)]))
sed.beta.5 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 5),], 
                lev = "beta", q = 1, boot = T, 
                boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                         design$order >= 5)]))

# Calculate Sediment Gamma Diversity
sed.gamma.1 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 1),], 
                 lev = "gamma", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 1)]))
sed.gamma.2 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 2),], 
                 lev = "gamma", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 2)]))
sed.gamma.3 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 3),], 
                 lev = "gamma", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 3)]))
sed.gamma.4 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 4),], 
                 lev = "gamma", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 4)]))
sed.gamma.5 <- d(OTUsREL.log[which(design$habitat == "sediment" & design$order == 5),], 
                 lev = "gamma", q = 1, boot = T, 
                 boot.arg = list(s.sizes = s.size[which(design$habitat == "sediment" & 
                                                          design$order >= 5)]))

# Create Vectors of Diversity Measures
water.a <- c(water.alpha.1$D.Value, water.alpha.2$D.Value, water.alpha.3$D.Value, 
             water.alpha.4$D.Value, water.alpha.5$D.Value)
water.a.se <- c(water.alpha.1$StdErr, water.alpha.2$StdErr, water.alpha.3$StdErr, 
                water.alpha.4$StdErr, water.alpha.5$StdErr)
water.b <- c(water.beta.1$D.Value, water.beta.2$D.Value, water.beta.3$D.Value, 
             water.beta.4$D.Value, water.beta.5$D.Value)
water.b.se <- c(water.beta.1$StdErr, water.beta.2$StdErr, water.beta.3$StdErr, 
                water.beta.4$StdErr, water.beta.5$StdErr)
water.g <- c(water.gamma.1$D.Value, water.gamma.2$D.Value, water.gamma.3$D.Value, 
             water.gamma.4$D.Value, water.gamma.5$D.Value)
water.g.se <- c(water.gamma.1$StdErr, water.gamma.2$StdErr, water.gamma.3$StdErr, 
                water.gamma.4$StdErr, water.gamma.5$StdErr)

sed.a <- c(sed.alpha.1$D.Value, sed.alpha.2$D.Value, sed.alpha.3$D.Value, 
           sed.alpha.4$D.Value, sed.alpha.5$D.Value)
sed.a.se <- c(sed.alpha.1$StdErr, sed.alpha.2$StdErr, sed.alpha.3$StdErr, 
              sed.alpha.4$StdErr, sed.alpha.5$StdErr)
sed.b <- c(sed.beta.1$D.Value, sed.beta.2$D.Value, sed.beta.3$D.Value, 
           sed.beta.4$D.Value, sed.beta.5$D.Value)
sed.b.se <- c(sed.beta.1$StdErr, sed.beta.2$StdErr, sed.beta.3$StdErr, 
              sed.beta.4$StdErr, sed.beta.5$StdErr)
sed.g <- c(sed.gamma.1$D.Value, sed.gamma.2$D.Value, sed.gamma.3$D.Value, 
           sed.gamma.4$D.Value, sed.gamma.5$D.Value)
sed.g.se <- c(sed.gamma.1$StdErr, sed.gamma.2$StdErr, sed.gamma.3$StdErr, 
              sed.gamma.4$StdErr, sed.gamma.5$StdErr)
orders <- c(1,2,3,4,5)
water.hill <- as.data.frame(cbind(orders, water.a, water.a.se, water.b, water.b.se, water.g, water.g.se))
sed.hill <- as.data.frame(cbind(orders, sed.a, sed.a.se, sed.b, sed.b.se, sed.g, sed.g.se))

# Plot by stream order to visualize
png(filename = "./figures/DiversityPartitioning_Order.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(water.hill$orders, water.hill$water.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "black",
     xaxt = "n", yaxt = "n", ylim = c(0, 12000),
     xlab = "", ylab = "")
points(water.hill$orders, water.hill$water.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "black")
points(sed.hill$orders, sed.hill$sed.a, type = "o", lty = 2, lwd = 2, col = "grey70", pch = 22)
points(sed.hill$orders, sed.hill$sed.g, type = "o", lty = 3, lwd = 2, col = "grey70", pch = 23)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
par(new = TRUE)
plot(water.hill$orders, water.hill$water.b, type = "o", lwd = 2, pch = 21, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(1,2.5))
points(sed.hill$orders, sed.hill$sed.b, type = "o", lwd = 2, col = "grey", pch = 21)
axis(side = 4, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(alpha," and ",gamma,"-diversity")), side = 2, line = 3.5, cex = 1.5)
mtext(expression(paste(beta,"-diversity")), side = 4, line = 3.5, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-alpha.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(water.hill$orders, water.hill$water.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "black",
     xaxt = "n", yaxt = "n", ylim = c(0, 12000),
     xlab = "", ylab = "")
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(alpha,"-diversity")), side = 2, line = 3.5, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-beta.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(water.hill$orders, water.hill$water.b, type = "o", lwd = 2, pch = 21, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(1,3))
error.bar(water.hill$orders, water.hill$water.b, upper = water.hill$water.b.se)
axis(side = 4, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(beta,"-diversity")), side = 4, line = 3.5, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-gamma.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(water.hill$orders, water.hill$water.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 12000))
error.bar(water.hill$orders, water.hill$water.g, upper = water.hill$water.g.se)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(gamma,"-diversity")), side = 2, line = 3.5, cex = 1.5)
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-water.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(water.hill$orders, water.hill$water.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "black",
     xaxt = "n", yaxt = "n", ylim = c(0, 12000),
     xlab = "", ylab = "")
points(water.hill$orders, water.hill$water.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "black")
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
par(new = TRUE)
plot(water.hill$orders, water.hill$water.b, type = "o", lwd = 2, pch = 21, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(1,3))
error.bar(water.hill$orders, water.hill$water.b, upper = water.hill$water.b.se)
axis(side = 4, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(alpha," and ",gamma,"-diversity")), side = 2, line = 3.5, cex = 1.5)
mtext(expression(paste(beta,"-diversity")), side = 4, line = 3.5, cex = 1.5)
legend(x = 3.7, y = 3, pch = c(22, 21, 23), lty = c(2, 2, 3), lwd = c(2,2,2),
       pt.bg = c("black", "black", "black"),
       legend = c(expression(paste(alpha,"-diversity")), 
                  expression(paste(beta,"-diversity")),
                  expression(paste(gamma,"-diversity"))), bty = "n")
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-sed.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(sed.hill$orders, sed.hill$sed.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "black",
     xaxt = "n", yaxt = "n", ylim = c(0, 9000),
     xlab = "", ylab = "")
points(sed.hill$orders, sed.hill$sed.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "black")
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
par(new = TRUE)
plot(sed.hill$orders, sed.hill$sed.b, type = "o", lwd = 2, pch = 21, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(1,2.5))
error.bar(sed.hill$orders, sed.hill$sed.b, upper = sed.hill$sed.b.se)
axis(side = 4, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(alpha," and ",gamma,"-diversity")), side = 2, line = 3.5, cex = 1.5)
mtext(expression(paste(beta,"-diversity")), side = 4, line = 3.5, cex = 1.5)
legend(x = 3.7, y = 2.5, pch = c(22, 21, 23), lty = c(2, 2, 3), lwd = c(2,2,2),
       pt.bg = c("black", "black", "black"),
       legend = c(expression(paste(alpha,"-diversity")), 
                  expression(paste(beta,"-diversity")),
                  expression(paste(gamma,"-diversity"))), bty = "n")
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-total.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(water.hill$orders, water.hill$water.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "black",
     xaxt = "n", yaxt = "n", ylim = c(0, 12000),
     xlab = "", ylab = "")
points(water.hill$orders, water.hill$water.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "black")
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
points(sed.hill$orders, sed.hill$sed.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "grey50", col = "grey50")
points(sed.hill$orders, sed.hill$sed.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "grey50", col = "grey50")

par(new = TRUE)
plot(water.hill$orders, water.hill$water.b, type = "o", lwd = 2, pch = 21, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(1,3))
error.bar(water.hill$orders, water.hill$water.b, upper = water.hill$water.b.se)
points(sed.hill$orders, sed.hill$sed.b, type = "o", lwd = 2, pch = 21, bg = "grey50", col = "grey50")
error.bar(sed.hill$orders, sed.hill$sed.b, upper = sed.hill$sed.b.se, col = "grey50")

axis(side = 4, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(alpha," and ",gamma,"-diversity")), side = 2, line = 3.5, cex = 1.5)
mtext(expression(paste(beta,"-diversity")), side = 4, line = 3.5, cex = 1.5)
legend(x = 3.7, y = 3, pch = c(22, 21, 23), lty = c(2, 2, 3), lwd = c(2,2,2),
       pt.bg = c("black", "black", "black"),
       legend = c(expression(paste(alpha,"-diversity")), 
                  expression(paste(beta,"-diversity")),
                  expression(paste(gamma,"-diversity"))), bty = "n")
dev.off()
graphics.off()

png(filename = "./figures/DiversityPartitioning_Order-sed.png",
    height = 1200, width = 1200, res = 96*2)
par(mar = c(5,5,5,5) + 0.5)
plot(sed.hill$orders, sed.hill$sed.a, type = "o", lty = 2, lwd = 2, pch = 22, bg = "black",
     xaxt = "n", yaxt = "n", ylim = c(0, 1000),
     xlab = "", ylab = "")
points(sed.hill$orders, sed.hill$sed.g, type = "o", lty = 3, lwd = 2, pch = 23, bg = "black")
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
par(new = TRUE)
plot(sed.hill$orders, sed.hill$sed.b, type = "o", lwd = 2, pch = 21, bg = "black",
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(1,2.5))
error.bar(sed.hill$orders, sed.hill$sed.b, upper = sed.hill$sed.b.se)
axis(side = 4, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
mtext("Stream Order", side = 1, line = 3, cex = 1.5)
mtext(expression(paste(alpha," and ",gamma,"-diversity")), side = 2, line = 3.5, cex = 1.5)
mtext(expression(paste(beta,"-diversity")), side = 4, line = 3.5, cex = 1.5)
legend(x = 3.7, y = 2.5, pch = c(22, 21, 23), lty = c(2, 2, 3), lwd = c(2,2,2),
       pt.bg = c("black", "black", "black"),
       legend = c(expression(paste(alpha,"-diversity")), 
                  expression(paste(beta,"-diversity")),
                  expression(paste(gamma,"-diversity"))), bty = "n")
dev.off()
graphics.off()

model.water.b <- lm(water.hill$water.b ~ water.hill$orders)
summary(model.water.b)
conf95.water.b <- predict.lm(object = model.water.b, interval="confidence")

model.sed.b <- lm(sed.hill$sed.b ~ sed.hill$orders)
summary(model.sed.b)
conf95.sed.b <- predict.lm(object = model.sed.b, interval="confidence")




### Calculate Bray-Curtis across stream orders
order.1.bray.w <- vegdist(OTUsREL.log[which(design$order==1 & design$habitat == "water"),], method = "bray")
order.2.bray.w <- vegdist(OTUsREL.log[which(design$order>=2 & design$habitat == "water"),], method = "bray")

order.1.bray.s <- vegdist(OTUsREL.log[which(design$order==1 & design$habitat == "sediment"),], method = "bray")
order.2.bray.s <- vegdist(OTUsREL.log[which(design$order>=2 & design$habitat == "sediment"),], method = "bray")

order.1.bray.w.list <- liste(order.1.bray.w)[,3]
order.2.bray.w.list <- liste(order.2.bray.w)[,3]

order.1.bray.s.list <- liste(order.1.bray.s)[,3]
order.2.bray.s.list <- liste(order.2.bray.s)[,3]

orders.bray.w <- as.data.frame(cbind(c(
  rep(1, length(order.1.bray.w.list)),
  rep(2, length(order.2.bray.w.list))),
  c(order.1.bray.w.list, order.2.bray.w.list)))
colnames(orders.bray.w) <- c("order", "dist")
model.orders.w <- lm(orders.bray.w$dist ~ orders.bray.w$order)
(modsum.w <- summary(model.orders.w))

orders.bray.s <- as.data.frame(cbind(c(
  rep(1, length(order.1.bray.s.list)),
  rep(2, length(order.2.bray.s.list))),
  c(order.1.bray.s.list, order.2.bray.s.list)))
colnames(orders.bray.s) <- c("order", "dist")
model.orders.s <- lm(orders.bray.s$dist ~ orders.bray.s$order)
(modsum.s <- summary(model.orders.s))
# model.orders.w <- aov(orders.bray.w$dist ~ as.factor(orders.bray.w$order))
# (modsum.w <- summary(model.orders.w))
# TukeyHSD(model.orders.w)
# 
# model.orders.s <- aov(orders.bray.s$dist ~ as.factor(orders.bray.s$order))
# (modsum.s <- summary(model.orders.s))
# TukeyHSD(model.orders.s)

#### Figures

### Figure 1: Rarefied diversity in sediment and surface-water communities
png(filename = "./figures/AlphaDiv.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 3, 2) + 0.3)
boxplot(S.rare ~ habitat, data = alpha.div, at = c(3,1),
        names = c("Sediment", "Water"), col = c("grey", "white"),
        xlim = c(0,4),
        yaxt = "n", ylab = "", xlab = "", cex.axis = 1.2)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 1, labels = F, at = c(1,3), lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd = 2)
text(3, 1700, "*", cex = 2)
mtext("Habitat Type", side = 1, line = 3, cex = 1.5)
mtext("Species Richness", side = 2, line = 4, cex = 1.5)
dev.off()
graphics.off()
#img <- readPNG("./figures/AlphaDiv.png")
#grid.raster(img)

# Figure 3: Beta-diversity is higher among headwaters than among higher order streams.
png(filename = "./figures/BetaDiv_BrayCurtis.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5, 5, 2, 2) + 0.1)
boxplot(dist ~ order, data = orders.bray.w, at = c(1,4), xlim = c(0,6), ylim = c(0.4, 0.9), 
        col = "white", xlab = "", xaxt = "n", yaxt = "n")
boxplot(dist ~ order, data = orders.bray.s, at = c(2,5), add = T, col = "grey",
        xlab = "", xaxt = "n", yaxt = "n")
axis(side = 1, labels = F, at = c(1.5, 4.5), lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
box(lwd=2)
mtext(side = 2, cex = 1.5, line = 3.5, expression(paste(beta,"-diversity")))
mtext(side = 1, cex = 1.5, line = 3, "Stream Order")
mtext(side = 1, at = c(1.5,4.5), line = 1, c("Headwaters","Higher-order"), cex = 1.2)

legend("bottomleft", c("Water", "Sediment"),
       pt.bg = c("white", "grey"), pch = c(22,22), cex = 1.5, bty = "n")
dev.off()
graphics.off()

img <- readPNG("./figures/BetaDiv_BrayCurtis.png")
grid.raster(img)



png(filename = "./figures/BetaDiv_NumEquiv.png",
    width = 1200, height = 1200, res = 96*2)
par(mar = c(5,5,2,1) + 0.1)

plot(water.hill$orders, water.hill$water.b, type = "l", lty = 1, lwd = 2,
     ylim = c(1,2.5), ylab = "", xlab = "", xaxt = "n", yaxt = "n")
#matlines(water.hill$orders, conf95.water.b, lty = c(1,0,0),
#         col = c("black", "gray50", "gray50"), lwd = c(2,1,1))
points(sed.hill$orders, sed.hill$sed.b, type = "l", lty = 2, lwd = 2)
#matlines(sed.hill$orders, conf95.sed.b, lty = c(1,0,0),
#         col = c("black", "gray50", "gray50"), lwd = c(2,1,1))
#points(water.hill$orders, (water.hill$water.b + water.hill$water.b.se), type = "l", lty = 3)
#points(water.hill$orders, (water.hill$water.b - water.hill$water.b.se), type = "l", lty = 3)
#points(sed.hill$orders, (sed.hill$sed.b + sed.hill$sed.b.se), type = "l", lty = 3)
#points(sed.hill$orders, (sed.hill$sed.b - sed.hill$sed.b.se), type = "l", lty = 3)
axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
box(lwd=2)
mtext(side = 2, cex = 1.5, line = 3, expression(paste(beta,"-diversity")))
mtext(side = 1, cex = 1.5, line = 3, "Stream Order")

legend(1, 1.35, c("Bacterioplankton", "Sediment-associated"), 
       lty = c(1,2), lwd = c(2,2), bty = "n", cex = 1.5)

dev.off()
graphics.off()

img <- readPNG("./figures/BetaDiv_NumEquiv.png")
grid.raster(img)