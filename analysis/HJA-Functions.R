# Functions used in hja-streams.Rmd

se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


# resample.comms <- function(comm = "", n = 100, ...){
#   # commsize <- dim(comm)
#   # comms <- array(data = NA, dim = c(commsize[1:2], n))
#   # for (i in 1:n) {
#   #   comms[,,i] <- decostand(
#   #     rrarefy(OTUs, sample = min(rowSums(OTUs))),
#   #     method = transformation)
#   #   
#   # }
#   coefs <- matrix(NA, nrow = n, ncol = 2)
#   
#   for(i in 1:n){
#     comm.i <- decostand(
#       rrarefy(OTUs, sample = min(rowSums(OTUs))),
#       method = "hellinger")
#     comm.dist <- vegdist(comm.i, method = "euclidean")
#     ddr <- lm(comm.dist ~ den.dists)
#     coefs[i, ] <- coefficients(ddr)
#   }
#   
#   
#   
#   return(coefs)
# }

####

run.pcoa <- function(comm = NULL, dist.metric = "bray", plot = T, ...){
  dist.matrix <- vegdist(comm, method = dist.metric)
  pcoa <- cmdscale(dist.matrix, eig = TRUE)
  var1 <- round(pcoa$eig[1] / sum(pcoa$eig), 3) * 100
  var2 <- round(pcoa$eig[2] / sum(pcoa$eig), 3) * 100
  cat("PCoA Axis 1 explains", var1, "percent of total variation.\n", sep = " ")
  cat("PCoA Axis 2 explains", var2, "percent of total variation.\n", sep = " ")
  if(plot) ordiplot(pcoa)
  return(list(pcoa = pcoa, var1 = var1, var2 = var2, dist.matrix = dist.matrix))
}

# Function to extract explained variation from an ordination axis
explain.var <- function(ord = "", axis = 1){
  return(round((eigenvals(ord)[axis]/sum(eigenvals(ord))*100), 3))
}

add.axes <- function(main = NULL, xlab = NULL, ylab = NULL, s1 = T, s2 = T, s3 = T, s4 = T, ...){
  if(s1) axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  if(s2) axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  if(s3) axis(side = 3, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  if(s4) axis(side = 4, labels = F, lwd.ticks = 2, cex.axis = 1.2, las = 1)
  box(lwd = 2)
  mtext(xlab, side = 1, line = 3, cex = 1.5)
  mtext(ylab, side = 2, line = 3, cex = 1.5)
  mtext(main, side = 3, line = 2, cex = 2)
}




DDR <- function(dists = NULL, comm = "otus", env = "env", space = "den"){
  out.models <- list()
  comm.dist <- dists[[comm]]
  env.dist <- dists[[env]]
  spa.dist <- dists[[space]]
  
  # spatial detrend
  spatial.lm <- lm(comm.dist ~ spa.dist)
  spa.detrended.comm.dist <- residuals(spatial.lm) + coef(spatial.lm)[1]
  
  # now DDR for env.
  ddr.env <- lm(spa.detrended.comm.dist ~ env.dist)
  
  # env detrend
  env.lm <- lm(comm.dist ~ env.dist)
  env.detrended.comm.dist <- residuals(env.lm) + coef(env.lm)[1]
  
  # now DDR for space
  ddr.space <- lm(env.detrended.comm.dist ~ spa.dist)
  
  out.models$spatial <- ddr.space
  out.models$env <- ddr.env
  
  plot(spa.detrended.comm.dist ~ env.dist, ylab = comm, xlab = env)
  abline(ddr.env)
  print(summary(ddr.env))
  plot(env.detrended.comm.dist ~ spa.dist, ylab = comm, xlab = space)
  abline(ddr.space)
  print(summary(ddr.space))
  
  detrended.dists <- list(env.detrend = env.detrended.comm.dist, 
                          spatial.detrend = spa.detrended.comm.dist,
                          env.dist = env.dist, spatial.dist = spa.dist)
  
  out.models$detrended.dists <- detrended.dists
  return(out.models)
}


nb.calc <- function(species = ""){
  Pij2 = species^2
  B = 1/sum(Pij2)
  return(B)
}


plot.DDRs <- function(models = NULL){
  with(models, {
    
       # plot env distance decay
       plot(x = detrended.dists$env.dist, y = detrended.dists$spatial.detrend, 
            xaxt = "n", yaxt = "n", xlab = "", ylab = "")
       
       # add r2 and p val
       r2 <- round(summary(models$env)$r.squared, 3)
       my.p <- round(summary(models$env)$coefficients[2,4],3)
       if(my.p < 0.001) my.p <- "< 0.001"
       if(my.p > 0.05) my.p <- "N.S."
       rp = vector('expression',2)
       rp[1] = substitute(expression(italic(r)^2 == MYVALUE), 
                          list(MYVALUE = format(r2,dig=3)))[2]
       rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                          list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
       legend('topright', legend = rp, bty = 'n')
       
       axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
       axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
       axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
       axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
       box(lwd = 2)
       mtext("Community Dissimilarity", side = 2, line = 3, cex = 1.5)
       mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
       
       if(my.p != "N.S.") abline(models$env, lwd = 2)
       
       # plot spatial distance decay
       plot(x = detrended.dists$spatial.dist, y = detrended.dists$env.detrend, 
            xaxt = "n", yaxt = "n", xlab = "", ylab = "")
       
       # add r2 and p val
       r2 <- round(summary(models$spatial)$r.squared, 3)
       my.p <- round(summary(models$spatial)$coefficients[2,4],3)
       if(my.p < 0.001) my.p <- "< 0.001"
       if(my.p > 0.05) my.p <- "N.S."
       rp = vector('expression',2)
       rp[1] = substitute(expression(italic(r)^2 == MYVALUE), 
                          list(MYVALUE = format(r2,dig=3)))[2]
       rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                          list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
       legend('topright', legend = rp, bty = 'n')
       
       axis(side=1, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
       axis(side=2, labels=T, lwd.ticks=2, cex.axis=1.2, las=1)
       axis(side=3, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
       axis(side=4, labels=F, lwd.ticks=2, cex.axis=1.2, las=1)
       box(lwd = 2)
       mtext("Community Dissimilarity", side = 2, line = 3, cex = 1.5)
       mtext("Dendritic Distance", side = 1, line = 3, cex = 1.5)
       
       if(my.p != "N.S.") abline(models$spatial, lwd = 2)
       
       })
}


fill.table <- function(lms.in = NULL, response.metric = response.matrix, ddr.summary = ddr.summary){
  
  spatial.df <- data.frame(
    metric = response.metric,
    dat = deparse(substitute(lms.in)),
    response = as.character(formula(lms.in$spatial)[2]),
    predictor = as.character(formula(lms.in$spatial)[3]),
    int = coef(lms.in$spatial)[1],
    slope = coef(lms.in$spatial)[2],
    r2 = summary(lms.in$spatial)$r.squared,
    p = summary(lms.in$spatial)$coefficients[2,4]
  )
  
  env.df <- data.frame(
    metric = response.metric,
    dat = deparse(substitute(lms.in)),
    response = as.character(formula(lms.in$env)[2]),
    predictor = as.character(formula(lms.in$env)[3]),
    int = coef(lms.in$env)[1],
    slope = coef(lms.in$env)[2],
    r2 = summary(lms.in$env)$r.squared,
    p = summary(lms.in$env)$coefficients[2,4]
  )
  
  return(rbind.data.frame(ddr.summary, rbind.data.frame(spatial.df, env.df)))
}
