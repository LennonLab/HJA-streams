# Functions used in hja-streams.Rmd

se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


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
  spa.detrended.comm.dist <- residuals(lm(comm.dist ~ spa.dist))
  spa.detrended.env.dist <- residuals(lm(env.dist ~ spa.dist))
  
  # now DDR for env.
  ddr.env <- lm(spa.detrended.comm.dist ~ spa.detrended.env.dist)
  
  # env detrend
  env.detrended.comm.dist <- residuals(lm(comm.dist ~ env.dist))
  env.detrended.spa.dist <- residuals(lm(spa.dist ~ env.dist))
  
  # now DDR for space
  ddr.space <- lm(env.detrended.comm.dist ~ env.detrended.spa.dist)
  
  out.models$spatial <- ddr.space
  out.models$env <- ddr.env
  
  plot(spa.detrended.comm.dist ~ spa.detrended.env.dist, ylab = comm, xlab = env)
  abline(ddr.env, lwd = 2)
  lines(smooth.spline(spa.detrended.env.dist, spa.detrended.comm.dist), col = "red", lwd = 2)
  print(summary(ddr.env))
  plot(env.detrended.comm.dist ~ env.detrended.spa.dist, ylab = comm, xlab = space)
  abline(ddr.space, lwd = 2)
  lines(smooth.spline(env.detrended.spa.dist, env.detrended.comm.dist), col = "red", lwd = 2)
  print(summary(ddr.space))
  
  detrended.dists <- list(env.detrended.comm.dist = env.detrended.comm.dist, 
                          spa.detrended.comm.dist = spa.detrended.comm.dist,
                          env.detrended.spa.dist = env.detrended.spa.dist,
                          spa.detrended.env.dist = spa.detrended.env.dist)
  
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
       plot(x = detrended.dists$spa.detrended.env.dist, 
            y = detrended.dists$spa.detrended.comm.dist, 
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
       mtext("Community Similarity", side = 2, line = 3.5, cex = 1.5)
       mtext("Environmental Distance", side = 1, line = 3, cex = 1.5)
       
       if(my.p != "N.S.") abline(models$env, lwd = 2)
       
       # plot spatial distance decay
       plot(x = detrended.dists$env.detrended.spa.dist, 
            y = detrended.dists$env.detrended.comm.dist, 
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
       mtext("Community Similarity", side = 2, line = 3.5, cex = 1.5)
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


normalize.matrix <- function(m){
  return((m - min(m)) / (max(m) - min(m)))
}

SimpE <- function(x = ""){
  x <- as.data.frame(x)
  D <- vegan::diversity(x, "inv")
  S <- sum((x > 0) * 1) 
  E <- (D)/S 
  return(E)
}


# Create dendritic distance matrix
make.dendritic.dists <- function(infile = "") {
  dend.dist.mat <-
    read.delim(
      file = infile,
      sep = ',',
      header = F,
      skip = 1,
      row.names = 1
    )
  colnames(dend.dist.mat) <-
    as.vector(lapply(strsplit(rownames(dend.dist.mat), " "), function(x)
      x[1]))
  rownames(dend.dist.mat) <- colnames(dend.dist.mat)
  empty.den.dist.mat <-
    matrix(NA, nrow = nrow(OTUs), ncol = nrow(OTUs))
  rownames(empty.den.dist.mat) <- env$sample
  colnames(empty.den.dist.mat) <- env$sample
  
  
  dend.dist.mat <- as.matrix(as.dist(dend.dist.mat))
  dend.dist.list <- simba::liste(dend.dist.mat)
  dend.dists <- simba::liste(vegdist(OTUs)) # to know which samples are remaining in dataset
  dend.dists[, 3] <- NA
  
  i <- 1
  
  for (this.row in rownames(empty.den.dist.mat)) {
    for (this.col in colnames(empty.den.dist.mat)) {
      row.unpack <- unlist(strsplit(this.row, "_"))
      col.unpack <- unlist(strsplit(this.col, "_"))
      
      row.site <- paste(row.unpack[1], "_", row.unpack[2], sep = "")
      col.site <- paste(col.unpack[1], "_", col.unpack[2], sep = "")
      
      empty.den.dist.mat[this.row, this.col] <-
        dend.dist.mat[row.site, col.site]
    }
    
  }
  
  den.dists <- empty.den.dist.mat
  den.dists <-
    den.dists[which(rownames(den.dists) %in% rownames(OTUs)),
              which(rownames(den.dists) %in% rownames(OTUs))]
  
  return(as.dist(den.dists))
}

# Rescale variables
scale_vec <- function(my_var){
  (my_var - mean(my_var)) / sd(my_var)
}

# ellipse functions 
# https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

calc.ellipse <- function(ord, ellipse){
  df_ell <- data.frame()
  for(g in levels(ord$group)){
    df_ell <- rbind(df_ell, 
                    cbind(as.data.frame(
                      with(ord[ord$group==g,], 
                           veganCovEllipse(ellipse[[g]]$cov, ellipse[[g]]$center, ellipse[[g]]$scale))), group = g))
  }
  return(df_ell)
  
}
