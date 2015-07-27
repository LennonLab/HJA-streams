rm(list=ls())
setwd("~/GitHub/HJA-Streams/")

# Import raw data
tn.samples <- read.csv("./data/water_chemistry/2015-07-17_HJA-TN.csv")
tn.stds <- read.csv("./data/water_chemistry/2015-07-17_HJA-TN_standards.csv")
tp.samples <- read.csv("./data/water_chemistry/2015-07-09_HJA-TP.csv")
tp.stds <- read.csv("./data/water_chemistry/2015-07-09_HJA-TP_standards.csv")

############ Calculate Total Nitrogen (TN) ##############

# Calculate TN standard curve
tn.stds.conc <- tn.stds[,1]
tn.stds.abs <- tn.stds[,2]
tn.stds.abs.2 <- tn.stds.abs^2
quad.model <- lm(tn.stds.conc ~ tn.stds.abs + tn.stds.abs.2)
summary(quad.model)
plot(tn.stds.abs, tn.stds.conc, pch=16, ylab='Total Nitrogen (µg/L)', 
     xlab='Absorbance')

# Plot TN standard curve
abs.values <- seq(0, .03, 0.0001)
predicted.counts <- predict(quad.model,
                    list(tn.stds.abs=abs.values, tn.stds.abs.2=abs.values^2))
plot(tn.stds.abs, tn.stds.conc, main = "TN Standard Curve",
     xlab = "Absorbance", ylab = "Total Nitrogen (µg/L)", 
     cex.lab = 1.3, pch=16, col = "blue")
lines(abs.values, predicted.counts, col = "darkgreen", lwd = 3)


