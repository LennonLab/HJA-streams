rm(list=ls())
setwd("~/GitHub/HJA-Streams/")

# This is a script to read in raw nutrient data from the HJ Andrews
# 2015 trip, calculate a standard curve and TN/TP from env samples.

# Import raw data
tn.samples <- read.csv("./data/water_chemistry/2015-07-17_HJA-TN.csv")
tn.stds <- read.csv("./data/water_chemistry/2015-07-17_HJA-TN_standards.csv")
tp.samples <- read.csv("./data/water_chemistry/2015-07-09_HJA-TP.csv")
tp.stds <- read.csv("./data/water_chemistry/2015-07-09_HJA-TP_standards.csv")

############ Calculate Total Nitrogen (TN) ##############

# Calculate TN standard curve (low values)
tn.stds.conc.low <- tn.stds[1:4,1]
tn.stds.abs.low <- tn.stds[1:4,2]
tn.linear <- lm(tn.stds.conc.low ~ tn.stds.abs.low)
summary(tn.linear)
plot(tn.stds.abs.low, tn.stds.conc.low, pch=16, main='TN Standard Curve',
     ylab='Total Nitrogen (µg/L)', xlab='Absorbance')
abline(tn.linear)

# Calculate TN standard curve (high values)
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

# Calculate TN in environmental samples
TN <- tn.samples[,2] * tn.linear$coefficients[2] + tn.linear$coefficients[1]
TN.out <- cbind(tn.samples, TN)

############## Calculate Total Phosphorus (TP) #################

# Calculate TP standard curve
tp.stds.conc <- tp.stds[,1]
tp.stds.abs <- tp.stds[,2]
tp.lm <- lm(tp.stds.conc ~ tp.stds.abs)
summary(tp.lm)
plot(tp.stds.abs, tp.stds.conc, pch=16, main='TP Standard Curve',
     ylab='Total Phosphorus (µg/L)', xlab='Absorbance')
abline(tp.lm)
text(1.5, 100, bquote(y == .(tp.lm$coefficients[2])*x + (.(tp.lm$coefficients[1]))))
text(1.5, 50, bquote('R'^2 == .(summary(tp.lm)$r.squared)))

# Calculate TP in environmental samples
TP <- tp.samples[,2] * tp.lm$coefficients[2] + tp.lm$coefficients[1]
TP.out <- cbind(tp.samples, TP)

############## Write data to File ####################
file.out <- cbind(TN.out[,3], TP.out[,3])
colnames(file.out) <- c("TN", "TP")
rownames(file.out) <- TN.out[,1]

read.csv('./data/')

write.table(file.out, './data/2015-07-27_water-chem-table.txt', sep='\t')
