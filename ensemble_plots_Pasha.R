######################################
# Arfa Pasha
# Ensemble plots 
# BISC 481
#10.25.2016
######################################

# Initialization
library(DNAshapeR)

# Extract sample sequences
fn <- system.file("extdata", "CGRsample.fa", package = "DNAshapeR")

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
heatShape(pred$ProT, 20)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$HelT)

