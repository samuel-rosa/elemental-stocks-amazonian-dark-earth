# title: Three-dimensional geostatistical modelling of elemental stocks in Amazonian Dark Earths
# author: Alessandro Samuel-Rosa

# Initial settings ############################################################################################

# Clean up and load user defined functions
rm(list = ls())
source("code/helper.R")

# Start GRASS
spgrass7::initGRASS(
  gisBase = gisBase, gisDbase = "data/GRASS", location = "caldeirao", mapset = "predictions", 
  override = TRUE, pid = Sys.getpid())

# Load and pre-process point soil data
# The coordinates of the sample points are given by `d` (depth), `x` (longitude), and `y` (latitude). The 
# depth `d` of the sample is equal to the center of the sampling layer. The identification of each sample 
# point is given by the variable `stake`.
pointData <- read.table("data/soil/caldeirao.csv", sep = "\t", header = T)
head(pointData)
length(unique(pointData$stake)) # number of unique points (n = 53)
dim(pointData)

# Round values
pointData$EXPH <- round(pointData$EXPH, 0)
pointData$EXCA <- round(pointData$EXCA, 1)
pointData$EXMG <- round(pointData$EXMG, 1)
pointData$TOOC <- round(pointData$TOOC, 1)
pointData$TOPH <- round(pointData$TOPH, 1)
pointData$TOCA <- round(pointData$TOCA, 1)

# Create variable CAMG = EXCA + EXMG
pointData$CAMG <- pointData$EXCA + pointData$EXMG

# Save point soil data
save(pointData, file = "data/R/pointData.rda")

# Check spatial distribution
sp::coordinates(pointData) <- ~ x + y
sp::proj4string(pointData) <- sp::CRS("+init=epsg:32720")
sp::plot(pointData, pch = 20, cex = 0.3)

# Check depth-wise empirical distribution
# The bw plots show that the concentration of both `C` and `Ca` decreases with depth, as well as its 
# variation. This suggests that there is a strong vertical trend in `C` and `Ca` concentrations. On the other
# hand, both the concentration and variation of `P` present a slight increase with depth. The bw plots also
# suggest that there is an increase in the asymmetry of the data with depth. This means that the data requires
# transformation before modelling. We choose the Box-Cox family of power transformations to make the data 
# empirical distribution closer to the normal. A specific Box-Cox transform is used for each depth.
p <- depth_bwplot(pts = pointData, layout = c(3, 1))
p$index.cond[[1]] <- c(2, 1, 3)
p$condlevels$ind <- c("Ca", "TOC", "P")
dev.off()
png(paste("res/fig/box-and-whisker.png", sep = ""), width = 1200, height = 480, res = 150)
p
dev.off()
rm(p)

