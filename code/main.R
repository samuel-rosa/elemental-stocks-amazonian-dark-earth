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
pointData <- read.table("data/soil/caldeirao.csv", sep = ";", header = T)
head(pointData)
length(unique(pointData$stake)) # number of unique points (n = 53)
dim(pointData)

# Round values
pointData$PH <- round(pointData$PH, 1)
pointData$CLAY <- round(pointData$CLAY, 0)
pointData$ECEC <- round(pointData$ECEC, 1)
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
factor.levels <- 
  c(expression(paste('Clay (g ',kg^-1,')', sep = '')),
    expression(paste('CEC (', cmol[c], " ", kg^-1,')', sep = '')),
    "pH",
    expression(paste('Ca (g ',kg^-1,')', sep = '')),
    expression(paste('C (g ',kg^-1,')', sep = '')),
    expression(paste('P (g ',kg^-1,')', sep = '')))
p <- depth_bwplot(
  pts = pointData, vars = c("TOOC", "TOCA", "TOPH", "PH", "CLAY", "ECEC"), layout = c(3, 2),
  strip = lattice::strip.custom(bg = "lightgray", factor.levels = factor.levels))
p$index.cond[[1]] <- c(2, 3, 1, 5, 4, 6)
dev.off()
png(paste("res/fig/box-and-whisker.png", sep = ""), width = 1200, height = 1200 / 1.6, res = 150)
p
dev.off()
rm(p)

# Deterministic component of spatial variation ################################################################
# We model the soil spatial variation using depth-wise linear models where the independed variable is the
# exponential environmental gradient.
# By fitting depth-wise regression models we already take into account the effect of the depth.

# Load covariate data
covar <- spgrass7::readRAST("past_landuse")
covar$past_landuse <- (covar$past_landuse - min(covar$past_landuse, na.rm = TRUE)) / 
  (max(covar$past_landuse, na.rm = TRUE) - min(covar$past_landuse, na.rm = TRUE))
covar$past_landuse <- exp(covar$past_landuse)
sp::spplot(covar, col.regions = soil.colors)

# Total organic carbon
soil_var <- prepare_soil_data(pointData = pointData, sv = "TOOC", covar = covar)
tooc_anova <- get_lm_output(soil_var = soil_var, sv = "TOOC")

# Total calcium
soil_var <- prepare_soil_data(pointData = pointData, sv = "TOCA", covar = covar)
toca_anova <- get_lm_output(soil_var = soil_var, sv = "TOCA")

# Total phosphorus
soil_var <- prepare_soil_data(pointData = pointData, sv = "TOPH", covar = covar)
toph_anova <- get_lm_output(soil_var = soil_var, sv = "TOPH")

# Save anova and lm results in a csv file
res_anova <- rbind(cbind("tooc",  tooc_anova), cbind("toca", toca_anova), cbind("toph", toph_anova))
write.csv(res_anova, file = "res/tab/anova.csv")
