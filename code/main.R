# title: Spatial modelling of total carbon, calcium and phosphorus stocks in Amazonian Dark Earths
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

# Create variable PRETIC according to WRB specifications.
# Prepare auxiliary plotting variables
pointData$PRETIC <-
  sapply(1:nrow(pointData), function (i) {
    ifelse(pointData[i, "EXPH"] >= 30 && pointData[i, "CAMG"] >= 2 && pointData[i, "TOOC"] >= 10, 1, 0)
  })
pointData$pch <- ifelse(pointData$PRETIC == 1, 21, 1)
pointData$col <- ifelse(pointData$PRETIC == 1, "red", "ivory")

# Save point soil data
save(pointData, file = "data/R/pointData.rda")
# load("data/R/pointData.rda")

# Check spatial distribution
sp::coordinates(pointData) <- ~ x + y
sp::proj4string(pointData) <- sp::CRS("+init=epsg:32720")
plot(pointData@coords, pch = 20, cex = 0.3)

# Compute the average distance between nearest neighbouring observations
x <- pointData@coords[seq(1, nrow(pointData@coords), 5), ]
d <- apply(as.matrix(dist(x)), 1, function (x) min(x[x > 0]))
mean(d)
range(d)

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

# Soil bulk density ###########################################################################################
# Soil density data is available at three sampling locations, more specificaly at three soil profiles 
# described in the year of 2011.
# Bulk density data is in grams per cubic centimetre.
# The approach we employ consists of fitting a spline function to the soil profile data and predicting the
# soil bulk density at the five standard depths (10, 30, 50, 70, and 90 cm). Predicted values and prediction
# error variance are used to compute the soil mass per square metre in each sampling layer (20 cm depth). 
# Soil mass data is in megagrams.

# load data
density <- read.table("data/soil/density.csv", sep = ";", header = TRUE, dec = ",")
density$BUDE <- round(density$BUDE, 2)
density <- density[order(density$profile), ]
head(density)

# fit linear model to bulk density using depth as explanatory variable
# We choose the number of degrees of freedom of the natural spline based on the RMSE and R-squared 
# returned by a leave-one-out cross validation. We find out that a compromise solution is to use five 
# degrees of freedom.
df <- 1:6
# fit_bude <- lapply(df, function (i)
#   caret::train(BUDE ~ splines::ns(depth, df = i), data = density, method = "lm",
#                trControl = caret::trainControl(method = "LOOCV")))
fit_bude <- list()
fit_bude[[1]] <- caret::train(BUDE ~ splines::ns(depth, df = 1), data = density, method = "lm",
                            trControl = caret::trainControl(method = "LOOCV"))
fit_bude[[2]] <- caret::train(BUDE ~ splines::ns(depth, df = 2), data = density, method = "lm",
                              trControl = caret::trainControl(method = "LOOCV"))
fit_bude[[3]] <- caret::train(BUDE ~ splines::ns(depth, df = 3), data = density, method = "lm",
                              trControl = caret::trainControl(method = "LOOCV"))
fit_bude[[4]] <- caret::train(BUDE ~ splines::ns(depth, df = 4), data = density, method = "lm",
                              trControl = caret::trainControl(method = "LOOCV"))
fit_bude[[5]] <- caret::train(BUDE ~ splines::ns(depth, df = 5), data = density, method = "lm",
                              trControl = caret::trainControl(method = "LOOCV"))
fit_bude[[6]] <- caret::train(BUDE ~ splines::ns(depth, df = 6), data = density, method = "lm",
                              trControl = caret::trainControl(method = "LOOCV"))
cv <- sapply(fit_bude, function (x) x$results)
c(which.min(cv["RMSE", ]), which.max(cv["Rsquared", ]))
fit_bude <- lm(BUDE ~ splines::ns(depth, df = 5), density)
# splines::ns(density$depth, df = 5)

# Make predictions
newdata <- data.frame(depth = seq(0, max(density$depth), 1))
pred_bude <- predict(fit_bude, newdata = newdata, se.fit = TRUE, interval = "prediction", level = 0.90)
pred_bude2 <- predict(fit_bude, newdata = newdata, se.fit = TRUE, interval = "confidence", level = 0.90)

# Prepare figure with BUDE predictions, confidence interval and prediction interval
p <- 
  lattice::xyplot(
    newdata$d ~ pred_bude$fit[, 1], type = "l", col = "black",
    xlab = expression(paste('Bulk density (Mg ',m^-3,')', sep = '')), ylab = "Depth (cm)",
    ylim = rev(extendrange(density$depth)),
    xlim = extendrange(c(density$BUDE, pred_bude$fit[, 2:3])),
    key = list(corner = c(1, 0.1),
               points = list(pch = c(20, unique(density$profile))), 
               text = list(c("Predictions", paste("Profile", unique(density$profile))))),
    panel = function (...) {
      lattice::panel.grid(v = -1, h = -1)
      lattice::panel.polygon(
        y = c(newdata$depth, rev(newdata$depth)), x = c(pred_bude$fit[, 2], rev(pred_bude$fit[, 3])), 
        col = "gray90", border = "gray90")
      lattice::panel.polygon(
        y = c(newdata$depth, rev(newdata$depth)), x = c(pred_bude2$fit[, 2], rev(pred_bude2$fit[, 3])), 
        col = "gray85", border = "gray85")
      lattice::panel.points(density$depth ~ density$BUDE, pch = density$profile, col = "gray25")
      lattice::panel.xyplot(...)
      lattice::panel.abline(
        h = attr(splines::ns(density$depth, df = 5), "knots"), lty = "dashed", col = "gray75")
      lattice::panel.points(
        newdata$depth[newdata$depth %in% seq(10, 90, 20)] ~ 
          pred_bude$fit[newdata$depth %in% seq(10, 90, 20), 1], pch = 20, col = "black")
    }
  )
# p$par.settings <- list(fontsize = list(text = 12))
dev.off()
png("res/fig/bude.png", width = 480 * 3, height = 480 * 3, res = 72 * 4)
p
dev.off()
rm(newdata, pred_bude2, pred_bude, p, density)

# predict BUDE at standard depths and calculate the prediction error variance
bude <- predict(fit_bude, newdata = data.frame(depth = seq(10, 90, 20)), se.fit = TRUE)
bude <-  
  data.frame(mean = bude$fit, sd = sqrt(1 + c(bude$se.fit / bude$residual.scale) ^ 2) * bude$residual.scale)

# Volume of coarse fragments ##################################################################################
p <- 
  lattice::histogram(
    pointData$FRAG, xlab = "Volume of ceramics (%)", nint = 50,
    panel = function (...) {
      lattice::panel.grid(h = -1, v = -1)
      lattice::panel.rug(..., col = "gray50")
      lattice::panel.histogram(..., col = "gray85")
    })
dev.off()
png("res/fig/ceramics.png", width = 480 * 3, height = 480 * 3, res = 72 * 4)
p
dev.off()
rm(p)

# Study area ##################################################################################################
# Prepare a figure of the study area using Google Earth imagery. The figure must show the study area and 
# the Solimões River so that the reader can have an idea of the past landscape setting, with the river bed 
# running very close to the study area.

# Get image from Google Maps
location <- c(-60.24, -3.26, -60.22, -3.253)
map <- ggmap::get_map(location, maptype = "satellite", color = "bw")

# Prepare polygon data
boundary <- raster::shapefile("data/QGIS/boundary.shp")
boundary@data <- data.frame(id = 1)
boundary <- sp::spTransform(boundary, sp::CRS("+init=epsg:4326"))
boundary@data$id <- rownames(boundary@data)
lab <- boundary@bbox
lab <- lab[, 1] + apply(lab, 1, diff) * 0.5
boundary <- ggplot2::fortify(boundary, region = "id")

# Prepare image
p <- 
  ggmap::ggmap(map) + 
  ggplot2::xlab("Longitude (°)") +
  ggplot2::ylab("Latitude (°)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black"),
                 axis.text.y = ggplot2::element_text(color = "black")) +
  ggplot2::geom_polygon(
    ggplot2::aes(x = long, y = lat), boundary, show.legend = FALSE, colour = "black", fill = NA, size = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = "Solimões River", x = -60.233, y = -3.261), size = 4) +
  ggplot2::geom_text(ggplot2::aes(label = "Caldeirão", x = lab[1], y = lab[2] + 0.0005), size = 4)

# Save image
dev.off()
png("res/fig/caldeirao.png", width = 480 * 4, height = 480 * 4, res = 72 * 4)
p
dev.off()
rm(p, location, map)

# Get regional-scale image
location <- c(-60.76, -3.76, -59.00, -2.50)
map <- ggmap::get_map(location, maptype = "terrain", color = "bw")

# Prepare regional-scale image
p <- 
  ggmap::ggmap(map) + 
  ggplot2::xlab("Longitude (°)") +
  ggplot2::ylab("Latitude (°)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black"),
                 axis.text.y = ggplot2::element_text(color = "black")) +
  ggplot2::geom_polygon(
    ggplot2::aes(x = long, y = lat), boundary, show.legend = FALSE, colour = "black", fill = NA, size = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = "Caldeirão", x = lab[1] - 0.11, y = lab[2] + 0.0005), size = 4)

# Save regional-scale image
dev.off()
png("res/fig/manaus.png", width = 480 * 4, height = 480 * 4, res = 72 * 4)
p
dev.off()
rm(p, location, map)

# Get country-scale image
location <- c(-74.5, -34.5, -29.5, 5.5)
map <- ggmap::get_map(location, maptype = "terrain", color = "bw")

# Prepare regional-scale image
p <- 
  ggmap::ggmap(map) + 
  ggplot2::xlab("Longitude (°)") +
  ggplot2::ylab("Latitude (°)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black"),
                 axis.text.y = ggplot2::element_text(color = "black")) +
  ggplot2::geom_polygon(
    ggplot2::aes(x = long, y = lat), boundary, show.legend = FALSE, colour = "black", fill = NA, size = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = "Caldeirão", x = lab[1] - 0.09, y = lab[2] + 0.0005), size = 4)

# Save country-scale image
dev.off()
png("res/fig/brasil.png", width = 480 * 4, height = 480 * 4, res = 72 * 4)
p
dev.off()

rm(p, location, map, lab, boundary)

# Deterministic component of spatial variation ################################################################
# We model the soil spatial variation using depth-wise linear models where the independed variable is the
# exponential environmental gradient.
# By fitting depth-wise regression models we already take into account the effect of the depth.

# Load covariate data
covar <- spgrass7::readRAST("past_landuse")
covar$past_landuse <- (covar$past_landuse - min(covar$past_landuse, na.rm = TRUE)) / 
  (max(covar$past_landuse, na.rm = TRUE) - min(covar$past_landuse, na.rm = TRUE))

# Save image of covariate
# We change the coordinates setting the origin to (0, 0) 
map <- covar
sp::gridded(map) <- FALSE
min <- apply(map@coords, 2, min)
map@coords[, 1] <- map@coords[, 1] - min[1]
map@coords[, 2] <- map@coords[, 2] - min[2]
map@bbox <- sp::bbox(map@coords)
sp::gridded(map) <- TRUE
pts <- pointData@coords[seq(1, nrow(pointData), 5), ]
pts[, 1] <- pts[, 1] - min[1]
pts[, 2] <- pts[, 2] - min[2]
profiles <- data.frame(x = c(807955, 808041, 808131) - min[1], y = c(9640044, 9640126, 9640150) - min[2])

# prepare figure
p <- sp::spplot(
  map, 
  # col.regions = soil.colors, 
  col.regions = gray.colors(100, start = 1, end = 0),
  colorkey = TRUE, scales = list(draw = TRUE),
  xlab = "Easting (m)", ylab = "Northing (m)",
  panel = function (...) {
    lattice::panel.grid(h = -1, v = -1)
    lattice::panel.levelplot(...)
    lattice::panel.points(pts, pch = 21, fill = "lightgray", col.symbol = "black", cex = 0.5)
    # Soil profiles with BUDE
    lattice::panel.text(profiles, labels = c("P1", "P2", "P3"), col = "black", cex = 0.75, pos = 1)
    lattice::panel.points(profiles, pch = 2, col = "black", cex = 0.75)
  })
dev.off()
png("res/fig/covar.png", width = 480 * 4, height = 480 * 4, res = 72 * 4)
p
dev.off()
rm(p, map, pts)

# Use exponential
covar$past_landuse <- exp(covar$past_landuse)

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

rm(soil_var, res_anova, tooc_anova, toca_anova, toph_anova)

# Map the pretic criteria #####################################################################################
# This first attempt is based on the use of nearest neighbour distances.
pretic_dist <- SpatialTools::dist2(sp::coordinates(covar), pointData@coords[pointData$d == 10, ])
pretic_dist <- apply(pretic_dist, 1, which.min)
pretic_dist <- pointData@data[pointData$d == 10, "PRETIC"][pretic_dist]
pretic_dist <- ifelse(is.na(covar$past_landuse), NA, pretic_dist)
pretic_dist <- cbind(sp::coordinates(covar), pretic_dist)
pretic_dist <- as.data.frame(pretic_dist)
pretic_dist$pretic_dist <- ifelse(pretic_dist$pretic_dist == 1, "Pretic", "Adjacent")
pretic_dist$pretic_dist <- as.factor(pretic_dist$pretic_dist)
sp::gridded(pretic_dist) <- ~ s1 + s2

# Save figure
p <- sp::spplot(
  pretic_dist,
  col.regions = rev(gray.colors(2)),
  colorkey = TRUE, scales = list(draw = TRUE),
  xlab = "Easting (m)", ylab = "Northing (m)",
  panel = function (...) {
    lattice::panel.grid(h = -1, v = -1)
    lattice::panel.levelplot(...)
    lattice::panel.points(
      pointData@coords[pointData$d == 10, ], pch = 21, fill = "lightgray", col.symbol = "black", cex = 0.5)
  })
dev.off()
png("res/fig/pretic-dist.png", width = 480 * 4, height = 480 * 4, res = 72 * 4)
p
dev.off()
rm(p, pretic_dist)

# Stochastic component of spatial variation ###################################################################

depth <- unique(pointData$d)
range <- c(50, 75, 100)

# TOTAL ORGANIC CARBON ----
sv <- "TOOC"
tooc_data <- prepare_soil_data(pointData = pointData, sv = sv, covar = covar, save4back = TRUE)
colnames(tooc_data@data) <- gsub("TOOC", "C", colnames(tooc_data@data))
colnames(tooc_data@data) <- gsub("10", "1", colnames(tooc_data@data))
colnames(tooc_data@data) <- gsub("30", "2", colnames(tooc_data@data))
colnames(tooc_data@data) <- gsub("50", "3", colnames(tooc_data@data))
colnames(tooc_data@data) <- gsub("70", "4", colnames(tooc_data@data))
colnames(tooc_data@data) <- gsub("90", "5", colnames(tooc_data@data))
sv <- "C"

# fit competing cross-variogram models
tooc_vario <- compute_sample_variogram(soil_data = tooc_data, sv = sv, cross = TRUE)
plot(tooc_vario$v, scales = list(relation = "same"), pch = 20, cex = 0.5)
tooc_cross <- parallel::mclapply(1:length(range), function (i)
  gstat::gstat(
    tooc_vario$g, id = paste(sv, ".", depth[1], sep = ""),
    model = gstat::vgm(psill = 0.8, model = "Exp", range = range[i], nugget = 0.2), fill.all = TRUE))
tooc_lmc <- parallel::mclapply(1:length(tooc_cross), function (i)
  gstat::fit.lmc(v = tooc_vario$v, g = tooc_cross[[i]], correct.diagonal = 1.01))

# cross-validation
tooc_cv <- parallel::mclapply(
  tooc_lmc, gstat::gstat.cv, nfold = length(unique(pointData$stake)), remove.all = TRUE, all.residuals = TRUE,
  boundaries = attr(tooc_vario$v, "boundaries"), correct.diagonal = 1.01)
tooc_cv <- round(do.call(rbind, lapply(tooc_cv, colMeans)), 4)
apply(abs(tooc_cv), 2, which.min)
apply(abs(tooc_cv), 2, which.max)
tooc_cv <- tooc_cv[1, ]
tooc_lmc <- tooc_lmc[[1]]
plot(tooc_vario$v, tooc_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)
round(slope(tooc_lmc), 4)

# prepare variogram plot
tooc_plot <-
  plot(tooc_vario$v, tooc_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, 
       col = "black", strip = lattice::strip.custom(bg = "lightgray"),  xlab = "Distance (m)", 
       ylab = "Semivariance (-)")
tmp <- tooc_plot + addGridLines
tooc_plot <- tmp + latticeExtra::as.layer(tooc_plot)
dev.off()
png("res/fig/tooc_cross_vario.png", width = 480 * 3, height = 480 * 3, res = 72 * 3)
tooc_plot
dev.off()

# save results
save(tooc_vario, tooc_lmc, file = "data/R/tooc_vario.rda")
# load("data/R/tooc_vario.rda")

# make spatial predictions
t0 <- proc.time()
tooc_pred <- predict(object = tooc_lmc, newdata = covar)
proc.time() - t0
tooc_pred <- back_transform(pred = tooc_pred, soil_data = tooc_data)
save(tooc_pred, file = "data/R/tooc_pred.rda")

# make spatial simulations
nsim <- 1000
t0 <- proc.time()
set.seed(2001)
tooc_sim <- gstat::krige(
  formula = C.1 ~ past_landuse, locations = tooc_data, newdata = covar, 
  model = tooc_lmc$model[[1]], nsim = nsim, nmax = 30, maxdist = tooc_lmc$model[[1]]$range[2], 
  debug.level = -1)
proc.time() - t0

# compute tooc limits to classify as a pretic horizon
tooc_lim <- attributes(tooc_data)[c("lambda", "sv_mean", "sv_sd")]
tooc_lim <- (car::bcPower(10, tooc_lim$lambda[1]) - tooc_lim$sv_mean[1]) / tooc_lim$sv_sd[1]

# compute the number of simulated values that are larger than the limits
tooc_ade <- rowSums(apply(tooc_sim@data, 2, function (x) x >= tooc_lim)) / nsim
tooc_ade <- data.frame(cbind(sp::coordinates(covar), tooc_ade))
sp::gridded(tooc_ade) <- ~ s1 + s2
sp::proj4string(tooc_ade) <- sp::proj4string(covar)

# save figure with pretic probability
p <- sp::spplot(
  tooc_ade, col.regions = rev(gray.colors(100)), colorkey = TRUE, scales = list(draw = FALSE),
  panel = function (...) {
    lattice::panel.grid(h = -1, v = -1)
    lattice::panel.levelplot(...)
    lattice::panel.points(
      pointData@coords[pointData$d == 10, ], pch = 21, fill = "lightgray", col.symbol = "black", cex = 0.5)
  })
dev.off()
png(filename = "res/fig/tooc-pretic.png", width = 480 * 4, height = 480 * 5, res = 72 * 5)
p
dev.off()
save(tooc_ade, file = "data/R/tooc_ade.rda")
rm(tooc_ade, tooc_sim, tooc_lim, p)
gc()

# compute stocks
# load("data/R/tooc_pred.rda")
tooc_pred@data[, seq(2, 10, 2)] <- 
  lapply(1:5, function (i) {
    0.2 * sqrt(bude$mean[i] ^ 2 * tooc_pred@data[, seq(2, 10, 2)][, i] ^ 2 + 
                 tooc_pred@data[, seq(1, 9, 2)][, i] ^ 2 * bude$sd[i] ^ 2)
  })
tooc_pred@data[, seq(1, 9, 2)] <- 
  lapply(1:5, function (i) tooc_pred@data[, seq(1, 9, 2)][, i] * bude$mean[i] * 0.2)
tooc_pred@data$stock <- rowSums(tooc_pred@data[, seq(1, 9, 2)])
tooc_pred@data$stock_var <- sqrt(rowSums(tooc_pred@data[, seq(2, 10, 2)] ^ 2))
sum(tooc_pred$stock, na.rm = TRUE)
mean(tooc_pred$stock, na.rm = TRUE)
range(tooc_pred$stock, na.rm = TRUE)

# save figure with depth-wise predictions
map <- layer_predictions(tooc_pred, "pred")
dev.off()
png(filename = "res/fig/tooc_pred.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
map
dev.off()
rm(map)

# save figure with depth-wise prediction error standard deviation
map <- layer_predictions(tooc_pred, "var")
dev.off()
png(filename = "res/fig/tooc_sd.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
map
dev.off()
rm(map)

# save figure with profile predictions
map <- profile_predictions(tooc_pred, "stock")
dev.off()
png(filename = "res/fig/tooc_pred_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
map
dev.off()
rm(map)

# save figure with profile prediction error standard deviation
map <- profile_predictions(tooc_pred, "stock_var", col.regions = uncertainty.colors)
dev.off()
png(filename = "res/fig/tooc_sd_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
map
dev.off()
rm(map)

# CAMG = EXCA + EXMG ----
sv <- "CAMG"
camg_data <- prepare_soil_data(pointData = pointData, sv = sv, covar = covar, save4back = TRUE)
# colnames(camg_data@data) <- gsub("camg", "C", colnames(camg_data@data))
# colnames(camg_data@data) <- gsub("10", "1", colnames(camg_data@data))
# colnames(camg_data@data) <- gsub("30", "2", colnames(camg_data@data))
# colnames(camg_data@data) <- gsub("50", "3", colnames(camg_data@data))
# colnames(camg_data@data) <- gsub("70", "4", colnames(camg_data@data))
# colnames(camg_data@data) <- gsub("90", "5", colnames(camg_data@data))
sv <- "CAMG"

# fit competing cross-variogram models
camg_vario <- compute_sample_variogram(soil_data = camg_data, sv = sv, cross = TRUE)
plot(camg_vario$v, scales = list(relation = "same"), pch = 20, cex = 0.5)
camg_cross <- parallel::mclapply(1:length(range), function (i)
  gstat::gstat(
    camg_vario$g, id = paste(sv, ".", depth[1], sep = ""),
    model = gstat::vgm(psill = 0.6, model = "Exp", range = range[i], nugget = 0.2), fill.all = TRUE))
camg_lmc <- parallel::mclapply(1:length(camg_cross), function (i)
  gstat::fit.lmc(v = camg_vario$v, g = camg_cross[[i]], correct.diagonal = 1.01))

# cross-validation
camg_cv <- parallel::mclapply(
  camg_lmc, gstat::gstat.cv, nfold = length(unique(pointData$stake)), remove.all = TRUE, all.residuals = TRUE,
  boundaries = attr(camg_vario$v, "boundaries"), correct.diagonal = 1.01)
camg_cv <- round(do.call(rbind, lapply(camg_cv, colMeans)), 4)
apply(abs(camg_cv), 2, which.min)
apply(abs(camg_cv), 2, which.max)
camg_cv <- camg_cv[3, ]
camg_lmc <- camg_lmc[[3]]
plot(camg_vario$v, camg_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)

# prepare variogram plot
# camg_plot <-
#   plot(camg_vario$v, camg_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, 
#        col = "black", strip = lattice::strip.custom(bg = "lightgray"),  xlab = "Distance (m)", 
#        ylab = "Semivariance (-)")
# tmp <- camg_plot + addGridLines
# camg_plot <- tmp + latticeExtra::as.layer(camg_plot)
# dev.off()
# png("res/fig/camg_cross_vario.png", width = 480 * 3, height = 480 * 3, res = 72 * 3)
# camg_plot
# dev.off()

# save results
# save(camg_vario, camg_lmc, file = "data/R/camg_vario.rda")
# load("data/R/camg_vario.rda")

# make spatial predictions
# t0 <- proc.time()
# camg_pred <- predict(object = camg_lmc, newdata = covar)
# proc.time() - t0
# camg_pred <- back_transform(pred = camg_pred, soil_data = camg_data)
# save(camg_pred, file = "data/R/camg_pred.rda")

# make spatial simulations
nsim <- 1000
t0 <- proc.time()
set.seed(2001)
camg_sim <- gstat::krige(
  formula = CAMG.10 ~ past_landuse, locations = camg_data, newdata = covar, 
  model = camg_lmc$model[[1]], nsim = nsim, nmax = 30, maxdist = camg_lmc$model[[1]]$range[2], 
  debug.level = -1)
proc.time() - t0

# compute camg limits to classify as a pretic horizon
camg_lim <- attributes(camg_data)[c("lambda", "sv_mean", "sv_sd")]
camg_lim <- (car::bcPower(2, camg_lim$lambda[1]) - camg_lim$sv_mean[1]) / camg_lim$sv_sd[1]

# compute the number of simulated values that are larger than the limits
camg_ade <- rowSums(apply(camg_sim@data, 2, function (x) x >= camg_lim)) / nsim
rm(camg_sim); gc()
camg_ade <- data.frame(cbind(sp::coordinates(covar), camg_ade))
sp::gridded(camg_ade) <- ~ s1 + s2
sp::proj4string(camg_ade) <- sp::proj4string(covar)

# save figure with pretic probability
p <- sp::spplot(
  camg_ade, col.regions = rev(gray.colors(100)), colorkey = TRUE, scales = list(draw = FALSE),
  panel = function (...) {
    lattice::panel.grid(h = -1, v = -1)
    lattice::panel.levelplot(...)
    lattice::panel.points(
      pointData@coords[pointData$d == 10, ], pch = 21, fill = "lightgray", col.symbol = "black", cex = 0.5)
  })
dev.off()
png(filename = "res/fig/camg-pretic.png", width = 480 * 4.1, height = 480 * 5, res = 72 * 5)
p
dev.off()
save(camg_ade, file = "data/R/camg_ade.rda")
rm(camg_ade, camg_lim, p)
gc()

# # compute stocks
# # load("data/R/camg_pred.rda")
# camg_pred@data[, seq(2, 10, 2)] <- 
#   lapply(1:5, function (i) {
#     0.2 * sqrt(bude$mean[i]^2 * camg_pred@data[, seq(2, 10, 2)][, i]^2 + 
#                  camg_pred@data[, seq(1, 9, 2)][, i]^2 * bude$sd[i]^2)
#   })
# camg_pred@data[, seq(1, 9, 2)] <- 
#   lapply(1:5, function (i) camg_pred@data[, seq(1, 9, 2)][, i] * bude$mean[i] * 0.2)
# # camg_pred@data[, seq(2, 10, 2)] <- camg_pred@data[, seq(2, 10, 2)] / camg_pred@data[, seq(1, 9, 2)]
# camg_pred@data$stock <- rowSums(camg_pred@data[, seq(1, 9, 2)])
# camg_pred@data$stock_var <- sqrt(rowSums(camg_pred@data[, seq(2, 10, 2)]^2))
# totals(camg_pred)
# range(camg_pred$stock, na.rm = TRUE)

# # save figure with depth-wise predictions
# map <- layer_predictions(camg_pred, "pred")
# dev.off()
# png(filename = "res/fig/camg_pred.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
# map
# dev.off()
# rm(map)

# # save figure with depth-wise prediction error standard deviation
# map <- layer_predictions(camg_pred, "var")
# dev.off()
# png(filename = "res/fig/camg_sd.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
# map
# dev.off()
# rm(map)

# # save figure with profile predictions
# map <- profile_predictions(camg_pred, "stock")
# dev.off()
# png(filename = "res/fig/camg_pred_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
# map
# dev.off()
# rm(map)

# # save figure with profile prediction error standard deviation
# map <- profile_predictions(camg_pred, "stock_var", col.regions = uncertainty.colors)
# dev.off()
# png(filename = "res/fig/camg_sd_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
# map
# dev.off()
# rm(map)

# Extratable phosphorus ----
sv <- "EXPH"
exph_data <- prepare_soil_data(pointData = pointData, sv = sv, covar = covar, save4back = TRUE)
# colnames(exph_data@data) <- gsub("exph", "C", colnames(exph_data@data))
# colnames(exph_data@data) <- gsub("10", "1", colnames(exph_data@data))
# colnames(exph_data@data) <- gsub("30", "2", colnames(exph_data@data))
# colnames(exph_data@data) <- gsub("50", "3", colnames(exph_data@data))
# colnames(exph_data@data) <- gsub("70", "4", colnames(exph_data@data))
# colnames(exph_data@data) <- gsub("90", "5", colnames(exph_data@data))
sv <- "exph"

# fit competing cross-variogram models
exph_vario <- compute_sample_variogram(soil_data = exph_data, sv = sv, cross = TRUE)
plot(exph_vario$v, scales = list(relation = "same"), pch = 20, cex = 0.5)
exph_cross <- parallel::mclapply(1:length(range), function (i)
  gstat::gstat(
    exph_vario$g, id = paste(sv, ".", depth[1], sep = ""),
    model = gstat::vgm(psill = 0.6, model = "Exp", range = range[i], nugget = 0.2), fill.all = TRUE))
exph_lmc <- parallel::mclapply(1:length(exph_cross), function (i)
  gstat::fit.lmc(v = exph_vario$v, g = exph_cross[[i]], correct.diagonal = 1.01))

# cross-validation
exph_cv <- parallel::mclapply(
  exph_lmc, gstat::gstat.cv, nfold = length(unique(pointData$stake)), remove.all = TRUE, all.residuals = TRUE,
  boundaries = attr(exph_vario$v, "boundaries"), correct.diagonal = 1.01)
exph_cv <- round(do.call(rbind, lapply(exph_cv, colMeans)), 4)
apply(abs(exph_cv), 2, which.min)
apply(abs(exph_cv), 2, which.max)
exph_cv <- exph_cv[1, ]
exph_lmc <- exph_lmc[[1]]
plot(exph_vario$v, exph_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)

# prepare variogram plot
# exph_plot <-
#   plot(exph_vario$v, exph_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, 
#        col = "black", strip = lattice::strip.custom(bg = "lightgray"),  xlab = "Distance (m)", 
#        ylab = "Semivariance (-)")
# tmp <- exph_plot + addGridLines
# exph_plot <- tmp + latticeExtra::as.layer(exph_plot)
# dev.off()
# png("res/fig/exph_cross_vario.png", width = 480 * 3, height = 480 * 3, res = 72 * 3)
# exph_plot
# dev.off()

# save results
# save(exph_vario, exph_lmc, file = "data/R/exph_vario.rda")
# load("data/R/exph_vario.rda")

# make spatial predictions
# t0 <- proc.time()
# exph_pred <- predict(object = exph_lmc, newdata = covar)
# proc.time() - t0
# exph_pred <- back_transform(pred = exph_pred, soil_data = exph_data)
# save(exph_pred, file = "data/R/exph_pred.rda")

# make spatial simulations
nsim <- 1000
t0 <- proc.time()
set.seed(2001)
exph_sim <- gstat::krige(
  formula = EXPH.10 ~ past_landuse, locations = exph_data, newdata = covar, 
  model = exph_lmc$model[[1]], nsim = nsim, nmax = 30, maxdist = exph_lmc$model[[1]]$range[2], 
  debug.level = -1)
proc.time() - t0

# compute exph limits to classify as a pretic horizon
# set limit accordingly!!!
exph_lim <- attributes(exph_data)[c("lambda", "sv_mean", "sv_sd")]
exph_lim <- (car::bcPower(30, exph_lim$lambda[1]) - exph_lim$sv_mean[1]) / exph_lim$sv_sd[1]

# compute the number of simulated values that are larger than the limits
exph_ade <- rowSums(apply(exph_sim@data, 2, function (x) x >= exph_lim)) / nsim
rm(exph_sim); gc()
exph_ade <- data.frame(cbind(sp::coordinates(covar), exph_ade))
sp::gridded(exph_ade) <- ~ s1 + s2
sp::proj4string(exph_ade) <- sp::proj4string(covar)

# save figure with pretic probability
p <- sp::spplot(
  exph_ade, col.regions = rev(gray.colors(100)), colorkey = TRUE, scales = list(draw = FALSE),
  panel = function (...) {
    lattice::panel.grid(h = -1, v = -1)
    lattice::panel.levelplot(...)
    lattice::panel.points(
      pointData@coords[pointData$d == 10, ], pch = 21, fill = "lightgray", col.symbol = "black", cex = 0.5)
  })
dev.off()
png(filename = "res/fig/exph-pretic.png", width = 480 * 4.1, height = 480 * 5, res = 72 * 5)
p
dev.off()
save(exph_ade, file = "data/R/exph_ade.rda")
rm(exph_ade, exph_lim, p)
gc()

# # compute stocks
# # load("data/R/exph_pred.rda")
# exph_pred@data[, seq(2, 10, 2)] <- 
#   lapply(1:5, function (i) {
#     0.2 * sqrt(bude$mean[i]^2 * exph_pred@data[, seq(2, 10, 2)][, i]^2 + 
#                  exph_pred@data[, seq(1, 9, 2)][, i]^2 * bude$sd[i]^2)
#   })
# exph_pred@data[, seq(1, 9, 2)] <- 
#   lapply(1:5, function (i) exph_pred@data[, seq(1, 9, 2)][, i] * bude$mean[i] * 0.2)
# # exph_pred@data[, seq(2, 10, 2)] <- exph_pred@data[, seq(2, 10, 2)] / exph_pred@data[, seq(1, 9, 2)]
# exph_pred@data$stock <- rowSums(exph_pred@data[, seq(1, 9, 2)])
# exph_pred@data$stock_var <- sqrt(rowSums(exph_pred@data[, seq(2, 10, 2)]^2))
# totals(exph_pred)
# range(exph_pred$stock, na.rm = TRUE)

# # save figure with depth-wise predictions
# map <- layer_predictions(exph_pred, "pred")
# dev.off()
# png(filename = "res/fig/exph_pred.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
# map
# dev.off()
# rm(map)

# # save figure with depth-wise prediction error standard deviation
# map <- layer_predictions(exph_pred, "var")
# dev.off()
# png(filename = "res/fig/exph_sd.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
# map
# dev.off()
# rm(map)

# # save figure with profile predictions
# map <- profile_predictions(exph_pred, "stock")
# dev.off()
# png(filename = "res/fig/exph_pred_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
# map
# dev.off()
# rm(map)

# # save figure with profile prediction error standard deviation
# map <- profile_predictions(exph_pred, "stock_var", col.regions = uncertainty.colors)
# dev.off()
# png(filename = "res/fig/exph_sd_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
# map
# dev.off()
# rm(map)

# PRETIC ----

# process data
load("data/R/tooc_ade.rda")
load("data/R/camg_ade.rda")
load("data/R/exph_ade.rda")
pretic <- raster::brick(list(raster::raster(tooc_ade), raster::raster(camg_ade), raster::raster(exph_ade)))
pretic <- min(pretic)
pretic <- as(pretic, "SpatialPixelsDataFrame")
rm(tooc_ade, camg_ade, exph_ade)
gc()
p <- 0.9
pretic$pretic <- as.factor(ifelse(pretic$layer >= p, "Pretic", "Adjacent"))
summary(pretic$pretic)

# prepare figure
min <- apply(pretic@coords, 2, min)
pretic@coords[, 1] <- pretic@coords[, 1] - min[1]
pretic@coords[, 2] <- pretic@coords[, 2] - min[2]
pretic@bbox <- sp::bbox(pretic@coords)
pts <- pointData@coords[seq(1, nrow(pointData), 5), ]
pts[, 1] <- pts[, 1] - min[1]
pts[, 2] <- pts[, 2] - min[2]
at <- seq(0, 1, by = 0.1)
col <- gray.colors(length(at) - 1, start = 1, end = 0)
p <- sp::spplot(
  pretic, 1, at = at, col.regions = col, colorkey = TRUE, scales = list(draw = TRUE), 
  xlab = "Easting (m)", ylab = "Northing (m)",
  panel = function (...) {
    lattice::panel.grid(h = -1, v = -1)
    lattice::panel.levelplot(...)
  }) + 
  # latticeExtra::as.layer(lattice::contourplot(
    # layer ~ x + y, as.data.frame(pretic), at = c(0.25, 0.50, 0.75), col = "black", labels = FALSE)) +
  # latticeExtra::layer(
    # lattice::panel.text(x = c(208, 262), y = c(350, 250), labels = "0.50", cex = 0.75)) +
  latticeExtra::layer(
    lattice::panel.points(pts, pch = 21, fill = "lightgray", col.symbol = "black", cex = 0.5))

# save figure
dev.off()
png(filename = "res/fig/pretic-prob.png", width = 480 * 4.1, height = 480 * 5, res = 72 * 5)
p
dev.off()
rm(p, pts)
gc()

# Stats for carbon stocks
pretic_stats(pretic, tooc_pred)

# TOTAL CALCIUM ----
sv <- "TOCA"
toca_data <- prepare_soil_data(pointData = pointData, sv = sv, covar = covar, save4back = TRUE)
colnames(toca_data@data) <- gsub("TOCA", "Ca", colnames(toca_data@data))
colnames(toca_data@data) <- gsub("10", "1", colnames(toca_data@data))
colnames(toca_data@data) <- gsub("30", "2", colnames(toca_data@data))
colnames(toca_data@data) <- gsub("50", "3", colnames(toca_data@data))
colnames(toca_data@data) <- gsub("70", "4", colnames(toca_data@data))
colnames(toca_data@data) <- gsub("90", "5", colnames(toca_data@data))
sv <- "Ca"

# fit competing cross-variogram models
toca_vario <- compute_sample_variogram(soil_data = toca_data, sv = sv, cross = TRUE)
plot(toca_vario$v, scales = list(relation = "same"), pch = 20, cex = 0.5)
toca_cross <- parallel::mclapply(1:length(range), function (i)
  gstat::gstat(
    toca_vario$g, id = paste(sv, ".", depth[1], sep = ""),
    model = gstat::vgm(psill = 0.8, model = "Exp", range = range[i], nugget = 0.2), fill.all = TRUE))
toca_lmc <- parallel::mclapply(1:length(toca_cross), function (i)
  gstat::fit.lmc(v = toca_vario$v, g = toca_cross[[i]], correct.diagonal = 1.01))

# cross-validation
toca_cv <- parallel::mclapply(
  toca_lmc, gstat::gstat.cv, nfold = length(unique(pointData$stake)), remove.all = TRUE, all.residuals = TRUE,
  boundaries = attr(toca_vario$v, "boundaries"), correct.diagonal = 1.01)
toca_cv <- round(do.call(rbind, lapply(toca_cv, colMeans)), 4)
apply(abs(toca_cv), 2, which.min)
apply(abs(toca_cv), 2, which.max)
toca_cv <- toca_cv[3, ]
toca_lmc <- toca_lmc[[3]]
plot(toca_vario$v, toca_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)
round(slope(toca_lmc), 4)

# prepare variogram plot
toca_plot <- 
  plot(toca_vario$v, toca_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, col = "black", 
       strip = lattice::strip.custom(bg = "lightgray"), xlab = "Distance (m)", ylab = "Semivariance")
tmp <- toca_plot + addGridLines
toca_plot <- tmp + latticeExtra::as.layer(toca_plot)
dev.off()
png("res/fig/toca_cross_vario.png", width = 480 * 3, height = 480 * 3, res = 72 * 3)
toca_plot
dev.off()

# save results
save(toca_vario, toca_lmc, file = "data/R/toca_vario.rda")

# make spatial predictions
t0 <- proc.time()
toca_pred <- predict(object = toca_lmc, newdata = covar)
proc.time() - t0
toca_pred <- back_transform(pred = toca_pred, soil_data = toca_data)
save(toca_pred, file = "data/R/toca_pred.rda")

# compute stocks
# load("data/R/toca_pred.rda")
toca_pred@data[, seq(2, 10, 2)] <- 
  lapply(1:5, function (i) {
    0.2 * sqrt(bude$mean[i]^2 * toca_pred@data[, seq(2, 10, 2)][, i]^2 + 
                 toca_pred@data[, seq(1, 9, 2)][, i]^2 * bude$sd[i]^2)
  })
toca_pred@data[, seq(1, 9, 2)] <- 
  lapply(1:5, function (i) toca_pred@data[, seq(1, 9, 2)][, i] * bude$mean[i] * 0.2)
toca_pred@data$stock <- rowSums(toca_pred@data[, seq(1, 9, 2)])
toca_pred@data$stock_var <- sqrt(rowSums(toca_pred@data[, seq(2, 10, 2)] ^ 2))

sum(toca_pred$stock, na.rm = TRUE)
mean(toca_pred$stock, na.rm = TRUE)
range(toca_pred$stock, na.rm = TRUE)

# Stats for stocks
pretic_stats(pretic, toca_pred)

# save figure with depth-wise predictions
map <- layer_predictions(toca_pred, "pred")
dev.off()
png(filename = "res/fig/toca_pred.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
map
dev.off()
rm(map)

# save figure with depth-wise prediction error standard deviation
map <- layer_predictions(toca_pred, "var")
dev.off()
png(filename = "res/fig/toca_sd.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
map
dev.off()
rm(map)

# save figure with profile predictions
map <- profile_predictions(toca_pred, "stock")
dev.off()
png(filename = "res/fig/toca_pred_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
map
dev.off()
rm(map)

# save figure with profile prediction error standard deviation
map <- profile_predictions(toca_pred, "stock_var", col.regions = uncertainty.colors)
dev.off()
png(filename = "res/fig/toca_sd_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
map
dev.off()
rm(map)

# TOTAL PHOSPHORUS ----
sv <- "TOPH"
toph_data <- prepare_soil_data(pointData = pointData, sv = sv, covar = covar, save4back = TRUE)
colnames(toph_data@data) <- gsub("TOPH", "P", colnames(toph_data@data))
colnames(toph_data@data) <- gsub("10", "1", colnames(toph_data@data))
colnames(toph_data@data) <- gsub("30", "2", colnames(toph_data@data))
colnames(toph_data@data) <- gsub("50", "3", colnames(toph_data@data))
colnames(toph_data@data) <- gsub("70", "4", colnames(toph_data@data))
colnames(toph_data@data) <- gsub("90", "5", colnames(toph_data@data))
sv <- "P"

# fit competing cross-variogram models
toph_vario <- compute_sample_variogram(soil_data = toph_data, sv = sv, cross = TRUE)
plot(toph_vario$v, scales = list(relation = "same"), pch = 20, cex = 0.5)
toph_cross <- parallel::mclapply(1:length(range), function (i)
  gstat::gstat(
    toph_vario$g, id = paste(sv, ".", depth[1], sep = ""),
    model = gstat::vgm(psill = 0.8, model = "Exp", range = range[i], nugget = 0.2), fill.all = TRUE))
toph_lmc <- parallel::mclapply(1:length(toph_cross), function (i)
  gstat::fit.lmc(v = toph_vario$v, g = toph_cross[[i]], correct.diagonal = 1.01))

# cross-validation
toph_cv <- parallel::mclapply(
  toph_lmc, gstat::gstat.cv, nfold = length(unique(pointData$stake)), remove.all = TRUE, all.residuals = TRUE,
  boundaries = attr(toph_vario$v, "boundaries"), correct.diagonal = 1.01)
toph_cv <- round(do.call(rbind, lapply(toph_cv, colMeans)), 4)
apply(abs(toph_cv), 2, which.min)
apply(abs(toph_cv), 2, which.max)
toph_cv <- toph_cv[3, ]
toph_lmc <- toph_lmc[[3]]
plot(toph_vario$v, toph_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)
round(slope(toph_lmc), 4)

# prepare variogram plot
toph_plot <- 
  plot(toph_vario$v, toph_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, col = "black", 
       strip = lattice::strip.custom(bg = "lightgray"), xlab = "Distance (m)", ylab = "Semivariance")
tmp <- toph_plot + addGridLines
toph_plot <- tmp + latticeExtra::as.layer(toph_plot)
dev.off()
png("res/fig/toph_cross_vario.png", width = 480 * 3, height = 480 * 3, res = 72 * 3)
toph_plot
dev.off()

# save results
save(toph_vario, toph_lmc, file = "data/R/toph_vario.rda")
# load("data/R/toph_vario.rda")

# make spatial predictions
t0 <- proc.time()
toph_pred <- predict(object = toph_lmc, newdata = covar)
proc.time() - t0
toph_pred <- back_transform(pred = toph_pred, soil_data = toph_data)
save(toph_pred, file = "data/R/toph_pred.rda")

# compute stocks
# load("data/R/toph_pred.rda")
toph_pred@data[, seq(2, 10, 2)] <- 
  lapply(1:5, function (i) {
    0.2 * sqrt(bude$mean[i]^2 * toph_pred@data[, seq(2, 10, 2)][, i]^2 + 
                 toph_pred@data[, seq(1, 9, 2)][, i]^2 * bude$sd[i]^2) 
  })
toph_pred@data[, seq(1, 9, 2)] <- 
  lapply(1:5, function (i) toph_pred@data[, seq(1, 9, 2)][, i] * bude$mean[i] * 0.2)
toph_pred@data$stock <- rowSums(toph_pred@data[, seq(1, 9, 2)])
toph_pred@data$stock_var <- sqrt(rowSums(toph_pred@data[, seq(2, 10, 2)]^2))

sum(toph_pred$stock, na.rm = TRUE)
mean(toph_pred$stock, na.rm = TRUE)
range(toph_pred$stock, na.rm = TRUE)

# Stats for stocks
pretic_stats(pretic, toph_pred)

# save figure with depth-wise predictions
map <- layer_predictions(toph_pred, "pred")
dev.off()
png(filename = "res/fig/toph_pred.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
map
dev.off()
rm(map)

# save figure with depth-wise prediction error standard deviation
map <- layer_predictions(toph_pred, "var")
dev.off()
png(filename = "res/fig/toph_sd.png", height = 480 * 1.4, width = 480 * 4, res = 72 * 3)
map
dev.off()
rm(map)

# save figure with profile predictions
map <- profile_predictions(toph_pred, "stock", col.regions = soil.colors)
dev.off()
png(filename = "res/fig/toph_pred_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
map
dev.off()
rm(map)

# save figure with profile prediction error standard deviation
map <- profile_predictions(toph_pred, "stock_var", col.regions = uncertainty.colors)
dev.off()
png(filename = "res/fig/toph_sd_profile.png", width = 480 * 4, height = 480 * 5, res = 72 * 6)
map
dev.off()
rm(map)
