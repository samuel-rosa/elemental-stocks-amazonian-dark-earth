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

# Create variable PRETIC according to WRB specifications
pointData$PRETIC <-
  apply(pointData, 1, function (x) 
    ifelse(x["EXPH"] >= 30 && x["CAMG"] >= 2 && x["TOOC"] >= 10, 1, 0))

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

rm(soil_var, res_anova, tooc_anova, toca_anova, toph_anova)

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
tooc_cv <- tooc_cv[[1]]
tooc_lmc <- tooc_lmc[[1]]
plot(tooc_vario$v, tooc_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)

# prepare variogram plot
tooc_plot <-
  plot(tooc_vario$v, tooc_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, 
       col = "black", strip = lattice::strip.custom(bg = "lightgray"),  xlab = "Distance (m)", 
       ylab = "Semivariance (-)")
tmp <- tooc_plot + addGridLines
tooc_plot <- tmp + latticeExtra::as.layer(tooc_plot)
dev.off()
png("res/fig/tooc_cross_vario.png", width = 480 * 2, height = 480 * 2, res = 72 * 2)
tooc_plot
dev.off()

# save results
save(tooc_vario, tooc_lmc, file = "data/R/tooc_vario.rda")

# make spatial predictions
t0 <- proc.time()
tooc_pred <- predict(object = tooc_lmc, newdata = covar)
proc.time() - t0
# tooc_pred <- back_transform(pred = tooc_pred, soil_data = tooc_data, depth = depth, n.sim = 10000)
save(tooc_pred, file = "data/R/tooc_pred.rda")

# make spatial simulations
nsim <- 1
pts <- 50000
simgrid <- make_simgrid(grid = covar, n = pts)
simgrid <- simgrid[sample(1:length(simgrid), size = length(simgrid)), ]
n <- round(seq(1, length(simgrid), length.out = pts / 250))
pts <- cbind(n[-length(n)], n[-1] + 1)
pts[nrow(pts), ncol(pts)] <- length(simgrid)
simgrid <- lapply(1:nrow(pts), function (i) simgrid[pts[i, 1]:pts[i, 2], ])
str(simgrid, 1)
t0 <- proc.time()
set.seed(1984)
tooc_sim <- parallel::mclapply(simgrid, function (x)
  predict(object = tooc_lmc, newdata = x, nsim = nsim, debug.level = -1, 
          nmax = 12, maxdist = tooc_lmc$model$C.1$range[2]))
# tooc_sim <- predict(object = tooc_lmc, newdata = simgrid, nsim = nsim, debug.level = -1, nmax = 4)
(proc.time() - t0) / 3600

tooc_sim <- do.call(rbind, tooc_sim)
simgrid <- do.call(rbind, simgrid)
str(tooc_sim)
str(simgrid)

sp::gridded(tooc_sim) <- TRUE
sp::gridded(simgrid) <- TRUE

sp::spplot(tooc_sim, col.regions = soil.colors(100))
sp::spplot(simgrid, col.regions = soil.colors(100))

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
toca_cv <- toca_cv[[3]]
toca_lmc <- toca_lmc[[3]]
plot(toca_vario$v, toca_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)

# prepare variogram plot
toca_plot <- 
  plot(toca_vario$v, toca_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, col = "black", 
       strip = lattice::strip.custom(bg = "lightgray"), xlab = "Distance (m)", ylab = "Semivariance")
dev.off()
png("res/fig/toca_cross_vario.png", width = 480 * 2, height = 480 * 2, res = 72 * 2)
toca_plot
dev.off()

# save results
save(toca_vario, toca_lmc, file = "data/R/toca_vario.rda")

# make spatial predictions
t0 <- proc.time()
toca_pred <- predict(object = toca_lmc, newdata = covar)
proc.time() - t0
# toca_pred <- back_transform(pred = toca_pred, soil_data = toca_data, depth = depth, n.sim = 10000)
save(toca_pred, file = "data/R/toca_pred.rda")

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
toph_cv <- toph_cv[[3]]
toph_lmc <- toph_lmc[[3]]
plot(toph_vario$v, toph_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5)

# prepare variogram plot
toph_plot <- 
  plot(toph_vario$v, toph_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, col = "black", 
       strip = lattice::strip.custom(bg = "lightgray"), xlab = "Distance (m)", ylab = "Semivariance")
dev.off()
png("res/fig/toph_cross_vario.png", width = 480 * 2, height = 480 * 2, res = 72 * 2)
toph_plot
dev.off()

# save results
save(toph_vario, toph_lmc, file = "data/R/toph_vario.rda")

# make spatial predictions
t0 <- proc.time()
toph_pred <- predict(object = toph_lmc, newdata = covar)
proc.time() - t0
# toph_pred <- back_transform(pred = toph_pred, soil_data = toph_data, depth = depth, n.sim = 10000)
save(toph_pred, file = "data/R/toph_pred.rda")

# Identification of ADE according to pretic criteria ##########################################################

# built linear model of coregionalization for Ca+Mg, P and C
pretic_data <- data.frame(
  stake = toph_data$stake, toph_data@coords, past_landuse = toph_data$past_landuse,
  EXPH = log(pointData$EXPH[pointData$d == 10]) + 1,
  CAMG = log(pointData$CAMG[pointData$d == 10] + 1),
  TOOC = log(pointData$TOOC[pointData$d == 10] + 1))
pretic_data_new <- apply(pretic_data[, c(5, 6, 7)], 2, function (x) (x - mean(x)) / sd(x))
pretic_data_new <- cbind(pretic_data[, 1:4], pretic_data_new)
sv <- c("EXPH", "CAMG", "TOOC")

# cross-variogram
sp::coordinates(pretic_data_new) <- ~ x + y
pretic_vario <- gstat::gstat(NULL, id = "EXPH", data = pretic_data_new, form = EXPH ~ past_landuse)
pretic_vario <- gstat::gstat(pretic_vario, id = "CAMG", data = pretic_data_new, form = CAMG ~ past_landuse)
pretic_vario <- gstat::gstat(pretic_vario, id = "TOOC", data = pretic_data_new, form = TOOC ~ past_landuse)
pretic_vario <- list(
  g = pretic_vario, 
  v = gstat::variogram(pretic_vario, boundaries = pedometrics::vgmLags(pretic_data_new@coords, n = 5)[-1]))
plot(pretic_vario$v, scales = list(relation = "same"), pch = 20, cex = 0.5)

# competing ranges
pretic_cross <- parallel::mclapply(1:length(range), function (i)
  gstat::gstat(
    pretic_vario$g, id = sv,
    model = gstat::vgm(psill = 0.4, model = "Exp", range = range[i], nugget = 0.2), fill.all = TRUE))
pretic_lmc <- parallel::mclapply(1:length(pretic_cross), function (i)
  gstat::fit.lmc(v = pretic_vario$v, g = pretic_cross[[i]], correct.diagonal = 1.01))

# cross-validation
pretic_cv <- parallel::mclapply(
  pretic_lmc, gstat::gstat.cv, nfold = length(unique(pointData$stake)), remove.all = TRUE, 
  all.residuals = TRUE, boundaries = attr(pretic_vario$v, "boundaries"), correct.diagonal = 1.01)
pretic_cv <- round(do.call(rbind, lapply(pretic_cv, colMeans)), 4)
apply(abs(pretic_cv), 2, which.min)
apply(abs(pretic_cv), 2, which.max)
pretic_cv <- pretic_cv[[3]]
pretic_lmc <- pretic_lmc[[1]]
plot(pretic_vario$v, pretic_lmc, scales = list(relation = "same"), pch = 20, cex = 0.5, ylim = c(0, 1.2))



d <- 10
pretic_data <- data.frame(
  stake = toph_data$stake, toph_data@coords, past_landuse = toph_data$past_landuse,
  EXPH = pointData$EXPH[pointData$d == d],
  CAMG = pointData$CAMG[pointData$d == d],
  TOOC = pointData$TOOC[pointData$d == d])

id <- sapply(1:nrow(pretic_data), function (i)
  pretic_data[i, "EXPH"] >= 30 && pretic_data[i, "CAMG"] >= 2 && pretic_data[i, "TOOC"] >= 10
  )
  
pretic_data_new$pretic <- id

tmp <- gstat::variogram(pretic ~ past_landuse, pretic_data_new, width = 50)
plot(tmp)
v.fit <- gstat::fit.variogram(tmp, gstat::vgm(0.1, "Exp", 50, 0.05))
plot(tmp, v.fit)
sp::proj4string(pretic_data_new) <- sp::proj4string(covar)
tmp <- gstat::krige(pretic ~ past_landuse, pretic_data_new, model = v.fit, newdata = covar)
tmp$var1.pred <- ifelse(tmp$var1.pred > 0.5, "TRUE", "FALSE")
tmp$var1.pred <-  as.factor(tmp$var1.pred)

dev.off()
png("res/fig/continuous_pretic.png")
sp::spplot(tmp, 1)

# plot(pretic_data[, 2:3])
points(tmp[, 2:3], col = "red", pch = 20)
text(tmp[, 2:3], labels = tmp$stake, pos = 4)
