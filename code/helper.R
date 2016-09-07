
# Start OS dependent GRASS GIS ----
if (.Platform$OS.type == "unix") {
  gisBase <- "/usr/lib/grass70/"
} else {
  gisBase <- "C:/Program Files (x86)/GRASS GIS 7.0.4-1"
}

# OS dependent function to run GRASS features ----
grassGis <- function (cmd) {
  if (.Platform$OS.type == "unix") {
    system(cmd)
  } else {
    shell(cmd)
  }
}

# Prepare depth-wise box-and-whisker plots -----
depth_bwplot <-
  function (pts = pointData, vars = c("TOOC", "TOCA", "TOPH"), ...) {
    
    # Pre-process SpatialPointsDataFrame
    tmp <- as.data.frame(pointData)[, c("d", vars)]
    
    # Stack soil variables
    tmp <- cbind(d = rep(tmp$d, length(vars)), stack(x = tmp, select = vars))
    
    # Create bwplot
    a <- lattice::bwplot(
      rev(d) ~ values | ind, data = tmp, ylab = "Sampling depth (cm)", 
      scales = list(y = "same", x = "free"), xlab = expression(paste('Content (g ',kg^-1,')', sep = '')),
      par.settings = list(
        fontsize = list(text = 14, points = 8), box.rectangle = list(col = "black"),
        box.umbrella = list(col = "black"), plot.symbol = list(col = "black", cex = 0.7),
        box.dot = list(cex = 0.7)),
      strip = lattice::strip.custom(bg = "lightgray"), ...)
    
    # Reverse the y axis
    a$y.limits <- rev(c("0-20", "20-40", "40-60", "60-80", "80-100"))
    
    # Output
    return (a)
  }

# Map colors ----
soil.colors <- 
  colorRampPalette(c("ivory", "orange3", "saddlebrown", "black"))
uncertainty.colors <- 
  colorRampPalette(c("olivedrab", "khaki", "maroon1"))

# Prepare soil data for modelling ----
prepare_soil_data <- 
  function (pointData, sv = "TOOC", covar, save4back = FALSE) {
    
    # Soil variables have a skeewed distribution and have to be transformed.
    # We add one unit to its values to avoid zeros. This is necessary to perform the Box-Cox transformations.
    soil_var <- pointData[, c("stake", "d", sv)]
    soil_var[, sv] <- soil_var[, sv] + 1 
    
    # Create five new variables based on sampling depth.
    soil_var <- reshape(soil_var, idvar = "stake", timevar = "d", direct = "wide")
    
    # Lambda values of the Box-Cox transformation are estimated for each depth-variable separately. Negative 
    # lambda values are set to zero to avoid problems with the back-transform (see Samuel-Rosa et al. (2015)
    # for more information).
    tmp <- soil_var[, -1]
    lambda <- lapply(tmp, car::powerTransform)
    lambda <- sapply(lambda, function (x) x$lambda)
    lambda[which(lambda < 0)] <- 0
    tmp <- lapply(1:length(lambda), function (i) car::bcPower(tmp[, i], lambda[i]))
    soil_var[, -1] <- tmp
    rm(tmp)
    
    # Standardize the Box-Cox transformed soil variable to zero mean and unit standard deviation. This is
    # needed to help fitting the linear mixed model. It puts all depth-variables in the same scale.
    sv_mean <- sapply(soil_var[, -1], mean)
    sv_sd <- sapply(soil_var[, -1], sd)
    soil_var[, -1] <- lapply(1:5, function (i) (soil_var[, -1][, i] - sv_mean[i]) / sv_sd[i])
    
    # Get coordinates and set the coordinate reference system
    soil_var <- cbind(soil_var, pointData[(1:nrow(soil_var)) * 5, c("x", "y")])
    sp::coordinates(soil_var) <- ~ x + y
    sp::proj4string(soil_var) <- sp::CRS("+init=epsg:32720")
    
    # Sample covariate at calibration points
    soil_var$past_landuse <- sp::over(soil_var, covar)$past_landuse
    
    # Save data for back-tranformation
    if (save4back) {
      a <- attributes(soil_var)
      if (!sv %in% c("PH", "CLAY")) {
        names(lambda) <- names(sv_mean)
        a$lambda <- lambda
      }
      a$sv_mean <- sv_mean
      a$sv_sd <- sv_sd
      attributes(soil_var) <- a
    }
    
    return(soil_var)
  }
