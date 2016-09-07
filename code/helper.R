
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
      scales = list(y = "same", x = "free", alternating = FALSE), 
      xlab = "Value",
      # xlab = expression(paste('Content (g ',kg^-1,')', sep = '')),
      par.settings = list(
        fontsize = list(text = 14, points = 8), box.rectangle = list(col = "black"),
        box.umbrella = list(col = "black"), plot.symbol = list(col = "black", cex = 0.7),
        box.dot = list(cex = 0.7)), ...)
    
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
    soil_var <- pointData@data[, c("stake", "d", sv)]
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
    soil_var <- cbind(soil_var, pointData@coords[(1:nrow(soil_var)) * 5, ])
    sp::coordinates(soil_var) <- ~ x + y
    sp::proj4string(soil_var) <- sp::proj4string(covar)
    
    # Sample covariate at calibration points
    soil_var$past_landuse <- sp::over(soil_var, covar)$past_landuse
    
    # Save data for back-tranformation
    if (save4back) {
      a <- attributes(soil_var)
      a$sv_mean <- sv_mean
      a$sv_sd <- sv_sd
      attributes(soil_var) <- a
    }
    
    return(soil_var)
  }

# Get lm fit output ----
get_lm_output <-
  function (soil_var = soil_var, sv = "TOOC", d = seq(10, 90, 20)) {
    lm_rss <- list()
    lm_coef <- list()
    lm_error <- list()
    
    for (i in 1:length(d)) {
      form <- as.formula(paste(sv, ".", d[i], " ~ past_landuse", sep = ""))
      
      # Fit linear models, run an analsis of variance, and extract the sum of squares of the fit,
      # the coefficient of determination, fitted coefficients and their standard error
      fit <- lm(form, data = soil_var@data)
      lm_rss[[i]] <- c(anova(fit)$`Sum Sq`[1], summary(fit)$r.squared)
      names(lm_rss[[i]]) <- c("SS", "R2")
      lm_coef[[i]] <- summary(fit)$coefficients[, 1]
      names(lm_coef[[i]]) <- c("Intercept", "exp(past)")
      lm_error[[i]] <- summary(fit)$coefficients[, 2]
      names(lm_error[[i]]) <- c("Intercept", "exp(past)")
    }
    # Put all data in a table so that we can check it.
    lm_rss <- round(do.call(rbind, lm_rss), digits = 2)
    lm_coef <- round(do.call(rbind, lm_coef), digits = 2)
    lm_error <- round(do.call(rbind, lm_error), digits = 2)
    
    tmp <- list(lm_rss, lm_coef, lm_error)
    tmp2 <- matrix(
      NA_character_, nrow = nrow(tmp[[2]]), ncol = ncol(tmp[[2]]), 
      dimnames = list(paste(seq(0, 80, 20), "-", seq(20, 100, 20), sep = ""), colnames(tmp[[2]])))
    
    for (i in 1:ncol(tmp[[2]])) {
      for (j in 1:nrow(tmp[[2]])) {
        tmp2[j, i] <- paste(tmp[[2]][j, i], " (", tmp[[3]][j, i], ")", sep = "")
      }
    }
    
    tmp <- cbind(tmp2, tmp[[1]])
    
    return (tmp)
  }
