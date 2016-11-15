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
        box.dot = list(cex = 0.7)), ...,
      panel = function (...) {
        lattice::panel.grid(v = -1, h = -1)
        lattice::panel.bwplot(...)
      })
    
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
      names(lambda) <- names(sv_mean)
      a$lambda <- lambda
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
      lm_rss[[i]] <- c(anova(fit)$`Sum Sq`[1], summary(fit)$r.squared, anova(fit)$`Pr(>F)`[1])
      names(lm_rss[[i]]) <- c("SS", "R2", "P")
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

# Compute empirical variogram ----
compute_sample_variogram <- 
  function (soil_data, sv = "TOOC", depth = seq(10, 90, 20), cross = FALSE) {
    
    # Define model formula
    # We model every soil property as a linear function of the covariate
    id <- colnames(soil_data@data)[2:6]
    
    # form <- lapply(1:length(depth), function (i) paste(sv, ".", depth[i], " ~ past_landuse", sep = ""))
    form <- lapply(1:length(depth), function (i) paste(id[i], " ~ past_landuse", sep = ""))
    form <- lapply(form, as.formula)
    
    # Compute lag-distance classes
    lags <- pedometrics::vgmLags(soil_data@coords, n = 5)[-1]
    
    if (cross) {
      # Create gstat object
      for (i in 1:length(depth)) {
        if (i == 1) {
          # g <- 
          # gstat::gstat(NULL, id = paste(sv, ".", depth[i], sep = ""), data = soil_data, form = form[[i]])
          g <- gstat::gstat(NULL, id = id[i], data = soil_data, form = form[[i]])
        } else {
          # g <- gstat::gstat(g, id = paste(sv, ".", depth[i], sep = ""), data = soil_data, form = form[[i]])
          g <- gstat::gstat(g, id = id[i], data = soil_data, form = form[[i]])
        }
      }
      
      # Compute empirical direct and cross-variograms
      v <- gstat::variogram(g, boundaries = lags)
      
      return (list(g = g, v = v))
      
    } else {
      # Compute empirical direct and cross-variograms
      v <- list()
      for (i in 1:length(depth)) {
        v[[i]] <- gstat::variogram(form[[i]], boundaries = lags, dat = soil_data)
      }
      return (v)
    }
  }

# Back-transform predictions ----
back_transform <- 
  function (pred, soil_data, depth = seq(10, 90, 20), n.sim = 10000) {
    
    # Identify the columns of 'predgrid' containing the predictions and prediction error variances
    id_mean <- seq(1, 9, 2)
    id_sd <- seq(2, 10, 2)
    
    # Rescale predictions
    pred@data[, id_mean] <- 
      lapply(1:length(depth), function (i) 
        (pred@data[, id_mean[i]] * attr(soil_data, "sv_sd")[i]) + attr(soil_data, "sv_mean")[i])
    
    # Rescale prediction error variance (now standard deviation)
    pred@data[, id_sd] <- 
      lapply(1:length(depth), function (i) 
        sqrt(pred@data[, id_sd[i]]) * attr(soil_data, "sv_sd")[i])
    
    for (i in 1:length(depth)) {

      # Check which cells are not NA
      n <- which(!is.na(pred@data[, id_mean[i]]))

      # Set a progress bar
      pb <- txtProgressBar(min = 1, max = length(n), style = 3)
      k <- 0

      # Start loop over prediction locations
      for (j in n) {
        k <- k + 1

        # Simulate n values at the j-th prediction locations
        sim <- geoR::rboxcox(
          n = n.sim, lambda = attr(soil_data, "lambda")[i],
          mean = pred@data[j, id_mean[i]], sd = pred@data[j, id_sd[i]])

        # Replace the predicted value with the mean of the n simulated values
        pred@data[j, id_mean[i]] <- mean(sim) - 1

        # Replace the prediction error standard deviation with the standard deviation
        # of the n simulated values
        pred@data[j, id_sd[i]] <- sd(sim)

        setTxtProgressBar(pb, k)
      }
      close(pb)
    }

    return (pred)
  }

# Add grid lines to lattice graphics aligned with the axis labels ----
addGridLines <- latticeExtra::layer(lattice::panel.grid(h = -1, v = -1))

# Compute total stocks ----
totals <- 
  function (x) {
    
    total <- sum(x@data$stock, na.rm = TRUE)
    total_sd <- sqrt(sum(x@data$stock ^ 2, na.rm = TRUE))
    
    res <- data.frame(total = total, total_sd = total_sd)
    return (res)
  }

# Prepare figure with depth-wise predictions ----
layer_predictions <-
  function (x, var, main = "") {
    
    if (var == "pred") {
      var <- seq(1, 9, 2)
      col.regions <- soil.colors
    } else {
      var <- seq(2, 10, 2)
      col.regions <- uncertainty.colors
    }
    
    sp::spplot(
      x, var, layout = c(5, 1), col.regions = col.regions, main = main,
      strip = lattice::strip.custom(factor.levels = paste(seq(10, 90, 20), "cm")),
      panel = function (...) {
        lattice::panel.grid(h = -1, v = -1)
        lattice::panel.levelplot(...)
        d <- depth[lattice::panel.number()]
        lattice::panel.points(
          pointData@coords[pointData$d == d, ], cex = 0.5, fill = pointData$col[pointData$d == d], 
          col = pointData$col[pointData$d == d], pch = pointData$pch[pointData$d == d])
      })
  }

# Prepare figure with profile predictions ----
profile_predictions <-
  function (x, var, col.regions = soil.colors, main = "") {
    sp::spplot(
      obj = x, zcol = var, col.regions = col.regions, 
      main = main,
      panel = function (...) {
        lattice::panel.grid(h = -1, v = -1)
        lattice::panel.levelplot(...)
        lattice::panel.points(
          pointData@coords, fill = pointData$col, col = pointData$col, pch = pointData$pch)
      })
  }

# Compute the slope of the exponential covariance function
slope <-
  function (lmc) {
    sapply(lmc$model, function (x) sum(x$psill) / (0.333 * x$range[2]))
  }
