
# Prepare depth-wise box-and-whisker plots ####
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
