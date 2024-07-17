#' Standardize a data set
#' @description
#' Use this function to "standardize" a data set. By default it will center and normalize (z-score transform) the data.
#' @param x the vector of values (e.g. timeseries) to be standardize
#' @param center T/F boolean to indicate if you want the dataset to be centered (subtract mean)
#' @param normalize T/F boolean to indicate if you want to normalize by the standard deviation
#' @author Thomas Bryce Kelly
#' @export
standardize = function(x, center = T, normalize = T) {
  xmean = 0
  if (center) {
    xmean = mean(x, na.rm = T)
  }
  deviation = 1
  if (normalize) {
    deviation = sd(x, na.rm = T)
  }
  (x - xmean) / deviation
}

#' Detrend a timeseries
detrend = function(time, signal) {
  l = order(time)
  time = time[l]
  signal = signal[l]
  
  model = lm(signal ~ time)
  
  signalDetrended = signal - (coefficients(model)[2] * time + coefficients(model)[1])
  
  ## Return
  signalDetrended
}


#' @title Calculate time integrated signal value.
#' @param time time-basis for the signal (in datetime stamps)
#' @param signal signal or vector to be integrated
#' @param tau window size
#' @param f function to apply to window entries (default: mean)
#' @param plot Flag to turn on/off plotting
#' @export
calculateIntegration = function(time,
                                signal,
                                tau = 365,
                                f = function(x) {mean(x, na.rm = T)},
                                plot = T) {
  
  l = order(time)
  time = time[l]
  signal = signal[l]
  
  signalStandardized = standardize(signal)
  signalIntegrated = rep(NA, length(signal))
  
  for (i in 2:length(signal)) {
      k = time <= time[i] & time > time[i] - tau * 86400
      signalIntegrated[i] = f(signalStandardized[k])
    }
  
  if (plot) {
    plot(time, signal, type = 'l', ylab = n, xlab = '', ylab = 'Signal')
    lines(time, signalIntegrated, col = 'red', lwd = 3)
  }
  
  ## Return
  signalIntegrated
}


#' @title Match driver and response by lagged window
#' @author Thomas Bryce Kelly
#' @description
#' Use this function to match up a given driver timeseries to a (typically lower resolution) response signal.
#' The value matched up can be lagged and a function applied over a window (e.g. 6 month lag to a 3 month running mean)
#' Note that all time units should be the same
#' @param driverTime times at which driver is given
#' @param driver the driver signal
#' @param responseTime times where you want the resulting driver signal computed
#' @param tau window size
#' @param lag number of entries to offset window
#' @param f function to apply to window entries (default: mean)
#' @export
laggedMatch = function(driverTime,
                       driver,
                       responseTime,
                       tau,
                       lag) {
  driverTime = as.numeric(driverTime)
  responseTime = as.numeric(responseTime)
  tau = as.numeric(tau)
  lag = as.numeric(lag)

  driverTime = driverTime + lag
  res = rep(NA, length(responseTime))
  for (i in 1:length(responseTime)) {
    if (min(driverTime) <= responseTime[i] - tau) {
      res[i] = mean(driver[driverTime > responseTime[i] - tau & driverTime <= responseTime[i]])
    }
  }
  res
}



#' @title Add Log Axis
#' @author Thomas Bryce Kelly
#' @description A helper function to add log axis to a plot.
#' @param side Side for the axis: 1 = bottom, 2 = left...
#' @param base the log bsae for the axis
#' @param col the color of the major axis labels
#' @param color.minor the color for the minor ticks and grid (if specified)
#' @param grid boolean, draw grid?
#' @export
addLogAxis = function(side = 1, at = NULL, labels = NULL, ticks = T, base = 10, col = 'black', color.minor = 'grey', grid = F, grid.major = F, ...) {
  at.default = c(-30:30)
  if(is.null(at)) {
    at = at.default
    if (ticks) {
      tick.pos = rep(1:(base-1), length(at)) * base^as.numeric(sapply(at, function(x) {rep(x, base-1)}))
    }
  } else {
    at = log(at, base)
    if (ticks) {
      tick.pos = rep(1:(base-1), length(at.default)) * base^as.numeric(sapply(at.default, function(x) {rep(x, base-1)}))
    }
  }
  
  ## Default labels are just the transformed values
  if (is.null(labels)) {
    labels = base^at
  }
  
  ## Draw Axis
  axis(side = side, at = log(tick.pos, base), tick = T, labels = F, col = color.minor, ...)
  axis(side = side, at = at, labels = labels, col = col, ...)
  
  if (grid & (side == 1 | side == 3)) {
    abline(v = log(tick.pos, base), col = color.minor, lty = 3)
  }
  if (grid & (side == 2 | side == 4)) {
    abline(h = log(tick.pos, base), col = color.minor, lty = 3)
  }
  if (grid.major & (side == 1 | side == 3)) {
    abline(v = at, col = color.minor, lty = 3)
  }
  if (grid.major & (side == 2 | side == 4)) {
    abline(h = at, col = color.minor, lty = 3)
  }
}


#' @title Calculate Simple Moving Average
#' @author Laura Whitmore
#' @param parameter signal or vector to be smoothed
#' @param k the integer number of entries to average over
#' @export
sma = function(parameter, k, f = mean) {
  if (k %% 2 == 1) {
    k = k - 1
  } else {
    stop('Please select an odd window size.')
  }
  k = k / 2
  sma = rep(NA, length(parameter)) ## need this vector to be the same length as parameter
  for (i in 1:length(parameter)) {
    if (i <= k | i >= length(parameter) - k + 1) {
      sma[i] = NA
    } else {
      series = c((i - k):(i + k))
      sma[i] = f(parameter[series])
    }
  }
  sma
}

