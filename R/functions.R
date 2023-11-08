#' Standardize a data set
#' @description
#' Use this function to "standardize" a data set. By default it will center and normalize (z-score transform) the data.
#' @param x the vector of values (e.g. timeseries) to be standardize
#' @param center T/F boolean to indicate if you want the dataset to be centered (subtract mean)
#' @param normalize T/F boolean to indicate if you want to normalize by the standard deviation
#' @author Thomas Bryce Kelly
#' @export
#'
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


#' @title Calculate time integrated signal value.
#' @author Laura Whitmore
#' @param parameter signal or vector to be integrated
#' @param tau window size
#' @param lag number of entries to offset window
#' @param f function to apply to window entries (default: mean)
#' @export
lagIntegrate = function(parameter, tau, lag, f = mean) {
  int = rep(NA, length(parameter)) ## need this vector to be the same length as parameter
  for (i in (tau+lag):length(parameter)) {
    int[i] = f(parameter[(i-tau+1):i])
  }
  int
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


