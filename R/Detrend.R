#' Detrend variables
#'
#' Detrends a vector of observations using the discrete wavelet transform
#' @param x a vector of observations
#' @param cutoff a frequency used as limit for the detrending. Default is 0.04 Hz
#' @param f sample rate of the vector of observations. Default is 4 Hz
#' @param wv wavelet to be passed into the \link[waveslim]{modwt} function. Default is
#'           d16
#' @param max_f maximum frequency to be accounted for in the detrending. default is 0.4 Hz
#'
#' @return The detrended vector of observations
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#'         
#' @export
#' @import waveslim
#'
#' @examples
#' data(Cardiovascular)
#' data <- InterpolateData(data)
#' RR <- DetrendByCutoff(data$RR) # The RR series is detrended
#' SBP <- DetrendByCutoff(data$SBP) # The SBP series is detrended
DetrendByCutoff <- function(x, cutoff = 0.04, f = 4, wv = "d16", 
                            max_f = 0.4){
  N = NROW(x)
  Nf = f/2
  level = 1
  repeat{
    target_f <- Nf / (2^level)
    if(target_f <= cutoff){
      break
    } else {
      level = level + 1
    }
  }
  level2 = 1
  repeat{
    target_m <- Nf / (2^level2)
    if(target_m <= max_f){
      break
    } else {
      level2 = level2 + 1
    }
  }
  dx = waveslim::modwt(x - mean(x), wf = wv, level)
  dx <- waveslim::universal.thresh.modwt(dx, level, hard = TRUE)
  for(n in level2:(level-1)) dx[[n]] <- double(N)
  trend <- waveslim::imodwt(dx)
  return(x - mean(x) - trend)
}


