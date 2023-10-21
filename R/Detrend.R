#' Detrend variables
#'
#' Detrends a vector of observations using the discrete wavelet transform
#' @param x A vector of observations
#' @param thr Manual threshold. Default is NULL.
#' @param per Boolean, use a percentage as a threshold value
#' @param cutoff A frequency used as limit for the detrending, or a reference
#'               frequency for calculating said limit. Default is 0.04 Hz
#' @param f Sample rate of the vector of observations. Default is 4 Hz
#' @param wv Wavelet to be passed into the function \link[waveslim]{modwt} of package
#'           \href{https://CRAN.R-project.org/package=waveslim}{waveslim}. Default is
#'           d16
#' @param max_f Maximum frequency to be accounted for in the detrending. default is 0.4 Hz
#' @param do.VHF Boolean. If true, Very High Frequency bands will also be processed.
#' @param use.universal Boolean. If TRUE, universal thresholding will be performed.
#' @param n.universal Number of scales in which to applied the universal thresholding, if
#'                  applied. If NULL (default), it will be computed automatically.
#' @param thr.type Perform hard or soft thresholdings
#'
#' @return The detrended vector of observations
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @details This function uses a wavelet detrending algorithm to select the aproppriate componets of a signal
#' and discard low frequency components. This function requires package \href{https://CRAN.R-project.org/package=waveslim}{waveslim}:
#' for more information regarding this package, please check the references section.
#' 
#' @references 
#' 
#' Brandon Whitcher (2020). waveslim: Basic Wavelet Routines for One-, Two-, and Three-Dimensional Signal
#' Processing. R package version 1.8.2. https://CRAN.R-project.org/package=waveslim
#' 
#' Li L, Liu C, Li K, Liu C. Comparison of Detrending Methods in Spectral Analysis of Heart Rate Variability.
#' Res J Appl Sci Eng Technol. 2011;3(9):1014-21
#'
#'         
#' @export
#' @import waveslim
#' 
#'
#' @examples
#' data(Cardiovascular)
#' int_data <- ResampleData(Cardiovascular)
#' RR <- DetrendWithCutoff(int_data$RR) # The RR series is detrended
#' SBP <- DetrendWithCutoff(int_data$SBP) # The SBP series is detrended
DetrendWithCutoff <-
  function(x,
           thr = NULL,
           per = FALSE,
           cutoff = 0.04,
           f = 4,
           wv = "d16",
           max_f = 0.4,
           do.VHF = FALSE,
           use.universal = FALSE,
           n.universal = NULL,
           thr.type = c("hard", "soft"))
  {
    thr.type <- match.arg(thr.type)
    N = NROW(x)
    Nf = f / 2
    level = 1
    repeat {
      target_f <- Nf / (2 ^ level)
      if (target_f <= cutoff) {
        break
      }
      else {
        level = level + 1
      }
    }
    level2 = 1
    if (do.VHF) {
      repeat {
        target_m <- Nf/(2^level2)
        if (target_m <= max_f) {
          break
        }
        else {
          level2 = level2 + 1
        }
      }
    }
    if (is.null(n.universal)){
      if(do.VHF){
        n.universal <- level2
      } else {
        n.universal <- level
      }
    }
    dx <-
      waveslim::modwt(x - mean(x), wf = wv, level, boundary = "reflection")
    if (use.universal)
      dx <-
      waveslim::universal.thresh.modwt(dx, n.universal, hard = (thr.type == "hard"))
    trend <- dx
    for (n in 1:level)
      trend[[n]] <- double(NROW(dx[[1]]))
    if (!is.null(thr)) {
      if (!per) {
        if (thr.type == "hard") {
          trend[[level + 1]] <-
            trend[[level + 1]] * (abs(trend[[level + 1]]) >= thr)
        } else {
          trend[[level + 1]] <-
            sign(trend[[level + 1]]) * (abs(trend[[level + 1]]) - thr) * (abs(trend[[level +
                                                                                       1]]) >= thr)
        }
      }
      if (per)
        trend[[level + 1]] <- trend[[level + 1]] * thr / 100
    }
    x = waveslim::imodwt(dx)
    trend = waveslim::imodwt(trend)
    return(x - trend)
  }
