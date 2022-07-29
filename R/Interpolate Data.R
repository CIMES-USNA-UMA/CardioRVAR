#' Interpolate data
#'
#' Interpolates data up to a particular sample rate
#' @param x a list containing data to be interpolated
#' @param f sample rate for the interpolated series. Default is 4 Hz
#'
#' @return A list containing the interpolated data.
#'         
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'
#' @export
#'
#' @examples
#' data(Cardiovascular)
#' interpolated <- InterpolateData(Cardiovascular)
#' plot(interpolated$Time, interpolated$RR, xlab = "Time (sec)", ylab = "RR intervals (ms)")

InterpolateData <- function(x, f = 4){
                   IntFunRR <- splinefun(x$Time, x$RR, 
                       method = "monoH.FC", ties = "ordered")
                   IntFunSBP <- splinefun(x$Time, x$SBP, 
                       method = "monoH.FC", ties = "ordered")
                   if(!is.null(x$DBP)){
                     IntFunDBP <- splinefun(x$Time, x$DBP, 
                        method = "monoH.FC", ties = "ordered")
                   } else {
                     IntFunDBP <- function(x) return(NULL)
                   }
                   if(!is.null(x$Resp)){
                     IntFunDBP <- splinefun(x$Time, x$Resp, 
                                            method = "monoH.FC", ties = "ordered")
                   } else {
                     IntFunResp <- function(x) return(NULL)
                   }
                   Time = seq(x$Time[1], x$Time[NROW(x$Time)], 1/f)
                   return(list(Time = Time, RR = IntFunRR(Time), 
                          SBP = IntFunSBP(Time) , DBP = IntFunDBP(Time),
                          Resp = IntFunResp(Time)))
}