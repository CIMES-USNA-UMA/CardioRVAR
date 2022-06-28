# Developed by Alvaro Chao-Ecija
#
# Function for interpolating the data. It uses some default parameters, like
# other packages such as RHRV.It also supports (if available) DBP and
# respiration data.


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