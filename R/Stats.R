# Developed by Alvaro Chao-Ecija

# Several functions have been developed, according to the proposed
# strategy by Seth for MVAR validation. The functions allow corrections
# for multiple comparisons offered by p.adjust.methods. The bonferroni
# correction has been selected as default, as suggested by Seth. 

# Therefore, these functions evaluate the stationarity of the 
# signals. They also test if the residuals of the MVAR model are
# white-noise, and evaluate the stability of the model.
#
# References: ---------------------------------------------------
#
#
#
# Seth AK. A MATLAB toolbox for Granger causal connectivity analysis. J 
# Neurosci Methods. 2010;186(2):262-73.
#
# Seth AK. Granger Causal Conectivity Analysis: A MATLAB Toolbox. 2011.
#
# Ding M, Bressler S, Yang W, Liang H. Short-window spectral analysis of cortical 
# event-related potentials by adaptative multivariate autoregressive modelling: data
# preprocessing, model validation, and variability assessment. Biol Cybern. 2000;83:35-45


#' Check stationarity
#'
#' Check stationarity of the data
#' @param x a matrix with the data to be checked. Number of columns should be equal
#'          to the number of variables.
#' @param alpha significance level for testing. Default is 0.05
#' @param lags specify a certain level of lags for the test. Default is NULL
#' @param warnings boolean, show warnings regarding the test results. Default is TRUE
#' @param verbose boolean, show results from the test in captions. Default is FALSE
#' @param correction choose a p-value correction method for multiple hypotheses. Default
#'                  is bonferroni. For other methods, please check \link[stats]{p.adjust} from
#'                  package \href{https://CRAN.R-project.org/package=stats}{stats}
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, depending on whether the arguments verbose and warning are TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#'         
#' @references 
#' Seth AK. A MATLAB toolbox for Granger causal connectivity analysis. J 
#' Neurosci Methods. 2010;186(2):262-73.
#'
#' Seth AK. Granger Causal Conectivity Analysis: A MATLAB Toolbox. 2011.
#'
#' Ding M, Bressler S, Yang W, Liang H. Short-window spectral analysis of cortical 
#' event-related potentials by adaptative multivariate autoregressive modelling: data
#' preprocessing, model validation, and variability assessment. Biol Cybern. 2000;83:35-45
#' @export
#' @import tseries
#'
#' @examples
#' data(DetrendedData)
#' 
#' CheckStationarity(DetrendedData)
#' CheckStationarity(DetrendedData, verbose = TRUE)
CheckStationarity <- function(x, alpha = 0.05, lags = NULL, warnings = TRUE,
  verbose = FALSE, correction = "bonferroni"){
   N <- dim(x)[1]
   M <- dim(x)[2]
   if(is.null(lags)){
      lags <- sqrt(N)
   }
   adf_p_vals <- double(M)
   kpss_p_vals <- double(M)
   w <- getOption("warn")   
   options(warn = -1)
   for(m in 1:M){
       adf_p_vals[m] <- tseries::adf.test(x[,m], k = lags)$p.value
       kpss_p_vals[m] <- tseries::kpss.test(x[,m])$p.value
   }
   options(warn = w)
   adf_p_vals <- p.adjust(adf_p_vals, correction)
   kpss_p_vals <- p.adjust(kpss_p_vals, correction)
   non_sig_adf <- (1:M)[adf_p_vals > alpha]
   sig_kpss <- (1:M)[kpss_p_vals <= alpha]
   if(sum(sig_kpss) != 0){
      if(warnings){
      warning(paste("\n", "Time series", sig_kpss, "is not stationary", "\n"))
      }
      return(FALSE)
   } else if(sum(non_sig_adf) != 0){
      if(warnings){
      warning(paste("\n", "ADF testing suggest non stationarity for time series", non_sig_adf,
       "\n"))
      }
      if(sum(sig_kpss) != 0){
         if(warnings){
         warning(paste("\n", "Time series", sig_kpss, "is not stationary", "\n"))
         }
         return(FALSE)
      } else {
      if(warnings){
      warning(paste("Time series are stationary according to KPSS test"))
      }
      return(TRUE)
      }
   } else {
      if(verbose){
         return("Time series are stationary")
      } else {
         return(TRUE)
      }
   }
}



#' Check white-noise nature of residuals
#'
#' Check that residuals from a VAR model are a white-noise process
#' @param model a VAR model to be checked
#' @param alpha significance level for testing. Default is 0.05
#' @param lags specify a certain level of lags for the test. Default is NULL
#' @param correction choose a p-value correction method for multiple hypotheses. Default
#'                  is bonferroni. For other methods, please check \link[stats]{p.adjust}
#' @param verbose boolean, show results from the test in captions. Default is FALSE
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, whether the argument verbose is TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' 
#'         
#' @references 
#' Seth AK. A MATLAB toolbox for Granger causal connectivity analysis. J 
#' Neurosci Methods. 2010;186(2):262-73.
#'
#' Seth AK. Granger Causal Conectivity Analysis: A MATLAB Toolbox. 2011.
#'
#' Ding M, Bressler S, Yang W, Liang H. Short-window spectral analysis of cortical 
#' event-related potentials by adaptative multivariate autoregressive modelling: data
#' preprocessing, model validation, and variability assessment. Biol Cybern. 2000;83:35-45
#' @export
#' @import lmtest
#'
#' @examples
#' data(DetrendedData)
#' 
#' model <- EstimateVAR(DetrendedData)
#' DiagnoseResiduals(model)
DiagnoseResiduals <- function(model, alpha = 0.05, correction = "bonferroni",
   verbose = FALSE){
    lms <- model$varresult
    p_vals <- double(length(lms))
    for(n in 1:length(lms)){
        mod <- lms[[n]]
        p_vals[n] <- lmtest::dwtest(mod)$p
    }
    p_vals <- p.adjust(p_vals, correction)
    if(verbose & (min(p_vals) <= alpha)){
       print("Residuals are not white processes")
    } else if(verbose & (min(p_vals) > alpha)){
       print("Model residuals are white processes")
    } else {
       return(min(p_vals) > alpha)
    }
}



#' Check stability of a VAR model
#'
#' Check that the estimated VAR model is stable
#' @param model a VAR model to be checked
#' @param do.plot boolean, plot a unitary circle showing the inverse roots of the
#'                model. Default is FALSE
#' @param verbose boolean, show results from the test in captions. Default is FALSE
#' @param col color for the roots inside the unit circle. Default is blue
#' @param col2 color for the roots outside the unit circle. Default is red
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, whether the argument verbose is TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' 
#'         
#' @references 
#' Seth AK. A MATLAB toolbox for Granger causal connectivity analysis. J 
#' Neurosci Methods. 2010;186(2):262-73.
#'
#' Seth AK. Granger Causal Conectivity Analysis: A MATLAB Toolbox. 2011.
#'
#' Ding M, Bressler S, Yang W, Liang H. Short-window spectral analysis of cortical 
#' event-related potentials by adaptative multivariate autoregressive modelling: data
#' preprocessing, model validation, and variability assessment. Biol Cybern. 2000;83:35-45
#' @export
#' @import vars
#'
#' @examples
#' data(DetrendedData)
#' 
#' model <- EstimateVAR(DetrendedData)
#' DiagnoseStability(model, do.plot = TRUE)
DiagnoseStability <- function(model, do.plot = FALSE, verbose = FALSE, col = "blue",
                              col2 = "red"){
    roots <- vars::roots(model)
    if(do.plot) PlotRoots(model)
    if(verbose & !(TRUE %in% (roots >= 1))){
       print("The model is stable")
    } else if(verbose & !(TRUE %in% (roots < 1))){
       print("The model is not stable")
    } else {
       return(!(TRUE %in% (roots >= 1)))
    }
}


#' Plot VAR model inverse roots
#'
#' Plot a unit circle with inverse roots computed from a VAR model
#' @param model a VAR model to be checked
#' @param col color for the roots inside the unit circle. Default is blue
#' @param col2 color for the roots outside the unit circle. Default is red
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, whether the argument verbose is TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @keywords internal
#'
#' @examples
#' data(DetrendedData)
#' 
#' model <- EstimateVAR(DetrendedData)
#' PlotRoots(model)
PlotRoots <- function(model, col = "blue", col2 = "red"){
  roots <- vars::roots(model, FALSE)
  roots_outside <- roots[abs(roots) > 1]
  plot(Re(roots), Im(roots), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), col = col)
  if(length(roots_outside) != 0) points(roots_outside, col = col2)
  lines(seq(-1, 1, len = 1000), sqrt(1 - seq(-1, 1, len = 1000)^2))
  lines(seq(-1, 1, len = 1000), -sqrt(1 - seq(-1, 1, len = 1000)^2))
}















