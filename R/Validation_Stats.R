
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
#'                  is bonferroni. For other methods, please check \link[stats]{p.adjust}
#'                  and \link[stats]{p.adjust.methods} from package
#'                  \href{https://CRAN.R-project.org/package=stats}{stats}
#' @param verbose.method method for the verbose option, either "c" (cat) or
#'                  "p" (print). Default is c.
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, depending on whether the arguments verbose and warning are TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @details 
#' Several functions are available in this package to statistically validate
#' properties of the evaluated time series and VAR models. These functions are based
#' on the validation algorithms proposed in Seth (2010) and Seth (2011) for MVAR validation.
#' The functions allow corrections for multiple comparisons offered by \link[stats]{p.adjust} 
#' and in \link[stats]{p.adjust.methods}. The bonferroni correction has been selected as
#' default, as suggested in Seth (2010). 
#'  
#' In this function, stationarity of the time series is evaluated according to the criteria
#' established in Seth (2010) and more specifically in Seth (2011).
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
  verbose = FALSE, correction = "bonferroni", verbose.method = c("c", "p")){
   verbose.method <- match.arg(verbose.method)
   catfun <- function(x) {cat(paste("\n", x, "\n"))}
   verbosefun <- ifelse(verbose.method == "c", catfun , print)
   N <- dim(x)[1]
   M <- dim(x)[2]
   if(is.null(lags)){
      lags <- sqrt(N) # As suggested in Seth (2011)
   }
   adf_p_vals <- double(M)
   kpss_p_vals <- double(M)
   w <- getOption("warn")   
   options(warn = -1)
   for(m in 1:M){
       adf_p_vals[m] <- tseries::adf.test(x[,m], k = lags)$p.value
       kpss_p_vals[m] <- tseries::kpss.test(x[,m], lshort = F)$p.value
   }
   options(warn = w)
   adf_p_vals <- p.adjust(adf_p_vals, correction)
   kpss_p_vals <- p.adjust(kpss_p_vals, correction)
   non_sig_adf <- (1:M)[adf_p_vals > alpha]
   sig_kpss <- (1:M)[kpss_p_vals <= alpha]
   if(sum(sig_kpss) != 0){
      if(warnings){
        warning(paste("\n", "Time series", sig_kpss, "is not stationary", "\n"))
        return(FALSE)
      } else if(!warnings & verbose & verbose.method == "p"){
        verbosefun(paste("Time series", sig_kpss, "is not stationary"))
      } else if(!warnings & verbose.method != "p"){
        if(verbose) verbosefun(paste("Time series", sig_kpss, "is not stationary"))
        return(FALSE)
      }
   } else if(sum(non_sig_adf) != 0){
      if(warnings){
      warning(paste("\n", "ADF testing suggest non stationarity for time series", non_sig_adf,
       "\n"))
      }
      if(sum(sig_kpss) != 0){
        if(warnings){
          warning(paste("\n", "Time series", sig_kpss, "is not stationary", "\n"))
          return(FALSE)
        } else if(!warnings & verbose){
          return(paste("\n", "Time series", sig_kpss, "is not stationary", "\n"))
        } else if(!warnings & !verbose){
          return(FALSE)
        }
      } else {
      if(warnings){
      message(paste("Time series are stationary according to KPSS test"))
      }
        if(verbose & (verbose.method == "p")){
          verbosefun("Time series are stationary")
        } else {
          if(verbose) verbosefun("Time series are stationary")
          return(TRUE)
        }
      }
   } else {
      if(verbose & (verbose.method == "p")){
        verbosefun("Time series are stationary")
      } else {
        if(verbose) verbosefun("Time series are stationary")
        return(TRUE)
      }
   }
}



#' Check white-noise nature of residuals
#'
#' Check that residuals from a VAR model are a white-noise process
#' @param model a VAR model to be checked
#' @param alpha significance level for testing. Default is 0.05
#' @param correction choose a p-value correction method for multiple hypotheses. Default
#'                  is bonferroni. For other methods, please check \link[stats]{p.adjust}
#' @param verbose boolean, show results from the test in captions. Default is FALSE
#' @param verbose.method method for the verbose option, either "c" (cat) or
#'                  "p" (print). Default is c.
#'
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, whether the argument verbose is TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @details 
#' Several functions are available in this package to statistically validate
#' properties of the evaluated time series and VAR models. These functions are based
#' on the validation algorithms proposed in Seth (2010) and Seth (2011) for MVAR validation.
#' The functions allow corrections for multiple comparisons offered by \link[stats]{p.adjust} 
#' and in \link[stats]{p.adjust.methods}. The bonferroni correction has been selected as
#' default, as suggested in Seth (2010). 
#'  
#' In this function, residuals obtained from the estimated models are evaluated according 
#' to the criteria established in Seth (2010) and Seth (2011).
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
    verbose = FALSE, verbose.method = c("c", "p")){
    verbose.method <- match.arg(verbose.method)
    catfun <- function(x) {cat(paste(x, "\n"))}
    verbosefun <- ifelse(verbose.method == "c", catfun , print)
    lms <- model$varresult
    p_vals <- double(length(lms))
    for(n in 1:length(lms)){
        mod <- lms[[n]]
        p_vals[n] <- lmtest::dwtest(mod)$p
    }
    p_vals <- p.adjust(p_vals, correction)
    if(verbose & (max(p_vals) <= alpha)){
       verbosefun("Residuals are not white processes")
       if(verbose.method != "p") return(FALSE)
    } else if(verbose & (min(p_vals) > alpha)){
       verbosefun("Model residuals are white noise processes")
       if(verbose.method != "p") return(TRUE)
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
#' @param verbose.method method for the verbose option, either "c" (cat) or
#'                  "p" (print). Default is c.
#'
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, whether the argument verbose is TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' #' @details 
#' Several functions are available in this package to validate
#' properties of the evaluated time series and VAR models. In this function, 
#' stability of the VAR model is evaluated according to the criteria
#' established in Ding (2000).
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
                              col2 = "red", verbose.method = c("c", "p")){
    verbose.method <- match.arg(verbose.method)
    catfun <- function(x) {cat(paste(x, "\n"))}
    verbosefun <- ifelse(verbose.method == "c", catfun , print)
    roots <- vars::roots(model)
    if(do.plot) PlotRoots(model)
    if(verbose & !(TRUE %in% (roots >= 1))){
       verbosefun("The model is stable")
       if(verbose.method != "p") return(TRUE)
    } else if(verbose & !(TRUE %in% (roots < 1))){
       verbosefun("The model is not stable")
      if(verbose.method != "p") return(FALSE)
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
PlotRoots <- function(model, col = "blue", col2 = "red"){
  roots <- vars::roots(model, FALSE)
  roots_outside <- roots[abs(roots) > 1]
  plot(Re(roots), Im(roots), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), col = col)
  if(length(roots_outside) != 0) points(roots_outside, col = col2)
  lines(seq(-1, 1, len = 1000), sqrt(1 - seq(-1, 1, len = 1000)^2))
  lines(seq(-1, 1, len = 1000), -sqrt(1 - seq(-1, 1, len = 1000)^2))
}















