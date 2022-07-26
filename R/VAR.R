

#' Estimate a VAR model
#'
#' Estimates and validates a VAR model
#' @param x a matrix with the data to be checked. Number of columns should be equal
#'          to the number of variables.
#' @param pmax the maximum order of the model to be estimated. Default is 22.
#' @param p a specific model order for the VAR model. Default is NULL
#' @param ic an information criterium to estimate the model order, if not previously
#'           specified. Default is AIC, the Akaike criterium
#' @param alpha significance level for testing. Default is 0.05
#' @param correction choose a p-value correction method for multiple hypotheses. Default
#'                  is bonferroni. For other methods, please check \link[stats]{p.adjust.methods}
#' @param type parameter to be passed into the \link[vars]{VAR} function, indicating the type of model.
#'             Default is none
#' @param exogen parameter to be passed into the \link[vars]{VAR} function, indicating 
#'               exogenous variables. Default is NULL
#' @param warnings boolean, show warnings regarding the test results. Default is FALSE
#'
#' @return A boolean variable indicating if the multivariate time series is stationary or not,
#'         or a caption indicating so, depending on whether the arguments verbose and warning are TRUE or
#'         FALSE
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @details 
#' This function uses the function \link[vars]{VAR} to estimate a particular VAR model. At the same
#' time, it checks the different properties of the variables and the model and thus informs about
#' the model's validity.
#'         
#' @export
#' @import vars
#'
#' @examples
#' data(DetrendedData)
#' 
#' CheckStationarity(DetrendedData)
#' CheckStationarity(DetrendedData, verbose = TRUE)

EstimateVAR <- function(x, pmax = 22, p = NULL, ic = "AIC", alpha = 0.05, 
    correction = "bonferroni", type = "none", exogen = NULL, warnings = TRUE){
                     if(!CheckStationarity(x, alpha = alpha,
                         warnings = warnings, correction = correction)){
                         #stop("Time series are not stationary")
                     } else {
                     model <- vars::VAR(x, type = type, exogen = exogen,
                       ic = ic, p = p, lag.max = pmax)
                     whiteness <- DiagnoseResiduals(model,
                       alpha = alpha, correction = correction)
                     stability <- DiagnoseStability(model)
                     if(!whiteness & warnings){
                        warning("Model residuals are not white noise. Non 
                          valid model")
                     }
                     if(!stability & warnings){
                        warning("Model is not stable. Non 
                          valid model")
                     }
                     return(model)
                     }
}



                         
               
               












                                
                
     

