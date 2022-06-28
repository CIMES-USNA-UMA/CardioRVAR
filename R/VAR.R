
# Developed by Alvaro Chao-Ecija

# This function automatically estimates and validates a particular
# VAR model. The VAR model is estimated thanks to the VAR function from
# vars package.

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


# References:

# vars package

                         
               
               












                                
                
     

