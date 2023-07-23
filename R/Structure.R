



GetModelStructure <- function(){
           output <- list()
           output$Data <- list()
           output$TimeDomain <- list()
           output$FreqDomain <- list()
           return(output)
}


BuildTimeDomain <- function(model, pmax = 22, p = NULL, ic = "AIC", alpha = 0.05, 
    bonferroni = FALSE, type = "none", exogen = NULL){
           x = model$Data
           TimeModel <- EstimateVAR(x, pmax, p, ic, alpha, bonferroni, type,
             exogen)
           order <- TimeModel$p
           coefs <- GetCoefs(TimeModel)
           resids <- resid(TimeModel)
           sigma <- summary(TimeModel)$cov
           output <- list()
           output$Coefs <- coefs
           output$Resids <- resids
           output$Sigma <- sigma
           output$Lags <- order
           model$TimeDomain$Original <- output
           return(model)
}

UpdateTimeDomain <- function(model){
           TimeModel <- model$TimeDomain$Original
           coefs <- TimeModel$Coefs
           sigma <- TimeModel$Sigma
           a0 <- GetA0Fun(sigma)$a0
           sigma <- GetA0Fun(sigma)$sigma
           coefs <- UpdateWithA0(a0, coefs)
           output <- list()
           output$Coefs <- coefs
           output$Sigma <- sigma
           output$A0 <- a0 #CMD check
           output$Lags <- TimeModel$Lags
           model$TimeDomain$Updated <- output
           return(model)
}


BuildFreqDomain <- function(model, len = 1000){
           coefs <- model$TimeDomain$Original$Coefs
           sigma <- model$TimeDomain$Original$Sigma
           output <- ParamFreqModel(coefs = coefs, sigma = sigma, len = len) 
           model$FreqDomain <- output
           return(model)
}










