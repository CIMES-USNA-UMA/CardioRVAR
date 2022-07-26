#' Parametric frequency model estimation
#'
#' Estimates frequency model parameters from a time domain VAR model.
#' @param model a time domain VAR model
#' @param len length of the vector of frequencies. Default is 1000
#' @param dt inverse of sample rate. Default is 0.25
#' @param A0 boolean; add 0 effects. Default TRUE
#' @param sigma include a specific covariance matrix; otherwise it is obtained
#'              from the VAR model. Default is NULL
#' @param coefs include VAR coefficients; otherwise these are obtained from the
#'              VAR model. Default is NULL
#'
#' @return A list with the estimated components of the frequency domain model.
#'         These components include closed-loop transfer functions, open-loop
#'         transfer functions, open-vs-closed-loop overestimations, noise transfer
#'         functions, and spectral estimations of the variables.
#'         
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' @export
#'
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
 
ParamFreqModel <- function(model, len = 1000, dt = 0.25, A0 = TRUE, sigma = NULL,
   coefs = NULL){
                   if(is.null(sigma) & is.null(coefs)){
                      coefs <- GetCoefs(model)
                      sigma <- summary(model)$cov
                   } else if(is.null(sigma) & !is.null(coefs)){
                             stop("Indicate coeficients")
                   } else if(is.null(sigma) & !is.null(coefs)){
                             stop("Indicate covariance matrix")
                   }
                   # For compatibility with other packages such as grangers.
                   freqs <- spec.pgram(model$y[,1], plot = FALSE, 
                                       pad = len/(model$totobs/2) - 1)$freq
                   B <- GetMatrixBfromCoefs(coefs, freqs) 
                   A <- GetMatrixAfromB(B)
                   H <- GetMatrixHfromA(A)
                   TFuns <- GetTransFuns(A)
                   sigma <- sigma * 2
                   S <- GetSpectra(H, sigma)
                   OpenTFuns <- GetOpenTFuns(S)
                   OpenVSClosed <- GetOpenVSClosedDif(OpenTFuns, TFuns)
                   a0 <- diag(1,2)
                   output <- list(Freqs = freqs / dt, Vars_Transfer_funs = B,
                      Transfer_Functions = TFuns, Open_Transfer_Functions = 
                        OpenTFuns, Open_Transfer_Functions_Overestimation =
                           OpenVSClosed, Noise_Transfer_fun = H,
                             Spectra = S, Noise_Spectra = sigma, a0 = a0)
                   if(A0){
                      output <- IncludeA0Effects(output)
                   }
                   return(output)
}


#' Include A0 effects 
#'
#' Estimates and includes A0 effects for a frequency domain model
#' @param system the model of the system whose instantaneous effects we want to include
#'
#' @return A list with the estimated components of the frequency domain model.
#'         These components include closed-loop transfer functions, open-loop
#'         transfer functions, open-vs-closed-loop overestimations, noise transfer
#'         functions, and spectral estimations of the variables.
#'         
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @export
#'
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model, A0 = FALSE)
#' 
#' #To include the A0 function:
#' freq_model <- IncludeA0Effects(freq_model)
IncludeA0Effects <- function(system){
                sigma <- system$Noise_Spectra 
                Noise_Transfer_fun <- system$Noise_Transfer_fun
                Vars_Transfer_fun <- system$Vars_Transfer_funs
                Open_Transfer_Functions <- system$Open_Transfer_Functions
                a0_Effects <- GetA0Fun(sigma)
                new_sigma <- a0_Effects$sigma
                a0 <- a0_Effects$a0
                I <- diag(c(1,1), 2, 2)
                B <- array(0, dim = dim(Vars_Transfer_fun))
                for(n in 1:(dim(Noise_Transfer_fun)[3])){
                    Noise_Transfer_fun[,,n] <- Noise_Transfer_fun[,,n] %*% 
                       solve(a0)
                }
                A <- GetMatrixAfromH(Noise_Transfer_fun)
                Vars_Transfer_fun <- GetMatrixBfromA(A, a0)
                TFuns <- GetTransFuns(A)
                OpenVSClosed <- GetOpenVSClosedDif(Open_Transfer_Functions, 
                    TFuns)
                system$Noise_Spectra <- new_sigma 
                system$Noise_Transfer_fun <- Noise_Transfer_fun
                system$Vars_Transfer_funs <- Vars_Transfer_fun
                system$Transfer_Functions <- TFuns
                system$Open_Transfer_Functions_Overestimation <- 
                     OpenVSClosed
                system$a0 <- a0
                return(system)
}





#' Estimate A0 function
#'
#' Estimates the A0 function from a covariance matrix by applying an LDLT decomposition
#' @param sigma a covariance matrix
#'
#' @return A list which includes the A0 function as well as the new covariance matrix
#' @author Alvaro Chao-Ecija
#' @author Marc Stefan Dawid-Milner
#'         
#' @keywords internal      
#'         
#' @references
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetA0Fun <- function(sigma){
            b <- t(chol(sigma))
            D <- diag(b)
            b <- b %*% diag(1/D, nrow = length(D))
            a0 <- solve(b)
            new_sigma <- diag(D^2, nrow = length(D))
            return(list(a0 = a0, sigma = new_sigma))
}

#' Update coefficients with A0 function
#'
#' Updates VAR model coefficients using an A0 function
#' @param A0 the A0 function
#' @param coefs coefficients to update
#'
#' @return The updated coefficients
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'   
#' @keywords internal
#'       
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
UpdateWithA0 <- function(A0, coefs){
            ncoefs <- coefs
            for(n in 1:length(coefs)){
                ncoefs[[n]] <- A0 %*% coefs[[n]]
            }  
            return(ncoefs)
}


#' Get coefs from VAR model
#'
#' Extracts VAR model coefficients for further analysis
#' @param var a VAR model
#'
#' @return The extracted coefficients from this VAR model
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @keywords internal
#'        
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetCoefs <- function(var){
            max.lag = var$p
            K <- var$K
            lags <- 1:max.lag
            coefs <- coef(var)
            sel_coefs <- NULL
            coef_list <- list()
            length(coef_list) <- max.lag
            for(k in 1:K){
                sel_coefs <- cbind(sel_coefs, coefs[[k]][,1])
            }
            for(lag in lags){
                coef_list[[lag]] <- t(sel_coefs[(1:K) + (K*(lag-1)),])
            }
            return(coef_list)
}



#' Calculate transfer function B from VAR model
#'
#' Calculates transfer function B from a specific VAR model
#' @param var a var model
#' @param freqs a vector of frequencies 
#' @param dt inverse sample rate. Default is 0.25
#'
#' @return The auxiliary matrix B
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'       
#' @keywords internal     
#'       
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetMatrixBfromVAR <- function(var, freqs, dt = 0.25){
              K <- var$K
              lags <- 1:var$p
              if(is.null(coefs)) coefs <- GetCoefs(var)
              B <- array(0, dim = c(K, K, NROW(freqs)))
              for(n in 1:K){
                  for(m in 1:K){
                      for(lag in lags){
                          B[n, m, ] <- B[n, m, ] +
                              coefs[[lag]][n, m] * exp(-2*pi*freqs*1i*lag)
                      }
                  }
              }
              return(B)
}


#' Calculate transfer function B from VAR coefficients
#'
#' Calculates transfer function B from a specific VAR model coefficients
#' @param coefs coefficients from a var model
#' @param freqs a vector of frequencies 
#' @param dt inverse sample rate. Default is 0.25
#'
#' @return The auxiliary matrix B
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' 
#' @keywords internal
#' 
#' 
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetMatrixBfromCoefs <- function(coefs, freqs, dt = 0.25){
              K <- ncol(coefs[[1]])
              lags <- 1:length(coefs)
              B <- array(0, dim = c(K, K, NROW(freqs)))
              for(n in 1:K){
                  for(m in 1:K){
                      for(lag in lags){
                          B[n, m, ] <- B[n, m, ] +
                               coefs[[lag]][n, m] * exp(-2*pi*freqs*1i*lag) 
                      }
                  }
              }
              return(B)
}


#' Calculate transfer function B from transfer function A
#'
#' Calculates transfer function B from a previously calculated transfer function A
#' @param B transfer function B
#'
#' @return Transfer function A
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'        
#' @keywords internal    
#'        
#'        
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetMatrixAfromB <- function(B){
              A <- array(0, dim(B))
              for(n in 1:dim(B)[3]){
                  A[ , , n] <- diag(c(1,1), dim(B)[1]) - B[ , , n]
              } 
              return(A)
}

#' Calculate noise transfer function from transfer function A
#'
#' Calculates noise transfer function H from a previously computed transfer function A
#' @param A transfer function A
#'
#' @return The noise transfer function H
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'     
#' @keywords internal 
#'     
#'        
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetMatrixHfromA <- function(A){
              H <- array(0, dim(A))
              for(n in 1:dim(A)[3]){
                  H[ , , n] <- solve(A[ , , n])
              } 
              return(H)
}

#' Calculate closed-loop transfer functions
#'
#' Calculates closed-loop transfer functions from a previously computed transfer function A
#' @param A transfer function A
#'
#' @return Closed-loop transfer functions
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'       
#' @keywords internal  
#'       
#'       
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetTransFuns <- function(A){
              TFuns <- array(0, dim(A))
              dims <- dim(A)
              for(n in 1:dims[2]){
                  for(m in 1:dims[2]){
                      TFuns[n,m,] <-  - A[n,m,] / A[n,n,]
                  }
              }
              return(TFuns)
}

#' Calculate spectral matrix
#'
#' Calculates spectral matrix from the noise transfer function and the covariance
#' matrix
#' 
#' @param H a noise transfer function
#' @param sigma a covariance matrix 
#'
#' @return The spectral matrix
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'    
#' @keywords internal
#'    
#'       
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetSpectra <- function(H, sigma){
              S <- array(0, dim(H))
              for(n in 1:dim(H)[3]){
                S[,,n] <- H[,,n] %*% sigma %*% t(Conj(H[,,n]))
              }
              return(S)
}




#' Calculate open-loop transfer functions
#'
#' Calculates open-loop transfer functions from the spectral matrix
#' 
#' @param S a spectral matrix
#' @param use.cross boolean; use cross-spectrum to compute the transfer functions.
#'                  default is TRUE 
#'
#' @return Open-loop transfer functions
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'     
#' @keywords internal 
#'     
#'     
#'       
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetOpenTFuns <- function(S, use.cross = TRUE){
  OpenTFuns <- array(0, dim(S))
  dims <- dim(S)
  for(n in 1:dims[2]){
    for(m in 1:dims[2]){
      if(use.cross){
        OpenTFuns[n,m,] <-  S[m,n,] / abs(S[m,m,])
      } else {
        OpenTFuns[m,n,] <-  sqrt(abs(S[m,m,]) / abs(S[n,n,]))
      }
    }
  }
  return(OpenTFuns)
}


#' Calculate open-vs-closed-loop estimation difference
#'
#' 
#' @param open Open-loop transfer functions
#' @param closed Closed loop transfer functions 
#'
#' @return Estimation difference between open and closed loop transfer functions
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'    
#' @keywords internal
#'    
#'      
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetOpenVSClosedDif <- function(open, closed){
              pers <- array(0, dim(open))
              dims <- dim(open)
              for(n in 1:dims[3]){
                  pers[,,n] <- abs(abs(open[,,n]) - abs(closed[,,n])) * 100 /
                     abs(open[,,n])
              }
              return(pers)
}


#' Calculate transfer function B from transfer function A
#'
#' 
#' @param A transfer function A
#' @param a0 an A0 function
#'
#' @return Transfer function B
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'    
#'    
#' @keywords internal
#'      
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetMatrixBfromA <- function(A, a0){
              B <- array(0, dim(A))
              for(n in 1:dim(B)[3]){
                  B[ , , n] <- diag(c(1,1), dim(B)[1]) - A[ , , n]
                  #A[ , , n] <- a0 - B[ , , n]
              } 
              return(B)
}

#' Calculate transfer function A from noise transfer function
#'
#' 
#' @param H noise transfer function
#' 
#' @return Transfer function A
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'     
#' @keywords internal
#'     
#'        
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' @examples
#' #ADD EXAMPLE
GetMatrixAfromH <- function(H){
              A <- array(0, dim(H))
              for(n in 1:dim(A)[3]){
                  A[ , , n] <- solve(H[ , , n])
              } 
              return(A)
}                        
                       



# Developed by Alvaro Chao-Ecija
#
# This function calculates the frequency domain components and transfer functions
# from the time domain MVAR model. For compatibility purposes, the frequency support
# and the covariance matrix are estimated such as other packages do, like package
# grangers (see references).
#
# The functions are estimated taking into account the following representations 
# of the model:
#
# X = BX + W
# AX = W
# X = HW
#
# So that:
#
# S = H * sigma * t(Conj(H))
# A * S * t(Conj(A)) = sigma
#
# By allowing immediate effects from SBP to RR, the baroreflex transfer function
# is calculated as:
#
# b21/(1-b22) = -a21/a22 = h21/h11
#
# References ---------------------------------------------------------------------
#
# Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
# Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#
# Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
# Carrault. Linear multivariate models for physiological signal analysis: applications. 
# Comput Methods Programs Biomed. 1996;51(1-2):121-30
#
#
# Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
# in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
# Communications and Software. 2011
#
# Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
# J Clin Monit Comput. 2006;20(1):101-8
#
# CITA GRANGERS!!




