#' Parametric frequency model estimation
#'
#' Estimates frequency model parameters from a time domain VAR model.
#' @param model A time domain VAR model
#' @param len Length of the vector of frequencies. Default is 1000
#' @param dt Inverse of sample rate. Default is 0.25
#' @param A0 Boolean; add 0 effects. Default TRUE
#' @param sigma Include a specific covariance matrix; otherwise it is obtained
#'              from the VAR model. Default is NULL
#' @param coefs Include VAR coefficients; otherwise these are obtained from the
#'              VAR model. Default is NULL
#' @param freq_support Method to compute the frequency support vector, either "pgram" or "seq".
#'                     Default is "pgram".
#' @param use.cross Boolean; use cross-spectrum to compute the open transfer functions.
#'                  default is TRUE. If FALSE, an alpha-index is used to compute the
#'                  open transfer functions
#' @param preserve.signals Boolean, include the input and output signals into the model list. Default
#'               is TRUE
#' @param preserve.p Boolean, include VAR model order p into the model list. Default
#'               is TRUE
#'
#' @return A list with the estimated components of the frequency domain model:
#' \item{Freqs}{Vector of frequency values}
#' \item{Vars_Transfer_funs}{Frequency-domain represantations of the VAR model blocks} 
#' \item{Transfer_Functions}{The computed closed-loop baroreflex transfer functions}
#' \item{Open_Transfer_Functions}{The computed open-loop baroreflex transfer functions}
#' \item{Open_Transfer_Functions_Overestimation}{How the open loop condition overestimates the results}
#' \item{Noise_Transfer_Functions}{The computed transfer functions from the noise source model}
#' \item{Spectra}{The computed auto-spectra and cross-spectra of the variables}
#' \item{Noise_Spectra}{The spectra of the white noise sources}
#' \item{a0}{The computed immediate effects of the model}
#'         
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' @details 
#'This function calculates the frequency domain components and transfer functions
#' from the time domain MVAR model. The functions are estimated taking into account the following representations 
#' of the model:
#'
#' \eqn{X = BX + W}
#' \eqn{AX = W}
#' \eqn{X = HW}
#'
#' So that:
#'
#' \eqn{S = H * sigma * t(Conj(H))}
#' \eqn{A * S * t(Conj(A)) = sigma}
#'
#' By allowing immediate effects from SBP to RR, the baroreflex transfer function
#' is calculated as:
#'
#' \eqn{b_{21}/(1-b_{22}) = -a_{21}/a_{22} = h_{21}/h_{11}}
#' 
#' 
#' Some of the modules integrated in these functions have been designed so that the outputs
#' have compatible characteristics (and are thus reproducible) with the ones obtained from 
#' other packages such as \href{https://CRAN.R-project.org/package=grangers}{package grangers}
#' (e.g., the freq_support option available in this function allows said compatibility). For
#' more information, please check the references section.
#'         
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
#'
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
#' Matteo Farne and Angela Montanari (2019). grangers: Inference on Granger-Causality 
#' in the Frequency Domain. R package version 0.1.0. https://CRAN.R-project.org/package=grangers
#' 
#' @export
#'
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
ParamFreqModel <- function(model, len = 1000, dt = 0.25, A0 = TRUE, sigma = NULL,
   coefs = NULL, freq_support = c("pgram", "seq"), use.cross = TRUE, preserve.signals = TRUE,
   preserve.p = TRUE){
                   if(is.null(sigma) & is.null(coefs)){
                      coefs <- GetCoefs(model)
                      sigma <- summary(model)$cov * 2 * dt 
                   } else if(is.null(sigma) & !is.null(coefs)){
                             stop("Indicate coeficients")
                   } else if(is.null(sigma) & !is.null(coefs)){
                             stop("Indicate covariance matrix")
                   }
                   freq_support <- match.arg(freq_support)
                   freqs <- GetFreqs(freq_support, model, len, dt)
                   B <- GetMatrixBfromCoefs(coefs, freqs, dt) 
                   A <- GetMatrixAfromB(B)
                   H <- GetMatrixHfromA(A)
                   TFuns <- GetTransFuns(A)
                   #sigma <- sigma * 2
                   S <- GetSpectra(H, sigma)
                   OpenTFuns <- GetOpenTFuns(S, use.cross = use.cross)
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
                   if(preserve.signals & !is.null(model$Signals)) output$Signals <-model$Signals
                   if(preserve.p) output$p <- model$p
                   return(output)
}


#' Include A0 effects 
#'
#' Estimates and includes A0 effects for a frequency domain model
#' @param system The model of the system whose instantaneous effects we want to include
#'
#' @return A list with the estimated components of the frequency domain model:
#' \item{Freqs}{Vector of frequency values}
#' \item{Vars_Transfer_funs}{Frequency-domain represantations of the VAR model blocks} 
#' \item{Transfer_Functions}{The computed closed-loop baroreflex transfer functions}
#' \item{Open_Transfer_Functions}{The computed open-loop baroreflex transfer functions}
#' \item{Open_Transfer_Functions_Overestimation}{How the open loop condition overestimates the results}
#' \item{Noise_Transfer_Functions}{The computed transfer functions from the noise source model}
#' \item{Spectra}{The computed auto-spectra and cross-spectra of the variables}
#' \item{Noise_Spectra}{The spectra of the white noise sources}
#' \item{a0}{The computed immediate effects of the model}
#'
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#' @references 
#' Barbieri R, Parati G, Saul JP. Closed- versus Open-Loop Assessment of Heart Rate
#' Baroreflex. IEEE Eng Med Biol Mag. 2001;20(2):33-42
#'
#' Korhonen I, Mainardi L, Baselli G, Bianchi A, Loula P, 
#' Carrault. Linear multivariate models for physiological signal analysis: applications. 
#' Comput Methods Programs Biomed. 1996;51(1-2):121-30
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
                B <- array(0, dim(Vars_Transfer_fun))
                for(n in 1:(dim(Noise_Transfer_fun)[3])){
                    Noise_Transfer_fun[,,n] <- Noise_Transfer_fun[,,n] %*% 
                       solve(a0)
                }
                A <- GetMatrixAfromH(Noise_Transfer_fun)
                Vars_Transfer_fun <- GetMatrixBfromA(A)
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
#' @param sigma A covariance matrix
#'
#' @return A list which includes the A0 function as well as the new covariance matrix
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
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
#' @param A0 The A0 function
#' @param coefs Coefficients to be updated
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
UpdateWithA0 <- function(A0, coefs){
            ncoefs <- coefs
            ncoefs[[1]] <- diag(1, nrow(A0)) - A0
            for(n in 2:length(coefs)){
                ncoefs[[n]] <- A0 %*% coefs[[n]]
            }  
            return(ncoefs)
}


#' Get coefs from VAR model
#'
#' Extracts VAR model coefficients for further analysis
#' @param var A VAR model
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
GetCoefs <- function(var){
            max.lag = var$p
            K <- var$K
            lags <- 0:max.lag
            coefs <- coef(var)
            for(k in 1:K){
              coefs[[k]] <- rbind(matrix(0, nrow = K, ncol = 4), coefs[[k]])
            }
            sel_coefs <- NULL
            coef_list <- list()
            length(coef_list) <- max.lag + 1
            for(k in 1:K){
                sel_coefs <- cbind(sel_coefs, coefs[[k]][,1])
            }
            indices <- as.vector(rbind(lags, lags))
            coef_list <-
              lapply(lapply(split(sel_coefs, indices), matrix, ncol = K), t)
            return(coef_list)
}



#' Calculate transfer function B from VAR model
#'
#' Calculates transfer function B from a specific VAR model
#' @param var A var model
#' @param freqs A vector of frequencies 
#' @param dt Inverse sample rate. Default is 0.25
#' @param a0 Boolean. Should the coefficients of the VAR model be adjusted with 
#'           instantaneous interactions? Default is FALSE
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
GetMatrixBfromVAR <- function(var, freqs, dt = 0.25, a0 = FALSE){
              coefs <- GetCoefs(var)
              if((a0)){
                a0 <- GetA0Fun(summary(model)$cov * 2 * dt)$a0
                coefs <- UpdateWithA0(a0, coefs)
              }
              B <- GetMatrixBfromCoefs(coefs, freqs, dt)
              return(B)
}


#' Calculate transfer function B from VAR coefficients
#'
#' Calculates transfer function B from a specific VAR model coefficients
#' @param coefs Coefficients from a var model
#' @param freqs A vector of frequencies 
#' @param dt Inverse sample rate. Default is 0.25
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
GetMatrixBfromCoefs <- function(coefs, freqs, dt = 0.25){
              if(max(freqs) > 0.5) freqs <- freqs * dt
              K <- ncol(coefs[[1]])
              lags <- (1:length(coefs)) - 1
              B <- array(0, c(K, K, NROW(freqs)))
              for(n in 1:K){
                  for(m in 1:K){
                      for(lag in 1:length(coefs)){
                          B[n, m, ] <- B[n, m, ] +
                               coefs[[lag]][n, m] * exp(-2*pi*freqs*1i*lags[lag]) 
                      }
                  }
              }
              return(B)
}


#' Calculate transfer function B from transfer function A
#'
#' Calculates transfer function B from a previously calculated transfer function A
#' @param B Transfer function B
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
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
#' @param A Transfer function A
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
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
#' @param H A noise transfer function
#' @param sigma A covariance matrix 
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
#' @param S A spectral matrix
#' @param use.cross Boolean; use cross-spectrum to compute the transfer functions.
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
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
#' @param A Transfer function A
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
GetMatrixBfromA <- function(A){
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
#' Faes L, Nollo G. Multivariate Frequency Domain Analysis of Causal Interactions
#' in Physiological Time Series. In Biomedical Engineering, Trends in Electronics,
#' Communications and Software. 2011
#'
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8
#' 
GetMatrixAfromH <- function(H){
              A <- array(0, dim(H))
              for(n in 1:dim(A)[3]){
                  A[ , , n] <- solve(H[ , , n])
              } 
              return(A)
}                        

#' Get vector of frequencies
#'
#' 
#' @param type Method to get the vector of frequencies
#' @param model A time domain VAR model
#' @param len Length of the vector of frequencies
#' @param dt Inverse of sample rate
#'
#' @return Vector of frequency values
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'
#' @details This function allows to obtain vectors of frequency values by either creating a sequence
#' or by spectral estimation. This second option was introduced to make the results compatible (and thus
#' reproducible) with the ones obtained from other packages such as package grangers (for more information
#' regarding this package, please check the reference section).
#'    
#' @keywords internal
#'      
#' @references 
#' Matteo Farne and Angela Montanari (2019). grangers: Inference on Granger-Causality in the Frequency
#' Domain. R package version 0.1.0. https://CRAN.R-project.org/package=grangers
#' 
GetFreqs <- function(type, model, len, dt){
  if(type == "pgram"){
    # For compatibility and reproducibility with other packages such as grangers.
    freqs <- spec.pgram(model$y[,1], plot = FALSE, 
                        pad = len/(model$totobs/2) - 1)$freq
  } else {
    #freqs <- seq(0, 1/(2*dt), len = len)
    freqs <- seq(0, 1/(2), len = len)
  }
  return(freqs)
}                      








