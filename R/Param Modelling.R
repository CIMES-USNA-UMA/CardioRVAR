
 




ParamFreqModel <- function(model, len = 1000, dt = 0.25, A0 = TRUE, sigma = NULL,
   coefs = NULL){
  #' Title Aqui titulo
  #'
  #' @param model Este es el modelo
  #' @param len Esto es el tamano del vector de frecuencias
  #' @param dt Hola
  #' @param A0 
  #' @param sigma 
  #' @param coefs 
  #'
  #' @return
  #' @export
  #'
  #' @examples
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






GetA0Fun <- function(sigma){
            b <- t(chol(sigma))
            D <- diag(b)
            b <- b %*% diag(1/D, nrow = length(D))
            a0 <- solve(b)
            new_sigma <- diag(D^2, nrow = length(D))
            return(list(a0 = a0, sigma = new_sigma))
}

UpdateWithA0 <- function(A0, coefs){
            ncoefs <- coefs
            for(n in 1:length(coefs)){
                ncoefs[[n]] <- A0 %*% coefs[[n]]
            }  
            return(ncoefs)
}



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




GetMatrixBfromVAR <- function(var, coefs, freqs, dt = 0.25){
              K <- var$K
              lags <- 1:var$p
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



GetMatrixAfromB <- function(B){
              A <- array(0, dim(B))
              for(n in 1:dim(B)[3]){
                  A[ , , n] <- diag(c(1,1), dim(B)[1]) - B[ , , n]
              } 
              return(A)
}


GetMatrixHfromA <- function(A){
              H <- array(0, dim(A))
              for(n in 1:dim(A)[3]){
                  H[ , , n] <- solve(A[ , , n])
              } 
              return(H)
}

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


GetSpectra <- function(H, sigma){
              S <- array(0, dim(H))
              for(n in 1:dim(H)[3]){
                S[,,n] <- H[,,n] %*% sigma %*% t(Conj(H[,,n]))
              }
              return(S)
}

GetOpenTFuns2 <- function(S){
              OpenTFuns <- array(0, dim(S))
              dims <- dim(S)
              for(n in 1:dims[2]){
                  for(m in 1:dims[2]){
                      OpenTFuns[m,n,] <-  sqrt(abs(S[m,m,]) / abs(S[n,n,]))
                  }
              }
              return(OpenTFuns)
}

# New version of GetOpenTFuns:

GetOpenTFuns <- function(S){
  OpenTFuns <- array(0, dim(S))
  dims <- dim(S)
  for(n in 1:dims[2]){
    for(m in 1:dims[2]){
      OpenTFuns[n,m,] <-  S[m,n,] / abs(S[m,m,])
    }
  }
  return(OpenTFuns)
}



GetOpenVSClosedDif <- function(open, closed){
              pers <- array(0, dim(open))
              dims <- dim(open)
              for(n in 1:dims[3]){
                  pers[,,n] <- abs(abs(open[,,n]) - abs(closed[,,n])) * 100 /
                     abs(open[,,n])
              }
              return(pers)
}
                            
GetMatrixBfromA <- function(A, a0){
              B <- array(0, dim(A))
              for(n in 1:dim(B)[3]){
                  B[ , , n] <- diag(c(1,1), dim(B)[1]) - A[ , , n]
                  #A[ , , n] <- a0 - B[ , , n]
              } 
              return(B)
}


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




