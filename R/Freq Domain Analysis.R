

#' Compute expected values
#' 
#' Compute expected values from a frequency domain model of a particular system
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param str boolean, if TRUE the output will be the structure of the returned R object.
#'            Default is FALSE
#' @param weight applies gaussian-weighting function to the estimates before
#'               calculating the expected values. Default is TRUE
#' @param use.coh boolean, if TRUE the coherence and coherence threshold will be
#'                used for the computation of expected values. Default is TRUE
#' @param thr a coherence threshold to ensure the validity of the estimates. Default is 0.5
#' @param coherence the spectral coherence between the analyzed variables            
#'
#' @return A list containing the expected values computed at the HF and LF bands; or
#'         the structure of this list.
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#' @references 
#' McLoone V, Ringwood JV. A system identification approach to baroreflex sensitivity
#' estimation. IET Irish Signals and Systems Conference (ISSC 2012); 2012. 
#' p. 1-6. doi: 10.1049/ic.2012.0219.
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' GetExpectedVals(freq_model, use.coh = FALSE)
#' ExpectedVals <- GetExpectedVals(freq_model, use.coh = FALSE, str = FALSE)
#' ExpectedVals
#' 
#' # If the coherence is to be used:
#' coherence <- CalculateCoherence(freq_model, 1, 2)
#' ExpectedVals <- GetExpectedVals(freq_model, weight = FALSE, coherence = coherence,
#' str = FALSE)
#' ExpectedVals
#' 
GetExpectedValues <- function(SM, VLF = 0.04, LF = 0.15, HF = 0.4, 
    str = TRUE, weight = TRUE, use.coh = TRUE, thr = 0.5, coherence){
    frequency <- SM$Freqs 
    F <- NROW(frequency)
    LF_f <- frequency[frequency >= VLF]
    LF_f <- LF_f[LF_f < LF]
    HF_f <- frequency[frequency >= LF]
    HF_f <- HF_f[HF_f <= HF] 
    LF_band <- match(min(LF_f), frequency):match(max(LF_f), frequency)
    HF_band <- match(min(HF_f), frequency):match(max(HF_f), frequency)
    if(use.coh){
      for(n in 2:7){
        for(m in 1:2){
          for(p in 1:2){
              SM[[n]][m,p, c(1:NROW(coherence))[(abs(coherence)^2) < thr]] <- NA
          }
        }
      }
    }
    if(weight){
       supLF <- seq(-1, 1, len = NROW(LF_band))
       supHF <- seq(-1, 1, len = NROW(HF_band))
       wLF <- dnorm(supLF, sd = sd(supLF))
       wHF <- dnorm(supHF, sd = sd(supHF))
    } else {
       wLF <- rep(1, NROW(LF_band)) 
       wHF <- rep(1, NROW(HF_band)) 
    }
    SM_LF <- SM_HF <- SM
    for(n in 2:7){
        SM_LF[[n]] <- SM_HF[[n]] <- array(0, dim = c(2,2,1))
        if( n == 7){
           for(m in 1:2){
               for(p in 1:2){
                   SM_LF[[n]][m, p,] = sum(abs(SM[[n]][m, p,][LF_band]), na.rm = TRUE)
                   SM_HF[[n]][m, p,] = sum(abs(SM[[n]][m, p,][HF_band]), na.rm = TRUE)
               }
           }
        } else {
           for(m in 1:2){
               for(p in 1:2){
                   SM_LF[[n]][m, p,] = weighted.mean(abs(SM[[n]][m, p,][LF_band]),
                     w = wLF, na.rm = TRUE)
                   SM_HF[[n]][m, p,] = weighted.mean(abs(SM[[n]][m, p,][HF_band]),
                     w = wHF, na.rm = TRUE)
               }
           }
        }
    }
    SM_LF$Freqs <- NULL
    SM_HF$Freqs <- NULL
    if(str){
       return(str(list(LF = SM_LF, HF = SM_HF)))
    } else {
       return(list(LF = SM_LF, HF = SM_HF))
    }
} 





#' Get estimates at maximum coherence
#' 
#' Use the coherence between two variables to get estimates of the frequency model
#' at each band at maximum coherence levels
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param str boolean, if TRUE the output will be the structure of the returned R object.
#'            Default is FALSE
#' @param coherence the spectral coherence between the analyzed variables            
#'
#' @return A list containing the estimates at the HF and LF bands; or
#'         the structure of this list.
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' coherence <- CalculateCoherence(freq_model, 1, 2)
#' GetEstimateAtMaxCoh(freq_model, coherence = coherence)
#' Estimates <- GetEstimateAtMaxCoh(freq_model, coherence = coherence,
#' str = FALSE)
#' Estimates
#' 
GetEstimateAtMaxCoh <- function(SM, VLF = 0.04, LF = 0.15, HF = 0.4, 
                              str = TRUE, coherence){
  frequency <- SM$Freqs 
  F <- NROW(frequency)
  LF_f <- frequency[frequency >= VLF]
  LF_f <- LF_f[LF_f < LF]
  HF_f <- frequency[frequency >= LF]
  HF_f <- HF_f[HF_f <= HF] 
  LF_band <- match(min(LF_f), frequency):match(max(LF_f), frequency)
  HF_band <- match(min(HF_f), frequency):match(max(HF_f), frequency)
  SM_LF <- SM_HF <- SM
  for(n in 2:7){
    SM_LF[[n]] <- SM_HF[[n]] <- array(0, dim = c(2,2,1))
      for(m in 1:2){
        for(p in 1:2){
          LFg <- abs(SM[[n]][m, p,][LF_band]) 
          HFg <- abs(SM[[n]][m, p,][HF_band]) 
          SM_LF[[n]][m, p,] <- LFg[match(max(abs(coherence[LF_band])^2), abs(coherence[LF_band])^2)]
          SM_HF[[n]][m, p,] <- HFg[match(max(abs(coherence[HF_band])^2), abs(coherence[HF_band])^2)]
        }
      }
   
  }
  SM_LF$Freqs <- NULL
  SM_HF$Freqs <- NULL
  if(str){
    return(str(list(LF = SM_LF, HF = SM_HF)))
  } else {
    return(list(LF = SM_LF, HF = SM_HF))
  }
} 


#' Get peaks for each band
#' 
#' Get the frequency values in which the maximum peaks are located for each band
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param str boolean, if TRUE the output will be the structure of the returned R object.
#'            Default is FALSE
#'
#' @return A list containing the estimates at the HF and LF bands; or
#'         the structure of this list.
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' GetPeaks(freq_model)
#' Peaks <- GetEstimateAtMaxCoh(freq_model, str = FALSE)
#' Peaks
#' 
GetPeaks <- function(SM, VLF = 0.04, LF = 0.15, HF = 0.4, 
                                 str = TRUE){
  frequency <- SM$Freqs 
  F <- NROW(frequency)
  LF_f <- frequency[frequency >= VLF]
  LF_f <- LF_f[LF_f < LF]
  HF_f <- frequency[frequency >= LF]
  HF_f <- HF_f[HF_f <= HF] 
  LF_band <- match(min(LF_f), frequency):match(max(LF_f), frequency)
  HF_band <- match(min(HF_f), frequency):match(max(HF_f), frequency)
  SM_LF <- SM_HF <- SM
  for(n in 2:7){
    SM_LF[[n]] <- SM_HF[[n]] <- array(0, dim = c(2,2,1))
    for(m in 1:2){
      for(p in 1:2){
        SM_LF[[n]][m, p,] = LF_f[match(max(abs(SM[[n]][m, p,][LF_band])), abs(SM[[n]][m, p,][LF_band]))]
        SM_HF[[n]][m, p,] = HF_f[match(max(abs(SM[[n]][m, p,][HF_band])), abs(SM[[n]][m, p,][HF_band]))]
      }
    }
    
  }
  SM_LF$Freqs <- NULL
  SM_HF$Freqs <- NULL
  if(str){
    return(str(list(LF = SM_LF, HF = SM_HF)))
  } else {
    return(list(LF = SM_LF, HF = SM_HF))
  }
}



