# Developed by Alvaro Chao-Ecija
# 
# These are several functions to calculate estimates from each frequency
# domain function. These functions can be weighted using a Gaussian functions,
# according to the procedure described by


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
       #tLF <- 1:NROW(LF_band)
       #tHF <- 1:NROW(HF_band)
       #wLF <- exp(-((tLF - mean(tLF))^2) / (2 * var(tLF)))
       #wHF <- exp(-((tHF - mean(tHF))^2) / (2 * var(tHF)))
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

# References

