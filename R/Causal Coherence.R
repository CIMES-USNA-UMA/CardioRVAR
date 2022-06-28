# Developed by Alvaro Chao-Ecija
#
# These functions have been designed for the purpose of calculating and displaying
# the causal coherence from X to Y and Y to X. The calculation is as indicated
# by Faes et al. 
#
# References: ------------------------------------------------------------
#
# Faes L, Porta A, Antolini R, Nollo G. Role of causality in the evaluation of coherence and 
# transfer function between heart period and systolic pressure in humans. Comput Cardiol.
# 2004;277-80.
#
# Faes L, Mase M, Nollo G, Chon KH, Florian JP. Measuring postural-related changes of spontaneous
# baroreflex sensitivity after repeated long-duration diving: frequency domain approaches. Auton Neurosci. 
# 2013 Nov;178(1-2):96-102.


CalculateCausalPhase <- function(ccoh){
  N <- length(ccoh)
  for(n in 1:N){
    ccoh[[n]] <- atan2(Im(ccoh[[n]]), Re(ccoh[[n]]))
  }
  return(ccoh)
}


CalculateCausalCoherence <- function(SM){
  S <- SM$Spectra
  H <- SM$Noise_Transfer_fun
  sigma <- SM$Noise_Spectra
  C1 <- abs(sigma[2,2]*abs(H[1,2,]^2)/S[1,1,]) #H21 = 0, Causalidad de S2 a S1 
  C2 <- abs(sigma[1,1]*abs(H[2,1,]^2)/S[2,2,]) #H12 = 0, Causalidad de S1 a S2
  # S2 = Input (RR)
  # S1 = Output (SBP)
  C3 <- abs(CalculateCoherence(SM, 1, 2))^2
  return(list(C1 = C1, C2 = C2, Cr = C3))
}




PlotCausalCoherence <- function(SM, CCoh, VLF = 0.04, LF = 0.15, HF = 0.4, xlim = NULL){
  freqs <- SM$Freqs 
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  if(length(CCoh) > 3){
  Max <- max(max(CCoh$C1[(freqs <= HF) & (freqs >= VLF)]), max(CCoh$C2[(freqs <= HF) & (freqs >= VLF)]),
             max(CCoh$C3[(freqs <= HF) & (freqs >= VLF)]), max(CCoh$Cr[(freqs <= HF) & (freqs >= VLF)]))
  } else {
    Max <- max(max(CCoh$C1[(freqs <= HF) & (freqs >= VLF)]), max(CCoh$C2[(freqs <= HF) & (freqs >= VLF)]), 
               max(CCoh$Cr[(freqs <= HF) & (freqs >= VLF)]))
  }
  plot(freqs, CCoh$C1, xlab = "Frequency (Hz)", ylab = "Causal Coherence", col = "red", ylim = c(0, Max), xlim = xlim,
       type = "p", pch = 4)
  points(freqs, CCoh$C2, col = "blue", pch = 2)
  lines(freqs, CCoh$Cr, col = "green")
  if(length(CCoh) > 3) points(freqs, CCoh$C3, col = "green", pch = 4)
  abline(v = VLF, col = "gray")
  abline(v = LF, col = "grey")
  abline(v = HF, col = "grey")
  
}
  

GetMeanCoherence <- function(SM, coherence, HF = 0.4, LF = 0.15, VLF = 0.04,
                             weight = TRUE){
  frequency <- SM$Freqs 
  F <- NROW(frequency)
  LF_f <- frequency[frequency >= VLF]
  LF_f <- LF_f[LF_f < LF]
  HF_f <- frequency[frequency >= LF]
  HF_f <- HF_f[HF_f <= HF] 
  LF_band <- match(min(LF_f), frequency):match(max(LF_f), frequency)
  HF_band <- match(min(HF_f), frequency):match(max(HF_f), frequency)
  HF <- coherence
  LF <- coherence
  if(weight){
    supLF <- seq(-1, 1, len = NROW(LF_band))
    supHF <- seq(-1, 1, len = NROW(HF_band))
    wLF <- dnorm(supLF, sd = sd(supLF))
    wHF <- dnorm(supHF, sd = sd(supHF))
  } else {
    wLF <- rep(1, NROW(LF_band)) 
    wHF <- rep(1, NROW(HF_band)) 
  }
  for(n in 1:length(HF)){
    HF[[n]] <- mean(coherence[[n]][HF_band])
    LF[[n]] <- mean(coherence[[n]][LF_band])
    HF[[n]] <- weighted.mean(coherence[[n]][HF_band],
                  w = wHF, na.rm = TRUE)
    LF[[n]] <- weighted.mean(coherence[[n]][LF_band],
                             w = wLF, na.rm = TRUE)
  }
  return(list(HF = HF, LF = LF))
}

GetMaxCoherence <- function(SM, coherence, HF = 0.4, LF = 0.15, VLF = 0.04){
  frequency <- SM$Freqs 
  F <- NROW(frequency)
  LF_f <- frequency[frequency >= VLF]
  LF_f <- LF_f[LF_f < LF]
  HF_f <- frequency[frequency >= LF]
  HF_f <- HF_f[HF_f <= HF] 
  LF_band <- match(min(LF_f), frequency):match(max(LF_f), frequency)
  HF_band <- match(min(HF_f), frequency):match(max(HF_f), frequency)
  HF <- coherence
  LF <- coherence
  for(n in 1:length(HF)){
    HF[[n]] <- max(coherence[[n]][HF_band])
    LF[[n]] <- max(coherence[[n]][LF_band])
  }
  return(list(HF = HF, LF = LF))
}




