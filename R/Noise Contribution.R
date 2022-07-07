# Developed by Alvaro Chao-Ecija

# References:
# Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
# J Clin Monit Comput. 2006;20(1):101-8
#


NoiseContribution <- function(SM, index1, index2, VLF = 0.04, LF = 0.15, HF = 0.4, use.coh  =TRUE,
                              thr = 0.5, coherence){
  frequency <- SM$Freqs 
  F <- NROW(frequency)
  LF_f <- frequency[frequency >= VLF]
  LF_f <- LF_f[LF_f < LF]
  HF_f <- frequency[frequency >= LF]
  HF_f <- HF_f[HF_f <= HF] 
  if(use.coh){
    for(n in 2:7){
      for(m in 1:2){
        for(p in 1:2){
          SM[[n]][m,p, c(1:NROW(coherence))[(abs(coherence)^2) < thr]] <- NA
        }
      }
    }
  }
  LF_band <- match(min(LF_f), frequency):match(max(LF_f), frequency)
  HF_band <- match(min(HF_f), frequency):match(max(HF_f), frequency)
  HF_cont <- sum((abs(SM$Noise_Transfer_fun[index1, index2, HF_band])^2) *
                   SM$Noise_Spectra[index2, index2], na.rm = TRUE) / sum(abs(SM$Spectra[index1,
                                                                          index1,
                                                                          HF_band]), na.rm = TRUE)
  LF_cont <- sum((abs(SM$Noise_Transfer_fun[index1, index2, LF_band])^2) *
                   SM$Noise_Spectra[index2, index2], na.rm = TRUE) / sum(abs(SM$Spectra[index1,
                                                                          index1,
                                                                          LF_band]), na.rm = TRUE)
  return(c(HF = HF_cont*100, LF = LF_cont*100))
}

PlotNoiseContribution <- function(cont1, cont2, label1 = "RR noise", label2 = "SBP noise"){
  par(mfrow = c(1,2))
  pie(c(cont1[1], cont2[1]), labels = c(label1, label2), main = "HF band")
  pie(c(cont1[2], cont2[2]), labels = c(label1, label2), main = "LF band")
}



PlotCausality <- function(SM, index, VLF = 0.04, LF = 0.15, HF = 0.4, xlim = NULL){
  freqs <- SM$Freqs 
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  G21 <- (abs(SM$Noise_Transfer_fun[1, 2, ])^2) *
    SM$Noise_Spectra[2, 2] / abs(SM$Spectra[1,1,])
  G12 <- (abs(SM$Noise_Transfer_fun[2, 1, ])^2) *
    SM$Noise_Spectra[1, 1] / abs(SM$Spectra[2,2,])
  Max <- max(max(G12[(freqs <= HF) & (freqs >= VLF)]), max(G21[(freqs <= HF) & (freqs >= VLF)]))
  G = list(G12, G21)
  if(index == 1) index1 <- 2
  if(index == 2) index1 <- 1
  plot(freqs, G[[index]], xlim = xlim, xlab = "Frequency",
       ylab = "Causality",
       main = "Causal flow (current branch in red)", type = "l",
       ylim = c(0, Max), col = "red")
  lines(freqs, G[[index1]])
  abline(v = VLF, col = "grey")
  abline(v = LF, col = "grey")
  abline(v = HF, col = "grey")
}






  