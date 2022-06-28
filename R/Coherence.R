# Developed by Alvaro Chao-Ecija
#
# Functions for calculating and displaying the coherence between X and Y, 
# as well as the causal threshold.




CalculateCoherence <- function(SM, index1, index2){
  coherence <- SM$Spectra[index1,index2,]  / sqrt(SM$Spectra[index1,index1,]*SM$Spectra[index2,index2,])
  return(coherence)
}


PlotCoherence <- function(SM, index1, index2, coherence, thr = 0.5, xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                          LF = 0.15, HF = 0.4, show.cols = TRUE, thr.col = "red", phase.col = "red"){
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  gain <- abs(coherence)^2
  phase <- GetOpenTransFunPhase(SM, index1, index2)
  freqs <- SM$Freqs
  par(mfrow = c(1,2))
  max_gain <- 1
  plot(freqs, gain, xlim = xlim, xlab = "Frequency",
       ylab = "",
       main = "Magnitude Squared Coherence", type = "l",
       ylim = c(0, max_gain))
  if(show.cols){
    polygon(x = c(freqs[(freqs <= HF) & (freqs >= VLF)], 
                  rev(freqs[(freqs <= HF) & (freqs >= VLF)])),
            y = c(gain[(freqs <= HF) & (freqs >= VLF)],
                  double(NROW(gain[(freqs <= HF) & (freqs >= VLF)]))),
            col = HFcolor)
    polygon(x = c(freqs[(freqs <= LF) & (freqs >= VLF)], 
                  rev(freqs[(freqs <= LF) & (freqs >= VLF)])),
            y = c(gain[(freqs <= LF) & (freqs >= VLF)],
                  double(NROW(gain[(freqs <= LF) & (freqs >= VLF)]))),
            col = LFcolor)
    abline(h = thr, col = thr.col)
  }
  plot(freqs, phase, xlim = xlim, xlab = "Frequency",
       ylab = "Phase (rads)",
       main = "Phase Spectral Coherence",
       ylim = c(-pi, pi), type = "l")
  if(show.cols){
    abline(h = 0, col = phase.col)
    abline(h = pi/2, col = phase.col)
    abline(h = -pi/2, col = phase.col)
    abline(h = pi, col = phase.col)
    abline(h = -pi, col = phase.col)
  }
}



  
  
