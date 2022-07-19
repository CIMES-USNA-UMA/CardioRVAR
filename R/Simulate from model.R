





SimulateWithModel <- function(SM, noises, a0, f, xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                              LF = 0.15, HF = 0.4, show.cols = TRUE, phase.col = "red" ){
  sigma <- diag(noises)
  H <- SM$Noise_Transfer_fun
  N = dim(H)[3]
  b <- diag(rep(1,dim(H)[1]))
  b[lower.tri(b)] <- a0
  a0 <- b
  S <- H
  for(n in 1:N){
    H[,,n] <- (H[,,n] %*% SM$a0) %*% solve(a0)
    S[,,n] <- H[,,n] %*% sigma %*% Conj(t(H[,,n]))
  }
  return(S/f)
}

PlotSimulatedS <- function(SM, S, index, unit = "ms2/Hz",
                           xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                           LF = 0.15, HF = 0.4, show.cols = TRUE, phase.col = "red"){
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  S <- S[index, index, ]
  gain <- abs(S) 
  freqs <- SM$Freqs
  max_gain <- max(gain[(freqs <= HF) & (freqs >= VLF)])
  plot(freqs, gain, xlim = xlim, xlab = "Frequency",
       ylab = unit,
       main = "PSD", type = "l",
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
  }
}


PSD <- function(x, p, 
                           xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                           LF = 0.15, HF = 0.4, show.cols = TRUE, phase.col = "red", f = 4, 
                    plot = TRUE, output = FALSE, n = 1500, unit = "ms2/Hz"){
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  S <- spec.ar(ts(x, start = 0, frequency = f), plot = FALSE, order = p, n.freq = n)
  gain <- S$spec * 2
  freqs <- S$freq
  gain <- gain[(freqs <= HF) & (freqs >= VLF)]
  freqs <- freqs[(freqs <= HF) & (freqs >= VLF)]
  n_freqs <- NROW(freqs)
  max_gain <- max(gain)
  if(plot){
  plot(freqs, gain, xlim = xlim, xlab = "Frequency",
       ylab = unit,
       main = "PSD", type = "l",
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
  }
  }
  if(output){
    fHF <- freqs[(freqs <= HF) & (freqs > LF)]
    fLF <- freqs[(freqs <= HF) & (freqs > LF)]
    HF <- gain[(freqs <= HF) & (freqs > LF)]
    LF <- gain[(freqs <= LF) & (freqs >= VLF)]
    funHF <- splinefun(fHF, HF, "monoH.FC")
    funLF <- splinefun(fLF, LF, "monoH.FC")
    gHF <- integrate(funHF, lower = min(fHF), upper = max(fHF))
    gLF <- integrate(funLF, lower = min(fLF), upper = max(fLF))
    return(c(HF = gHF,
             LF = gLF,
             peakHF = freqs[match(max(gain[(freqs <= HF) & (freqs > LF)]), gain)],
             peakLF = freqs[match(max(gain[(freqs <= LF) & (freqs >= VLF)]), gain)]))
  }
}

