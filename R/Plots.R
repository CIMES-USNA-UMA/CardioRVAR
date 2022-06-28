# Developed by Alvaro Chao-Ecija
#
# Several functions to plot the calculated transfer functions and spectra from the
# frecuency domain models.
#
#

PlotTransferFun <- function(SM, index1, index2, unit1 = "ms", unit2 = "mmHg",
                            xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                            LF = 0.15, HF = 0.4, show.cols = TRUE, phase.col = "red", open = FALSE,
                            plot.phase = TRUE){
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  if(open){
    tf <- SM$Open_Transfer_Functions
    tf <- Conj(tf)
  } else {
    tf <- SM$Transfer_Functions
  }
  tf <- tf[index1, index2, ]
  if(index1 == index2){
    stop("Index values must be different")
  }
  if(index1 > index2){
    unit <- paste(unit1, "/", unit2, sep = "")
  } else {
    unit <- paste(unit2, "/", unit1, sep = "")
  }
  gain <- abs(tf)
  if(open){
    phase <- GetOpenTransFunPhase(SM, index1, index2)
  } else{
  phase <- GetTransFunPhase(SM, index1, index2)
  }
  freqs <- SM$Freqs
  
  if(plot.phase)  par(mfrow = c(1,2))
  max_gain <- max(gain[(freqs <= HF) & (freqs >= VLF)])
  plot(freqs, gain, xlim = xlim, xlab = "Frequency",
       ylab = paste("Gain (", unit, ")", sep = ""),
       main = "Gain Transfer Function", type = "l",
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
  if(plot.phase){
    plot(freqs, phase, xlim = xlim, xlab = "Frequency",
         ylab = "Phase (rads)",
         main = "Phase Transfer Function",
         ylim = c(-pi, pi), type = "l")
    if(show.cols){
      abline(h = 0, col = phase.col)
      abline(h = pi/2, col = phase.col)
      abline(h = -pi/2, col = phase.col)
      abline(h = pi, col = phase.col)
      abline(h = -pi, col = phase.col)
    }
  }
  
  
}


PlotBlockTransferFun <- function(SM, index1, index2, unit1 = "ms", unit2 = "mmHg",
                            xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                            LF = 0.15, HF = 0.4, show.cols = TRUE, phase.col = "red",
                            plot.phase = TRUE){
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  tf <- SM$Vars_Transfer_funs
  tf <- tf[index1, index2, ]
  if(index1 > index2){
    unit <- paste(unit1, "/", unit2, sep = "")
  } else if (index1 < index2) {
    unit <- paste(unit2, "/", unit1, sep = "")
  } else {
    unit <- "adimensional"
  }
  gain <- abs(tf)
  phase <- atan2(Im(tf), Re(tf))
  freqs <- SM$Freqs
  if(plot.phase) par(mfrow = c(1,2))
  max_gain <- max(gain[(freqs <= HF) & (freqs >= VLF)])
  plot(freqs, gain, xlim = xlim, xlab = "Frequency",
       ylab = paste("Gain (", unit, ")", sep = ""),
       main = "Gain Transfer Function", type = "l",
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
  if(plot.phase){
    plot(freqs, phase, xlim = xlim, xlab = "Frequency",
         ylab = "Phase (rads)",
         main = "Phase Transfer Function",
         ylim = c(-pi, pi), type = "l")
    if(show.cols){
      abline(h = 0, col = phase.col)
      abline(h = pi/2, col = phase.col)
      abline(h = -pi/2, col = phase.col)
      abline(h = pi, col = phase.col)
      abline(h = -pi, col = phase.col)
    }
  }
  
}


PlotNoiseTransferFun <- function(SM, index1, index2, unit1 = "ms", unit2 = "mmHg",
                                 xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
                                 LF = 0.15, HF = 0.4, show.cols = TRUE, phase.col = "red",
                                 plot.phase = TRUE){
  if(is.null(xlim)){
    xlim = c(VLF, HF)
  }
  tf <- SM$Noise_Transfer_fun
  tf <- tf[index1, index2, ]
  if(index1 > index2){
    unit <- paste(unit1, "/", unit2, sep = "")
  } else if (index1 < index2) {
    unit <- paste(unit2, "/", unit1, sep = "")
  } else {
    unit <- "adimensional"
  }
  gain <- abs(tf)
  phase <- atan2(Im(tf), Re(tf))
  freqs <- SM$Freqs
  if(plot.phase) par(mfrow = c(1,2))
  max_gain <- max(gain[(freqs <= HF) & (freqs >= VLF)])
  plot(freqs, gain, xlim = xlim, xlab = "Frequency",
       ylab = paste("Gain (", unit, ")", sep = ""),
       main = "Gain Transfer Function", type = "l",
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
  if(plot.phase){
    plot(freqs, phase, xlim = xlim, xlab = "Frequency",
         ylab = "Phase (rads)",
         main = "Phase Transfer Function",
         ylim = c(-pi, pi), type = "l")
    if(show.cols){
      abline(h = 0, col = phase.col)
      abline(h = pi/2, col = phase.col)
      abline(h = -pi/2, col = phase.col)
      abline(h = pi, col = phase.col)
      abline(h = -pi, col = phase.col)
    }
  }
  
}


