

#' Calculate spectral coherence
#' 
#' Calculates the spectral coherence between two variables
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#'
#' @return The spectral coherence between the two variables
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' coherence <- CalculateCoherence(freq_model, 1, 2)
#' coherence
#' 
CalculateCoherence <- function(SM, index1, index2){
  coherence <- SM$Spectra[index1,index2,]  / sqrt(SM$Spectra[index1,index1,]*SM$Spectra[index2,index2,])
  return(coherence)
}


#' Plot spectral coherence
#' 
#' Plots the spectral coherence between two variables
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#' @param coherence a vector with the spectral coherence values
#' @param thr a specific coherence threshold. Default is 0.5
#' @param xlim specific limits for the frequency axis. Default is NULL
#' @param HFcolor color for the HF band. Default is yellow
#' @param LFcolor color for the LF band. Default is green
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param show.cols boolean, show colors in the plot. Default is TRUE
#' @param thr.col color for the coherence threshold. Default is red
#' @param phase.col color for the phase difference plot. Default is red
#' 
#' 
#'
#' @return None
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' coherence <- CalculateCoherence(freq_model, 1, 2)
#' PlotCoherence(freq_model, 1, 2, coherence)
#' 
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



  
  
