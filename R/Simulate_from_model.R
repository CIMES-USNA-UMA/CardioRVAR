




#' Simulate frequency domain model
#' 
#' Simulates the behavior of a frequency domain model after choosing several variables
#'
#' @param SM A computed frequency domain model
#' @param noises Which noise spectra to be included in the simulation
#' @param a0 Which immediate effects to be included in the simulation
#' @param xlim Specific limits for the frequency axis. Default is NULL
#' @param HFcolor Color for the HF band. Default is yellow
#' @param LFcolor color for the LF band. Default is green
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param show.cols boolean, show colors in the plot. Default is TRUE
#' @param phase.col color for the phase difference plot. Default is red
#' 
#'
#' @return A modified version of the original model according to the chosen parameters
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' new_model <- SimulateWithModel(freq_model, c(2,3), a0 = 2, f = 4)
SimulateWithModel <- function(SM, noises, a0, xlim = NULL, HFcolor = "yellow", LFcolor = "green", VLF = 0.04,
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
  return(S)
}


#' Plot simulated model
#' 
#' Generates a graphical representation of the noise transfer functions
#'
#' @param SM The computed frequency domain model from the system to be analyzed
#' @param S Simulated spectral matrix from the model to be analyzed
#' @param index Numeric index to specify a variable from the model
#' @param unit Unit of the output variable. Default is ms2/Hz
#' @param xlim Specific limits for the frequency axis. Default is NULL
#' @param HFcolor Color for the HF band. Default is yellow
#' @param LFcolor Color for the LF band. Default is green
#' @param VLF Maximum limit for the VLF band. Default is 0.04
#' @param LF Maximum limit for the LF band. Default is 0.15
#' @param HF Maximum limit for the HF band. Default is 0.4
#' @param show.cols Boolean, show colors in the plot. Default is TRUE
#' @param phase.col Color for the phase difference plot. Default is red
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
#' new_model <- SimulateWithModel(freq_model, c(2,3), a0 = 2, f = 4)
#' PlotSimulatedS(freq_model, new_model, 2)
#' 
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

#' Compute and plot spectral density
#' 
#' Computes and plots spectral density estimates from a certain variable
#'
#' @param x a vector with observations to be analyzed
#' @param p model order for the autoregressive model
#' @param xlim specific limits for the frequency axis. Default is NULL
#' @param HFcolor color for the HF band. Default is yellow
#' @param LFcolor color for the LF band. Default is green
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param show.cols boolean, show colors in the plot. Default is TRUE
#' @param phase.col color for the phase difference plot. Default is red
#' @param f sample rate of the vector of observations. Default is 4 Hz
#' @param plot boolean, plot the PSD function. Default is TRUE
#' @param output boolean, return the PSD function. Default is FALSE
#' @param n length of the vector of frequencies. Default is 1500
#' @param unit unit of the output variable. Default is ms2/Hz
#' 
#'
#' @return If the output argument is TRUE, a list with the computed power at
#'         LF and HF bands, as well as the peak frequencies, is returned.
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' RR <- DetrendedData[,"RR"]
#' PSD(RR, 21)
#' HRV <- PSD(RR, 21, plot = FALSE, output = TRUE)
#' HRV
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
    fLF <- freqs[(freqs <= LF) & (freqs >= VLF)]
    gHF <- gain[(freqs <= HF) & (freqs > LF)]
    gLF <- gain[(freqs <= LF) & (freqs >= VLF)]
    funHF <- splinefun(fHF, gHF, "monoH.FC")
    funLF <- splinefun(fLF, gLF, "monoH.FC")
    igHF <- integrate(funHF, lower = min(fHF), upper = max(fHF))
    igLF <- integrate(funLF, lower = min(fLF), upper = max(fLF))
    return(c(HF = igHF$value,
             LF = igLF$value,
             peakHF = freqs[match(max(gain[(freqs <= HF) & (freqs > LF)]), gain)],
             peakLF = freqs[match(max(gain[(freqs <= LF) & (freqs >= VLF)]), gain)]))
  }
}

