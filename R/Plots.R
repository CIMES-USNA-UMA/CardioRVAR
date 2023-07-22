#' Plot transfer functions
#' 
#' Generates a graphical representation of the closed-loop and open-loop transfer
#' functions
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#' @param unit1 unit of the output variable. Default is ms
#' @param unit2 unit of the input variable. Default is mmHg
#' @param xlim specific limits for the frequency axis. Default is NULL
#' @param HFcolor color for the HF band. Default is yellow
#' @param LFcolor color for the LF band. Default is green
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param show.cols boolean, show colors in the plot. Default is TRUE
#' @param phase.col color for the phase difference plot. Default is red
#' @param open boolean, show open-loop transfer functions. Default is FALSE
#' @param plot.phase boolean, plot phase difference calculated from the transfer function.
#'                   Default is TRUE
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
#' PlotTransferFun(freq_model, 1, 2) # Closed-loop
#' PlotTransferFun(freq_model, 1, 2, open  = TRUE) # Open-loop
#' 
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


#' Plot noise transfer functions
#' 
#' Generates a graphical representation of the noise transfer functions
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#' @param unit1 unit of the output variable. Default is ms
#' @param unit2 unit of the input variable. Default is mmHg
#' @param xlim specific limits for the frequency axis. Default is NULL
#' @param HFcolor color for the HF band. Default is yellow
#' @param LFcolor color for the LF band. Default is green
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param show.cols boolean, show colors in the plot. Default is TRUE
#' @param phase.col color for the phase difference plot. Default is red
#' @param plot.phase boolean, plot phase difference calculated from the transfer function.
#'                   Default is TRUE
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
#' PlotNoiseTransferFun(freq_model, 1, 2)
#' 
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


