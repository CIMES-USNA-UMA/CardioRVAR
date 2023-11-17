
#' Noise source contribution
#' 
#' Estimates the noise source contribution for the variables of the analyzed system
#'
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param use.coh boolean, if TRUE the coherence and coherence threshold will be
#'                used for the computation of expected values. Default is TRUE
#' @param thr a specific coherence threshold. Default is 0.5
#' @param coherence a vector with the spectral coherence values 
#' @param print.flow Boolean. If TRUE, prints a message describing the analysed
#'                   contribution
#' 
#'
#' @return None
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#' @references
#' Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
#' J Clin Monit Comput. 2006;20(1):101-8    
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' noise_con <- NoiseContribution(freq_model, 1, 2, use.coh = FALSE)
#' noise_con
#' 
#' # The coherence can be used as a method for identifying reliable estimates
#' coherence <- CalculateCoherence(freq_model, 1, 2)
#' noise_con_thr <- NoiseContribution(freq_model, 1, 2, coherence = coherence)
#' noise_con_thr
#' 
NoiseContribution <- function(SM, index1, index2, VLF = 0.04, LF = 0.15, HF = 0.4, use.coh  =TRUE,
                              thr = 0.5, coherence, print.flow = FALSE){
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
  if(print.flow){
    variable_names <- colnames(SM$a0)
    variable_S <- variable_names[index1]
    variable_Noise <- variable_names[index2]
    print(paste("Noise source contribution of", variable_Noise, 
                "to the spectrum of", variable_S))
  }
  return(c(HF = HF_cont*100, LF = LF_cont*100))
}




#' Plot noise source constribution
#' 
#' Generates pie plots showing the estimated noise source contribution from two noise sources
#'
#' @param cont1 estimated contribution from one noise source
#' @param cont2 estimated contribution from a second noise source
#' @param label1 label for the first noise source. Default is "RR noise"
#' @param label2 label for the second noise source. Default is "SBP noise"
#' 
#' @return None
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' 
#' noise_con1 <- NoiseContribution(freq_model, 1, 2, use.coh = FALSE, print.flow = TRUE)
#' noise_con2 <- NoiseContribution(freq_model, 1, 1, use.coh = FALSE, print.flow = TRUE)
#' noise_con1
#' noise_con2
#' noise_con1 + noise_con2
#' PlotNoiseContribution(noise_con1, noise_con2, label1 = "RR noise", label2 = "SBP noise")
#' 
PlotNoiseContribution <- function(cont1, cont2, label1 = "RR noise", label2 = "SBP noise"){
  par(mfrow = c(1,2))
  pie(c(cont1[1], cont2[1]), labels = c(label1, label2), main = "HF band")
  pie(c(cont1[2], cont2[2]), labels = c(label1, label2), main = "LF band")
}


#' Plot causal flow
#' 
#' Generates frequency domain plots indicating the behavior of the causal flow of a certain branch
#' 
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index an numeric value indicating which branch should be highlited
#' @param VLF maximum limit for the VLF band. Default is 0.04
#' @param LF maximum limit for the LF band. Default is 0.15
#' @param HF maximum limit for the HF band. Default is 0.4
#' @param xlim specific limits for the frequency axis. Default is NULL
#' @return None
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#' @export
#' 
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' PlotCausality(freq_model, 1)
#' PlotCausality(freq_model, 2)
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






  