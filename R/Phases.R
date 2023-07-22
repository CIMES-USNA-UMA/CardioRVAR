

#' Calculate closed-loop phase transfer function
#'
#' Estimates the phase transfer function out of the closed-loop transfer functions
#' of a particular system
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#'
#' @return A vector containing the phase transfer function obtained from the analyzed
#'         closed-loop transfer function
#'         
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'
#' @export
#'
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' Phase <- GetTransFunPhase(freq_model, 1, 2)
#' Phase[1:100]
GetTransFunPhase <- function(SM, index1, index2){
  TF <- SM$Transfer_Functions[index1, index2,]
  phase <- atan2(Im(TF), Re(TF))
  return(phase)
}


#' Calculate open-loop phase transfer function
#'
#' Estimates the phase transfer function out of the open-loop transfer functions
#' of a particular system
#' @param SM the computed frequency domain model from the system to be analyzed
#' @param index1 numeric index to specify a variable from the model
#' @param index2 numeric index to specify a variable from the model
#'
#' @return A vector containing the phase transfer function obtained from the analyzed
#'         open-loop transfer function
#'         
#' @author Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
#'         
#'
#' @export
#'
#' @examples
#' data(DetrendedData)
#' model <- EstimateVAR(DetrendedData)
#' freq_model <- ParamFreqModel(model)
#' Phase_open <- GetOpenTransFunPhase(freq_model, 1, 2)
#' Phase_open[1:100]
GetOpenTransFunPhase <- function(SM, index1, index2){
  S <- GetMatrixAfromH(SM$Spectra)
  chosenS <- S[index1, index2,]
  phase <- atan2(Im(chosenS), Re(chosenS))
  return(phase)
}