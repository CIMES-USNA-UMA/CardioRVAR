

# Get Phases


GetTransFunPhase <- function(SM, index1, index2){
  TF <- SM$Transfer_Functions[index1, index2,]
  phase <- atan2(Im(TF), Re(TF))
  return(phase)
}

GetOpenTransFunPhase <- function(SM, index1, index2){
  S <- GetMatrixAfromH(SM$Spectra)
  chosenS <- S[index1, index2,]
  phase <- atan2(Im(chosenS), Re(chosenS))
  return(phase)
}