% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Causal_Coherence.R
\name{GetMeanCoherence}
\alias{GetMeanCoherence}
\title{Get mean coherence}
\usage{
GetMeanCoherence(
  SM,
  coherence,
  HF = 0.4,
  LF = 0.15,
  VLF = 0.04,
  weight = FALSE
)
}
\arguments{
\item{SM}{the computed frequency domain model from the system to be analyzed}

\item{coherence}{the coherence function to be decomposed}

\item{HF}{maximum limit for the HF band. Default is 0.4}

\item{LF}{maximum limit for the LF band. Default is 0.15}

\item{VLF}{maximum limit for the VLF band. Default is 0.04}

\item{weight}{applies gaussian-weighting function to the estimates before
calculating the expected values. Default is FALSE}
}
\value{
A list with the computed expected values at the LF and HF bands for the 
        coherence
}
\description{
Computes the expected coherence values for the HF and LF bands
}
\examples{
data(DetrendedData)
model <- EstimateVAR(DetrendedData)
freq_model <- ParamFreqModel(model)
coh <- CalculateCoherence(freq_model, 1, 2)
GetMeanCoherence(freq_model, coh, weight = FALSE) # No Gaussian weights
GetMeanCoherence(freq_model, coh, weight = TRUE) # Gaussian weights are applied
}
\author{
Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
}
