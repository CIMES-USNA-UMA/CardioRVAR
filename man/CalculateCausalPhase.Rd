% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Causal_Coherence.R
\name{CalculateCausalPhase}
\alias{CalculateCausalPhase}
\title{Calculate phase difference from causal coherence}
\usage{
CalculateCausalPhase(ccoh)
}
\arguments{
\item{ccoh}{A vector containing causal coherence values}
}
\value{
The computed phase difference
}
\description{
Calculate phase difference from causal coherence
}
\examples{
data(DetrendedData)
model <- EstimateVAR(DetrendedData)
freq_model <- ParamFreqModel(model)
ccoh <- CalculateCausalCoherence(freq_model, Mod = FALSE)
phase <- CalculateCausalPhase(ccoh)
}
\author{
Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
}
