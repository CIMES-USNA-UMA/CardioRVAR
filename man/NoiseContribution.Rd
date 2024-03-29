% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Noise_Contribution.R
\name{NoiseContribution}
\alias{NoiseContribution}
\title{Noise source contribution}
\usage{
NoiseContribution(
  SM,
  index1,
  index2,
  VLF = 0.04,
  LF = 0.15,
  HF = 0.4,
  use.coh = TRUE,
  thr = 0.5,
  coherence,
  print.flow = FALSE
)
}
\arguments{
\item{SM}{the computed frequency domain model from the system to be analyzed}

\item{index1}{numeric index to specify a variable from the model}

\item{index2}{numeric index to specify a variable from the model}

\item{VLF}{maximum limit for the VLF band. Default is 0.04}

\item{LF}{maximum limit for the LF band. Default is 0.15}

\item{HF}{maximum limit for the HF band. Default is 0.4}

\item{use.coh}{boolean, if TRUE the coherence and coherence threshold will be
used for the computation of expected values. Default is TRUE}

\item{thr}{a specific coherence threshold. Default is 0.5}

\item{coherence}{a vector with the spectral coherence values}

\item{print.flow}{Boolean. If TRUE, prints a message describing the analysed
contribution}
}
\value{
None
}
\description{
Estimates the noise source contribution for the variables of the analyzed system
}
\examples{
data(DetrendedData)
model <- EstimateVAR(DetrendedData)
freq_model <- ParamFreqModel(model)
noise_con <- NoiseContribution(freq_model, 1, 2, use.coh = FALSE)
noise_con

# The coherence can be used as a method for identifying reliable estimates
coherence <- CalculateCoherence(freq_model, 1, 2)
noise_con_thr <- NoiseContribution(freq_model, 1, 2, coherence = coherence)
noise_con_thr

}
\references{
Hytti H, Takalo R, Ihalainen H. Tutorial on Multivariate Autoregressive Modelling.
J Clin Monit Comput. 2006;20(1):101-8
}
\author{
Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
}
