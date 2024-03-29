% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Detrend.R
\name{DetrendWithCutoff}
\alias{DetrendWithCutoff}
\title{Detrend variables}
\usage{
DetrendWithCutoff(
  x,
  thr = NULL,
  per = FALSE,
  cutoff = 0.04,
  f = 4,
  wv = "d16",
  max_f = 0.4,
  do.VHF = FALSE,
  use.universal = FALSE,
  n.universal = NULL,
  thr.type = c("hard", "soft")
)
}
\arguments{
\item{x}{A vector of observations}

\item{thr}{Manual threshold. Default is NULL.}

\item{per}{Boolean, use a percentage as a threshold value}

\item{cutoff}{A frequency used as limit for the detrending, or a reference
frequency for calculating said limit. Default is 0.04 Hz}

\item{f}{Sample rate of the vector of observations. Default is 4 Hz}

\item{wv}{Wavelet to be passed into the function \link[waveslim]{modwt} of package
\href{https://CRAN.R-project.org/package=waveslim}{waveslim}. Default is
d16}

\item{max_f}{Maximum frequency to be accounted for in the detrending. default is 0.4 Hz}

\item{do.VHF}{Boolean. If true, Very High Frequency bands will also be processed.}

\item{use.universal}{Boolean. If TRUE, universal thresholding will be performed.}

\item{n.universal}{Number of scales in which to applied the universal thresholding, if
applied. If NULL (default), it will be computed automatically.}

\item{thr.type}{Perform hard or soft thresholdings}
}
\value{
The detrended vector of observations
}
\description{
Detrends a vector of observations using the discrete wavelet transform
}
\details{
This function uses a wavelet detrending algorithm to select the aproppriate componets of a signal
and discard low frequency components. This function requires package \href{https://CRAN.R-project.org/package=waveslim}{waveslim}:
for more information regarding this package, please check the references section.
}
\examples{
data(Cardiovascular)
int_data <- ResampleData(Cardiovascular)
RR <- DetrendWithCutoff(int_data$RR) # The RR series is detrended
SBP <- DetrendWithCutoff(int_data$SBP) # The SBP series is detrended
}
\references{
Brandon Whitcher (2020). waveslim: Basic Wavelet Routines for One-, Two-, and Three-Dimensional Signal
Processing. R package version 1.8.2. https://CRAN.R-project.org/package=waveslim

Li L, Liu C, Li K, Liu C. Comparison of Detrending Methods in Spectral Analysis of Heart Rate Variability.
Res J Appl Sci Eng Technol. 2011;3(9):1014-21
}
\author{
Alvaro Chao-Ecija, Marc Stefan Dawid-Milner
}
