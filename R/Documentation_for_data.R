
#' @docType data
#' @name Cardiovascular
#' @title Raw Cardiovascular Data
#' @description
#' Data object containing systolic and diastolic blood pressure values, as well as RR intervals
#' and heart rate values.
#' @format A list with 5 variables:
#' \describe{
#'   \item{Time}{Time values}
#'   \item{RR}{RR intervals}
#'   \item{SBP}{Beat-to-beat systolic blood pressure}
#'   \item{DBP}{Beat-to-beat dyastolic blood pressure}
#'   \item{HR}{Heart rate values}
#' }
#' @usage data(Cardiovascular)
#' @keywords dataset
#' @examples
#' data(Cardiovascular)
#' 
"Cardiovascular"

#' @docType data
#' @name DetrendedData
#' @title Detrended Cardiovascular Data
#' @description
#' Detrended RR intervals and systolic blood pressure from \link[CardioRVAR]{Cardiovascular}
#' @format A matrix with 2 columns:
#' \describe{
#'   \item{SBP}{Detrended systolic blood pressure series}
#'   \item{RR}{Detrended RR intervals series}
#' }
#' @usage data(DetrendedData)
#' @keywords dataset
#' @examples
#' data(DetrendedData)
#' 
"DetrendedData"




