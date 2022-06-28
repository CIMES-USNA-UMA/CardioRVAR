





PreprocessData2 <- function(data){
  names <- names(data)
  pos <- match(c("RR", "SBP"), names)
  for(n in pos){
    data[[n]] <- data[[n]] - mean(data[[n]])
    model <- summary(lm(data[[n]] ~ data$Time))$coefficients
    l_trend <- model[2,1] * data$Time + model[1,1]
    data[[n]] <- data[[n]] - l_trend
    data[[n]] <- dplR::pass.filt(data[[n]], 0.02, "high", n = 6)
    data[[n]] <- dplR::pass.filt(data[[n]], 0.4, "low", n = 6)
    data[[n]] <- data[[n]] - mean(data[[n]])
  }
  return(data)
}