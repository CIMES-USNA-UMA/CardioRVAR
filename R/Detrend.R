# Developed by Alvaro Chao-Ecija
#
# Function for wavelet-detrending the cardiovascular signals. 


DetrendByCutoff <- function(x, cutoff = 0.04, f = 4, wv = "d16", 
                            max_f = 0.4){
  N = NROW(x)
  Nf = f/2
  level = 1
  repeat{
    target_f <- Nf / (2^level)
    if(target_f <= cutoff){
      break
    } else {
      level = level + 1
    }
  }
  level2 = 1
  repeat{
    target_m <- Nf / (2^level2)
    if(target_m <= max_f){
      break
    } else {
      level2 = level2 + 1
    }
  }
  dx = waveslim::modwt(x - mean(x), wf = wv, level)
  #for(n in 3:5) dx[[n]] <- double(N)
  dx <- waveslim::universal.thresh.modwt(dx, level, hard = TRUE)
  for(n in level2:(level-1)) dx[[n]] <- double(N)
  trend <- imodwt(dx)
  return(x - mean(x) - trend)
}


