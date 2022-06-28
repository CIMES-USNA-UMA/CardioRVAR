# Developed by Alvaro Chao-Ecija

# A series of functions to plot closed-loop diagrams, thanks to package
# diagram. It allows a better comprehension of the analyzed closed-loop.



PlotCLDiagram <- function(branch = 1){
  names <- c("Wsbp", "SBP", "RR", "Wrr")
  M <- matrix(nrow = 4, ncol = 4, byrow = TRUE, data = 0)
  M[2,1] <- "1/A[22]"
  M[2,3] <- "-A[12]/A[11]"
  M[3,2] <- "-A[21]/A[22]"
  M[3,4] <- "1/A[11]"
  col <- M
  col[] <- "black"
  if(branch == 1){
    col[2,3] <- "red"
    col2 <- c("lightgreen", "red", "lightblue", "lightgreen")
  }
  if(branch == 2){
    col[3,2] <- "red"
    col2 <- c("lightgreen", "lightblue", "red", "lightgreen")
  } 
  diagram::plotmat(M, pos = c(4), name = names, lwd = 1,
          box.lwd = 2, cex.txt = 2, box.cex = 2, box.size = 0.05,
          box.type = "circle", box.prop = 0.5, curve = 0.6, shadow.size = 0, 
          arr.col = col, box.col = col2, main = "Closed Loop Model", 
          cex.main = 2,
          arr.lcol = col)
}



PlotNoiseDiagram <- function(branch = 1, immediate = 1){
  M <- matrix(nrow = 4, ncol = 4, byrow = TRUE, data = 0)
  if(immediate == 1){
    names <- c("W RR", "W SBP", "RR", "SBP")
  } else {
    names <- c("W SBP", "W RR", "SBP", "RR")
  }
  M[3,1] <- "h[22]"
  M[4,1] <- "h[12]"
  M[3,2] <- "h[21]"
  M[4,2] <- "h[11]"
  col <- M
  col[] <- "black"
  if(branch == immediate){
    col[3,1] <- "red"
    col[3,2] <- "red"
    col2 <- c("lightgreen", "lightgreen", "red", "lightblue")
  }
  if(branch != immediate){
    col[4,1] <- "red"
    col[4,2] <- "red"
    col2 <- c("lightgreen", "lightgreen", "lightblue", "red")
  } 
  diagram::plotmat(M, pos = c(2,2), name = names, lwd = 1,
                   box.lwd = 2, cex.txt = 2, box.cex = 2, box.size = 0.05,
                   box.type = "circle", box.prop = 0.5, curve = 0, shadow.size = 0, 
                   arr.col = col, box.col = col2, main = "Noise Model",
                   cex.main = 2, arr.lcol = col, arr.pos = 0.4)
}


