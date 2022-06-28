


GenerateReport <- function(var = NULL, index1 = 2, index2 = 1, freq_model = NULL, include.plot = TRUE,
          use.var = TRUE, name = NULL){
          folder <- tcltk::tk_choose.dir()
          setwd(folder)
          doc <- officer::read_docx()
          doc <- officer::body_add_par(doc, "PhysioMVAR Report",
              style = "heading 1")
          doc <- officer::body_add_par(doc, paste("", sep = ""))
          par1 <- paste("This report was generated on date:", date())
          doc <- officer::body_add_par(doc, par1, style = "Normal")
          if(!is.null(var)){
             names <- colnames(var$y)
             coefs <- GetCoefs(var)
             sigma <- summary(var)$cov
             GetA0 <- GetA0Fun(sigma)
             nsigma <- GetA0$sigma
             d <- dim(coefs[[1]])[1]
             A0 <- GetA0$a0
             I <- diag(1, d)
             ncoefs <- UpdateWithA0(A0, coefs)
             models <- models2 <- list()
             length(models) <- length(models2) <- d
             for(n in 1:d){
                 model <- matrix(0, nrow = d, ncol = length(coefs) + 1)
                 model2 <- matrix(0, nrow = d, ncol = length(coefs) + 1)
                 model[,1] <- I[n,]
                 model2[,1] <- A0[n,] 
                 for(m in 1:length(coefs)){
                     model[,m+1] <- coefs[[m]][n,] 
                     model2[,m+1] <- ncoefs[[m]][n,] 
                 }
                 rownames(model) <- rownames(model2) <- names
                 colnames(model) <- colnames(model2) <- paste("Lag", 0:length(coefs))
                 models[[n]] <- cbind(Variable = names, round(model, 3))
                 models2[[n]] <- cbind(Variable = names, round(model2, 3))
             }
             sigma <- sigma 
             nsigma <- nsigma 
             A0 <- A0 
             stability <- DiagnoseStability(var)
             white <- DiagnoseResiduals(var)
             doc <- officer::body_add_par(doc, paste("", sep = ""))
             doc <- officer::body_add_par(doc, "VAR Model", style = "heading 1")
             doc <- officer::body_add_par(doc, paste("", sep = ""))
             parVAR1 <- "A VAR model has been submitted. Its characteristics are shown below:"
             doc <- officer::body_add_par(doc, parVAR1, style = "Normal")
             if(stability & white){
                valid <- "Yes"
             } else {
                valid <- "No"
             }
             if(stability){
                stability <- "Yes"
             } else {
                stability <- "No"
             }
             if(white){
                white <- "Yes"
             } else {
                white <- "No"
             }
             characteristics <- data.frame(Order = var$p, Criterion = names(var$p)[1], Stable = stability,
                White = white, Valid  =valid)
             colnames(characteristics) <- c("Order", "Criterion", "Stable", "White Noise Residuals",
               "Valid")
             sigma <- data.frame(cbind(Variables = names, round(sigma,3)))
             nsigma <- data.frame(cbind(Variables = names, round(nsigma,3)))
             A0 <- data.frame(cbind(Variables = names, round(A0,3)))
             colnames(sigma) <- colnames(nsigma) <- colnames(A0) <- c("Variable", names)
             doc <- officer::body_add_par(doc, paste("", sep = ""))
             doc <- officer::body_add_par(doc, "Model Characteristics", style = "table title")
             doc <- officer::body_add_table(doc, characteristics, style = "table_template")
             if(valid == "No"){
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                parVAR2 <- "The model is not valid. No further results will be estimated"
                doc <- officer::body_add_par(doc, parVAR2, style = "Normal")
             } else {
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                parVAR2 <- paste("The model has in total", d, "branches. Bellow, each model branch is shown:")
                doc <- officer::body_add_par(doc, parVAR2, style = "Normal")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                for(n in 1:length(models)){
                    doc <- officer::body_add_par(doc, paste("Branch ", n, " (", names[n], ")", sep = "" ), style = "table title")
                    tab <- data.frame(models[[n]])
                    colnames(tab) <- c("Variable", paste("Lag", 0:length(coefs)))
                    doc <- officer::body_add_table(doc, tab, style = "table_template")
                    doc <- officer::body_add_par(doc, paste("", sep = ""))
                }
                parVAR3 <- "Bellow, the noise covariance matrix is shown:"
                doc <- officer::body_add_par(doc, parVAR3, style = "Normal")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                doc <- officer::body_add_par(doc, "Noise Covariance Matrix", style = "table title")
                doc <- officer::body_add_table(doc, sigma, style = "table_template")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                doc <- officer::body_add_par(doc, "PhysioMVAR Report",
                  style = "heading 2")
                parVAR4 <- "The following no-delay effects have been detected:"
                doc <- officer::body_add_par(doc, parVAR4, style = "Normal")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                doc <- officer::body_add_par(doc, "No-Delay Effects", style = "table title")
                doc <- officer::body_add_table(doc, A0, style = "table_template")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                parVAR5 <- "The adjusted model for no-delay effects is shown:"
                doc <- officer::body_add_par(doc, parVAR5, style = "Normal")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                for(n in 1:length(models2)){
                    doc <- officer::body_add_par(doc, paste("Branch ", n, " (", names[n], ")", sep = "" ), style = "table title")
                    tab <- data.frame(models2[[n]])
                    colnames(tab) <- c("Variable", paste("Lag", 0:length(ncoefs)))
                    doc <- officer::body_add_table(doc, tab, style = "table_template")
                    doc <- officer::body_add_par(doc, paste("", sep = ""))
                }
                doc <- officer::body_add_par(doc, paste("", sep = ""))
                doc <- officer::body_add_par(doc, "Noise Covariance Matrix", style = "table title")
                doc <- officer::body_add_table(doc, nsigma, style = "table_template")
                doc <- officer::body_add_par(doc, paste("", sep = ""))
             }  
           }
           if(!is.null(var) & use.var & (valid == "Yes")){
              freq_model <- PhysioMVAR::ParamFreqModel(var)
           }
           n_names <- names[c(index1, index2)]
           doc <- officer::body_add_par(doc, "System Characteristics",
              style = "heading 1")
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           system_char <- data.frame(Input = n_names[1], Output = n_names[2])
           doc <- officer::body_add_par(doc, "Analyzed System", style = "table title")
           doc <- officer::body_add_table(doc, system_char, style = "table_template")
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, "Closed Loop Feedforward Transfer Function",
              style = "heading 1")
           evals <- GetExpectedValues(freq_model, str = FALSE)
           c_evals <- rbind(c(evals$HF$Transfer_Functions[index2,index1,1],
               evals$HF$Transfer_Functions[index1,index2,1]), c(evals$LF$Transfer_Functions[index2,index1,1],
                evals$LF$Transfer_Functions[index1,index2,1]))
           o_evals <- rbind(c(evals$HF$Open_Transfer_Functions[index2,index1,1],
               evals$HF$Open_Transfer_Functions[index1,index2,1]), c(evals$LF$Open_Transfer_Functions[index2,index1,1],
                evals$LF$Open_Transfer_Functions[index1,index2,1]))
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           tf <- PlotTransferFun(freq_model, index2, index1, tem = TRUE)
           doc <- officer::body_add_img(doc, tf, style = "centered", width = 6, height = 6)
           unlink(tf)
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, paste("", sep = ""))
        doc <- officer::body_add_table(doc, data.frame(LF = round(c_evals[2,1],3), 
           HF = round(c_evals[1,1],3)), style = "table_template")
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, "Closed Loop Feedback Transfer Function",
              style = "heading 1")
           tf <- PlotTransferFun(freq_model, index1, index2, tem = TRUE)
           doc <- officer::body_add_img(doc, tf, style = "centered", width = 6, height = 6)
           unlink(tf)
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, paste("", sep = ""))
        doc <- officer::body_add_table(doc, data.frame(LF = round(c_evals[2,2],3), 
           HF = round(c_evals[1,2],3)), style = "table_template")
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, "Open Loop Feedforward Transfer Function",
              style = "heading 1")
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           otf <- PlotTransferFun(freq_model, index2, index1, open = TRUE, tem = TRUE)
           doc <- officer::body_add_img(doc, otf, style = "centered", width = 6, height = 6)
           unlink(otf)
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, paste("", sep = ""))
         doc <- officer::body_add_table(doc, data.frame(LF = round(o_evals[2,1],1), 
           HF = round(o_evals[1,1],3)), style = "table_template")
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, "Open Loop Feedback Transfer Function",
              style = "heading 1")
           otf <- PlotTransferFun(freq_model, index1, index2, open = TRUE, tem = TRUE)
           doc <- officer::body_add_img(doc, otf, style = "centered", width = 6, height = 6)
           unlink(otf)
           doc <- officer::body_add_par(doc, paste("", sep = ""))
           doc <- officer::body_add_par(doc, paste("", sep = ""))
        doc <- officer::body_add_table(doc, data.frame(LF = round(o_evals[2,2],3), 
          HF = round(o_evals[1,2],3)), style = "table_template")
          if(is.null(name)) name <- "PhysioMVAR Report"
          print(doc, target = paste(name, ".docx", sep = ""))
}
