test_that("Estimating causal coherence should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  ccoh <- CalculateCausalCoherence(freq_model, Mod = FALSE)
  expect_true(is.list(ccoh))
})

test_that("Estimating causal coherence's absolute values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  ccoh <- CalculateCausalCoherence(freq_model, Mod = TRUE)
  cond1 <- class(ccoh$C1) == "numeric"
  cond2 <- class(ccoh$C2) == "numeric"
  cond3 <- class(ccoh$Cr) == "numeric"
  expect_true(cond1 & cond2 & cond3)
})

test_that("Estimating causal phase should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  ccoh <- CalculateCausalCoherence(freq_model, Mod = FALSE)
  phase <- CalculateCausalPhase(ccoh)
  expect_true(is.list(phase))
})

test_that("Getting mean coherence should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coh <- CalculateCoherence(freq_model,1,2)
  means <- GetMeanCoherence(freq_model, coh, weight = FALSE)
  expect_true(is.list(means))
})


test_that("Getting weighted mean coherence should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coh <- CalculateCoherence(freq_model,1,2)
  means <- GetMeanCoherence(freq_model, coh, weight = TRUE)
  expect_true(is.list(means))
})

test_that("Getting maximum coherence should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coh <- abs(CalculateCoherence(freq_model, 1, 2))^2
  max_coh <- GetMaxCoherence(freq_model, coh)
  expect_true(is.list(max_coh))
})

 



