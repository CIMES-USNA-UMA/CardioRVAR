test_that("Estimating coherence values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coherence <- CalculateCoherence(freq_model, 1, 2)
  expect_true(class(coherence) == "complex")
})

test_that("Plotting coherence values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coherence <- CalculateCoherence(freq_model, 1, 2)
  expect_true(is.null(PlotCoherence(freq_model, 1, 2, coherence)))
})

