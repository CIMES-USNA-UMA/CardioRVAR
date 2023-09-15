test_that("Getting expected values without coherence values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  ExpectedVals <- GetExpectedValues(freq_model, use.coh = FALSE, str = FALSE)
  expect_true(is.list(ExpectedVals))
})

test_that("Getting expected values with coherence values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coherence <- CalculateCoherence(freq_model, 1, 2)
  ExpectedVals <- GetExpectedValues(freq_model, coherence = coherence, str = FALSE)
  expect_true(is.list(ExpectedVals))
})

test_that("Getting expected values at maximum coherence should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coherence <- CalculateCoherence(freq_model, 1, 2)
  GetEstimateAtMaxCoh(freq_model, coherence = coherence)
  Estimates <- GetEstimateAtMaxCoh(freq_model, coherence = coherence,
                                   str = FALSE)
  expect_true(is.list(Estimates))
})

test_that("Getting peaks should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  Peaks <- GetPeaks(freq_model, str = FALSE)
  expect_true(is.list(Peaks))
})
 
