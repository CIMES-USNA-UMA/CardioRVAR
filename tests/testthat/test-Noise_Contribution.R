test_that("Estimating noise source contribution without coherence values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  noise_con <- NoiseContribution(freq_model, 1, 2, use.coh = FALSE)
  expect_true(is.vector(noise_con))
})

test_that("Estimating noise source contribution with coherence values should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  coherence <- CalculateCoherence(freq_model, 1, 2)
  noise_con_thr <- NoiseContribution(freq_model, 1, 2, coherence = coherence)
  expect_true(is.vector(noise_con_thr))
})

test_that("Total noise source contribution should be 100%", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  noise_con1 <- NoiseContribution(freq_model, 1, 2, use.coh = FALSE)
  noise_con2 <- NoiseContribution(freq_model, 1, 1, use.coh = FALSE)
  total <- noise_con1 + noise_con2
  cond1 <- all.equal(total[1], c(HF = 100))
  cond2 <- all.equal(total[2], c(LF = 100))
  expect_true(cond1 & cond2)
})

test_that("Plotting noise source contribution should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  noise_con1 <- NoiseContribution(freq_model, 1, 2, use.coh = FALSE)
  noise_con2 <- NoiseContribution(freq_model, 1, 1, use.coh = FALSE)
  expect_true(is.null(PlotNoiseContribution(noise_con1, noise_con2)))
})

test_that("Plotting causality should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  expect_true(is.null(PlotCausality(freq_model, 1)))
})

