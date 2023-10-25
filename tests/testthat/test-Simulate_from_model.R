test_that("Simulating model should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  new_model <- SimulateWithModel(freq_model, c(2,3), a0 = 2)
  expect_true(is.array(new_model))
})

test_that("Computing PSD should work", {
  data(DetrendedData)
  RR <- DetrendedData[,"RR"]
  PSD(RR, 21)
  HRV <- PSD(RR, 21, plot = FALSE, output = TRUE)
  expect_true(is.vector(HRV))
})

