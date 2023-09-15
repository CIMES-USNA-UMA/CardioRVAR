test_that("Plotting closed-loop transfer function should work", {
   data(DetrendedData)
   model <- EstimateVAR(DetrendedData)
   freq_model <- ParamFreqModel(model)
   expect_true(is.null(PlotTransferFun(freq_model, 1, 2)))
})

test_that("Plotting open-loop transfer function should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  expect_true(is.null(PlotTransferFun(freq_model, 1, 2, open = TRUE)))
})

test_that("Plotting noise transfer function should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  expect_true(is.null(PlotNoiseTransferFun(freq_model, 1, 2)))
})

