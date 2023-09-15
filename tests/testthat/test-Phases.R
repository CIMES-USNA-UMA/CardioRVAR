test_that("Getting phase differences should work", {
   data(DetrendedData)
   model <- EstimateVAR(DetrendedData)
   freq_model <- ParamFreqModel(model)
   Phase <- GetTransFunPhase(freq_model, 1, 2)
   expect_true(is.numeric(Phase))
})


test_that("Getting phase differences for open loop should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model <- ParamFreqModel(model)
  Phase_open <- GetOpenTransFunPhase(freq_model, 1, 2)
  expect_true(is.numeric(Phase_open))
})