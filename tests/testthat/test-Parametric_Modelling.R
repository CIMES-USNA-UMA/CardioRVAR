test_that("Modelling should work", {
   data(DetrendedData)
   model <- EstimateVAR(DetrendedData)
   freq_model <- ParamFreqModel(model)
   expect_true(is.list(freq_model))
})

test_that("Non inncluding a0 effects should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model_a0 <- ParamFreqModel(model)
  freq_model_no_a0 <- ParamFreqModel(model, A0 = FALSE)
  expect_true(sum(all.equal(freq_model_a0$a0, freq_model_no_a0$a0) == TRUE) == 0)
})


test_that("Including a0 effects manually should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  freq_model_a0 <- ParamFreqModel(model)
  freq_model_no_a0 <- ParamFreqModel(model, A0 = FALSE)
  freq_model_pos_a0 <- IncludeA0Effects(freq_model_no_a0)
  expect_true(all.equal(freq_model_a0$a0, freq_model_pos_a0$a0))
})


