test_that("Estimation on stationary data should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  expect_true(is.list(model))
})

test_that("Estimation on non-stationary data should give a warning", {
  data(Cardiovascular)
  Data <- ResampleData(Cardiovascular)
  Use_data <- cbind(Data$SBP, Data$RR)
  expect_warning(EstimateVAR(Use_data))
})

test_that("Estimation on non-stationary data should not give an output", {
  data(Cardiovascular)
  Data <- ResampleData(Cardiovascular)
  Use_data <- cbind(Data$SBP, Data$RR)
  model <- EstimateVAR(Use_data, warnings = FALSE)
  expect_true(is.null(model))
})
