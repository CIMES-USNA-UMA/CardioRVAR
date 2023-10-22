test_that("Stationarity check should work", {
  data(DetrendedData)
  expect_true(CheckStationarity(DetrendedData, warnings = F))
})

test_that("Stationarity check should work and should give a statement if required", {
  data(DetrendedData)
  expect_true(CheckStationarity(DetrendedData, warnings = F))
  First <- CheckStationarity(DetrendedData, warnings = F)
  Second <- CheckStationarity(DetrendedData, verbose = TRUE, verbose.method = "p",
                              warnings = F)
  expect_true(First & is.character(Second))
})

test_that("Stationarity check should recognize non-stationarity", {
  data(Cardiovascular)
  Data <- ResampleData(Cardiovascular)
  Use_data <- cbind(Data$SBP, Data$RR)
  expect_true(!CheckStationarity(Use_data, warnings = FALSE))
})

test_that("Stationarity check should give a warning when data is non-stationary", {
  data(Cardiovascular)
  Data <- ResampleData(Cardiovascular)
  Use_data <- cbind(Data$RR)
  expect_warning(CheckStationarity(Use_data))
})


test_that("Stability check should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  expect_true(DiagnoseStability(model))
})

test_that("Residuals check should work", {
  data(DetrendedData)
  model <- EstimateVAR(DetrendedData)
  expect_true(DiagnoseResiduals(model))
})



