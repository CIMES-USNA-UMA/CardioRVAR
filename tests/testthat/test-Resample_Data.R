test_that("Resampling should work", {
   data(Cardiovascular)
   interpolated <- ResampleData(Cardiovascular)
   expect_true(is.list(interpolated))
})
