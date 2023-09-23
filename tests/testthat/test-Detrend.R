test_that("Detrending should work", {
   data(Cardiovascular)
   data(DetrendedData)
   int_data <- ResampleData(Cardiovascular)
   RR <- DetrendWithCutoff(int_data$RR) 
   SBP <- DetrendWithCutoff(int_data$SBP) 
   expect_true(all.equal(RR, DetrendedData[,"RR"]) & all.equal(SBP, DetrendedData[,"SBP"]))
})
