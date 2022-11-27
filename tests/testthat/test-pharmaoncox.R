testthat::test_that("Total records in database", {
  testthat::expect_gt(
    NROW(
      pharmaOncoX::get_drugs(cache_dir = "~/Downloads")$records), 45000)
})
