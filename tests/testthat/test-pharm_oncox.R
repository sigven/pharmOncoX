testthat::test_that("Total records in database", {
  testthat::expect_gt(
    NROW(
      pharmOncoX::get_drugs(cache_dir = "~/Downloads")$records), 45000)
})
