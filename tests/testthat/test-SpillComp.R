
load(file=file.path(system.file("data", package="CytoSpill"), "Levine32_example.Rdata"))
results <- SpillComp(data = data_Levine32[1:1000,5:36], n = 0, cols = NULL, threshold = 0.1, flexrep = 5, neighbor = 1)
results_precomputed <- readRDS("Levine32_example_output.rds")

test_that("Compensation works", {
  expect_equal(flowCore::exprs(results[[1]]), flowCore::exprs(results_precomputed[[1]]))
  expect_equal(results[[2]], results_precomputed[[2]])
  expect_equal(results[[3]], results_precomputed[[3]])
})
