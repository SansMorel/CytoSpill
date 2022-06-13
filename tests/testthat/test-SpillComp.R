
load(file=file.path(system.file("data", package="CytoSpill"), "Levine32_example.Rdata"))
results <- SpillComp(data = data_Levine32[1:1000,5:36], n = 0, cols = NULL, threshold = 0.1, flexrep = 5, neighbor = 1)
results_precomputed <- readRDS("Levine32_example_output.rds")

test_that("Compensation works", {
  expect_equal(flowCore::exprs(results$compensated_fcs), flowCore::exprs(results_precomputed$compensated_fcs))
  expect_equal(results$spillover_matrix, results_precomputed$spillover_matrix)
  expect_equal(results$cutoffs, results_precomputed$cutoffs)
})
