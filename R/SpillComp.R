#' @export
SpillComp <- function(data, n = 1e4, cols = NULL, output = NULL, threshold = 0.1, flexrep = 10, neighbor = 1, seed = 42, n_threads = 8) {

  stopifnot("n must be non-negative" = n>=0)
  if(n > 0) {
    n <- min(nrow(data),n)
    rows <- sample(nrow(data), n)
  } else {
    rows <- 1:nrow(data)
  }
  if(is.null(cols)) cols <- 1:ncol(data)

  spillmat_results <- GetSpillMat(data[, cols], rows, threshold = threshold, flexrep = flexrep, neighbor = neighbor, seed = seed, n_threads = n_threads)
  names(spillmat_results) <- c("spillover_matrix", "cutoffs")
  rownames(spillmat_results$spillover_matrix) <- colnames(spillmat_results$spillover_matrix) <- colnames(data[, cols])
  set.seed(seed)
  data[,cols] <- t(apply(data[,cols], 1, function(row) nnls::nnls(t(spillmat_results$spillover_matrix), row)$x))
  spillmat_results$compensated_fcs <- flowCore::flowFrame(data)
  spillmat_results <- list(compensated_fcs = spillmat_results$compensated_fcs,
                           spillover_matrix = spillmat_results$spillover_matrix,
                           cutoffs = spillmat_results$cutoffs)

  if (!is.null(output)) {flowCore::write.FCS(spillmat_results$compensated_fcs, filename = output)}
  return(spillmat_results)
}
