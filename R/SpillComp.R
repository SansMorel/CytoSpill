#' Estimate spillover matrix and correct spillover
#' @param data A numerical matrix.
#' @param n Number of sampled events to be used for spillover estimation. n=0 uses all events.
#' @param cols Columns to be used for spillover estimation and spillover correction.
#' @param output Output FCS filename. Set output = NULL  to not write FCS file.
#' @param threshold Spillover max boundary. See https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.24298.
#' @param flexrep Passed to flexmix. Number of times flexmix is run. See flexmix nrep.
#' @param neighbor Abundance spillover channels. Can be 1 or 2.
#' @param seed seed for reproducibility.
#' @param n_threads Number of threads used for calculation of cutoffs.
#' @references
#'  Miao, Q, Wang, F, Dou, J, et al. Ab initio spillover compensation in mass cytometry data. Cytometry. 2021; 99: 899– 909. https://doi.org/10.1002/cyto.a.24298
#' @export
SpillComp <- function(data, n = 1e4, cols = NULL, output = NULL, threshold = 0.1, flexrep = 10, neighbor = 1, quantile = 0.1, seed = 42, n_threads = 8) {

  stopifnot("n must be non-negative" = n>=0)
  if(n > 0) {
    n <- min(nrow(data),n)
    rows <- sample(nrow(data), n)
  } else {
    rows <- 1:nrow(data)
  }
  if(is.null(cols)) cols <- 1:ncol(data)

  spillmat_results <- GetSpillMat(data[, cols], rows, threshold = threshold, flexrep = flexrep, neighbor = neighbor, quantile = quantile, seed = seed, n_threads = n_threads)
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
