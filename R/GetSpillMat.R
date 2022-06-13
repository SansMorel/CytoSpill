
#' Get spillover matrix and cutoffs
#'
#' @param data A numerical matrix.
#' @param rows Rows in data to be used for calculating cutoffs.
#' @param threshold Spillover max boundary.
#' @param flexrep Passed to flexmix. Number of times flexmix is run. See flexmix nrep.
#' @param neighbor Abundance spillover channels. Can be 1 or 2.
#' @param quantile Quantile used for cutoff calculation.
#' @param seed Seed for reproducibility.
#' @param n_threads Number of threads used for calculation of cutoffs. Passed to .DeriveCutoffsHelper.
#'
#' @return List of spillover matrix and cutoffs.
#' @export
#'
GetSpillMat <- function(data, rows, threshold, flexrep, neighbor, quantile, seed, n_threads) {
  stopifnot("`data` must be a numeric matrix." = is.matrix(data) & is.numeric(data))
  set.seed(seed)
  cutoffs <- .DeriveCutoffsHelper(x = data[rows,], quantile = quantile, flexrep = flexrep, seed = seed, n_threads = n_threads)
  model <- .EstimateSpill(data, cutoffs, upperbound = threshold, neighbor = neighbor)
  estimates <- model[[1]]
  xcols <- model[[2]]
  spillmat <- diag(length(xcols))
  for (i in 1:length(xcols)) {
    if (!is.na(xcols[[i]][1])) {
      for (j in 1:length(xcols[[i]])) {
        if (!is.na(estimates[[i]][j])){
          spillmat[xcols[[i]][j],i] <- min(estimates[[i]][j],threshold)
        }
      }
    }
  }
  return(list(spillmat,cutoffs))
}
