#' @export
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
