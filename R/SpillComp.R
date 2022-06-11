#' @export
SpillComp <- function(data, n = 1e4, cols = NULL, output = NULL, threshold = 0.1, flexrep = 10, neighbor = 1, seed = 42) {

  stopifnot("n must be non-negative" = n>=0)
  if(n > 0) {
    n <- min(nrow(data),n)
    rows <- sample(nrow(data), n)
  } else {
    rows <- 1:nrow(data)
  }
  if(is.null(cols)) cols <- 1:ncol(data)


  spillmat_results <- GetSpillMat(data[, cols], rows, threshold = threshold, flexrep = flexrep, neighbor = neighbor, seed = seed)
  spillmat <- spillmat_results[[1]]
  cutoffs <- spillmat_results[[2]]
  set.seed(seed)
  data[,cols] <- t(apply(data[,cols], 1, function(row) nnls::nnls(t(spillmat), row)$x))
  if (!is.null(output)) {flowCore::write.FCS(flowCore::flowFrame(data), filename = output)}
  return(list(flowCore::flowFrame(data), spillmat, cutoffs))
}
