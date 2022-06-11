#' @export
SpillComp <- function(data, output = NULL, threshold = 0.1, flexrep = 10, neighbor = 1, seed = 42) {
  spillmat_results <- GetSpillMat(data, threshold = threshold, flexrep = flexrep, neighbor = neighbor, seed = seed)
  spillmat <- spillmat_results[[1]]
  cutoffs <- spillmat_results[[2]]
  set.seed(seed)
  data_compensated <- t(apply(data, 1, function(row) nnls::nnls(t(spillmat), row)$x))
  colnames(data_compensated) <- colnames(data)
  if (!is.null(output)) {flowCore::write.FCS(flowCore::flowFrame(data_compensated), filename = output)}
  return(list(flowCore::flowFrame(data_compensated), spillmat, cutoffs))
}
