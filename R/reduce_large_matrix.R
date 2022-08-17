#' Reduce a large matrix by randomly sampling a cluster of closely related individuals
#'
#' An obtained sample (through contact tracing or random sampling) may be larger than needed.
#' This function extracts a subset of closely related individuals with a randomly selected "center".
#' This respects the `spike_root` option, and if `spike_root = TRUE`, will return the root individual
#' in the last column.
#'
#' @export
reduce_large_matrix <- function(oversampled_matrix, subsample_size,
                                spike_root = FALSE, root_position = 0) {
  # N_expected - how many sequences to consider. Do not include the root infection
  # oversampled_matrix - matrix of full population

  # Get dimensions of problem
  data_length <- dim(oversampled_matrix)[2]  # how many sequences in the dataset, [# of nrows, ncols]
  # For more elaborate subsampling of the matrix, construct a boolean mask
  mask <- rep(TRUE, data_length)

  # Don't sample the root as a possible center for the cluster.
  if (spike_root) mask[root_position] <- FALSE
  # Get the sample
  sample <- sample(x = data_length, size = 1, replace = FALSE, prob = mask/sum(mask))
  new_order <- sort(oversampled_matrix[sample, ], index.return = TRUE)
  closest_k <- new_order$ix[1:subsample_size]

  # If we need the root, place it at the end
  if (spike_root) closest_k <- as.vector(append(closest_k, data_length))
  # Now subset the matrix
  data_matrix <- oversampled_matrix[closest_k, closest_k]
  return(data_matrix)

}
