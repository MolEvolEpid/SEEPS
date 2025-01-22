#' Reduce a large matrix by randomly sampling a cluster of closely related individuals
#'
#' An obtained sample (through contact tracing or random sampling) may be larger than needed.
#' This function extracts a subset of closely related individuals with an
#' optionally selected "center". This respects the `spike_root` option, and if
#' `spike_root = TRUE`, will return the root individual in the last column.
#'
#' @param oversampled_matrix A matrix to subsample rows/columns from.
#' @param subsample_size The number of individuals to keep in the reduced matrix.
#' @param spike_root If TRUE, include the root individual in the reduced matrix.
#' @param root_position The index of the root individual in the matrix.
#' @param index_id if not -1, use this index as the center of the cluster. If -1,
#' randomly sample a new center.
#' @param sort_order If NULL (default), sort the resulting matrix by the distance
#' from the center. If not null, sorts using this as an index vector. Selects
#' the first `subsample_size` by position.
#' @export
reduce_large_matrix <- function(oversampled_matrix, subsample_size,
                                spike_root = FALSE, root_position = 0,
                                index_id = -1, sort_order = NULL) {
    # N_expected - how many sequences to consider. Do not include the root infection
    # oversampled_matrix - matrix of full population

    # Get dimensions of problem
    data_length <- dim(oversampled_matrix)[2] # how many sequences in the dataset, [# of nrows, ncols]
    # For more elaborate subsampling of the matrix, construct a boolean mask
    mask <- rep(TRUE, data_length)

    # Don't sample the root as a possible center for the cluster.
    if (spike_root) mask[root_position] <- FALSE
    # Get the sample
    if (index_id == -1) {
        sample <- sample(x = data_length, size = 1, replace = FALSE, prob = mask / sum(mask))
    } else {
        # If we are going to reduce multiple matrices from the same sample,
        # We want to use the same center/index case for each matrix.
        # This assumes the initial draw is done correctly (above logic)
        sample <- index_id
    }
    if (is.null(sort_order)) {
        new_order <- sort(oversampled_matrix[sample, ], index.return = TRUE)
        ix <- new_order$ix
    } else {
        # Expect a vector if indices to sort by
        ix <- sort_order
    }
    closest_k <- ix[1:subsample_size]

    # If we need the root, place it at the end
    if (spike_root) closest_k <- as.vector(append(closest_k, data_length))
    # Now subset the matrix
    data_matrix <- oversampled_matrix[closest_k, closest_k]
    return(list(
        "matrix" = data_matrix, "keep_indices" = closest_k,
        "sampled_index" = sample, "sort_indices" = ix
    ))
}
