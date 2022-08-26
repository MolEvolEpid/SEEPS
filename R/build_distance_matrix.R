#' Reduce a geneology to a pairwise distance matrix.
#'
#' @export
geneology_to_distance_matrix <- function(geneology, spike_root = FALSE) {
  # Convert a geneology to a pairwise difference matrix
  # spike_root option adds a row and column for a node at the root

  M = length(geneology[, 1]) / 2

  # if (spike_root) {
  #   M <- M + 1
  # }  # 1 more root node for multicluster data
  mrca_mat <- matrix(0, M, M)
  pair_diff_mat <- matrix(0, M, M, )
  ind_par_mat <- geneology
  # if (spike_root) { M <- M - 1 }
  ## Obtain an MRCA matrix
  for (k in 1:(M - 1)) {
    ancestors = rep(1e+5, M)
    pp = k
    ii = 1
    while (pp != 0) {
      ancestors[ii] = pp
      ii = ii + 1
      pp = geneology[pp, 2]
    }
    for (l in (k + 1):M) {
      pp = l

      while (!(pp %in% ancestors)) {
        pp = geneology[pp, 2]
      }

      mrca_mat[k, l] <- pp
      mrca_mat[l, k] <- pp
    }
  }

  for (k in 1:(M - 1)) {
    for (l in (k + 1):M) {
      finally <- mrca_mat[k, l]

      holder <- ind_par_mat[k, 2]
      sumhold <- ind_par_mat[k, 5]
      while (holder != finally) {
        sumhold <- sumhold + ind_par_mat[holder, 5]
        holder <- ind_par_mat[holder, 2]

      }

      holder <- ind_par_mat[l, 2]
      sumhold <- sumhold + ind_par_mat[l, 5]
      while (holder != finally) {
        sumhold <- sumhold + ind_par_mat[holder, 5]
        holder <- ind_par_mat[holder, 2]

      }
      pair_diff_mat[k, l] <- sumhold
      pair_diff_mat[l, k] <- sumhold
    }
  }
  return(pair_diff_mat)
}


###########################################################
#
# For BPB - need to read nodes off.
# Todo: merge into above function.
#
###########################################################

#' Reduce a geneology to a pairwise distance matrix.
#'
#' @export
geneology_to_distance_matrix_bpb <- function(geneology, spike_root = FALSE) {
  # Convert a geneology to a pairwise difference matrix
  # spike_root option adds a row and column for a node at the root
  # print(geneology)
  M <- sum(geneology[,6])  # number of leaf nodes

  inds <- which(geneology[,6] != 0) # Get the index of the nonzeros
  # Todo: Optimize this.
  # We need to convert global indexing (1:# of tracked infections) to local indexing (1:M)
  global_to_local <- rep(0, max(inds))
  global_to_local[inds] <- seq_along(M)  # So g_2_l[global] = local for O(1) lookup
  # if (spike_root) {
  #   M <- M + 1
  # }  # 1 more root node for multicluster data
  mrca_mat <- matrix(0, M, M)
  pair_diff_mat <- matrix(0, M, M, )
  ind_par_mat <- geneology
  # if (spike_root) { M <- M - 1 }
  ## Obtain an MRCA matrix
  for (k in inds) {
    ancestors = rep(-1, M)
    pp = k
    ii = 1
    while (pp != 0) {
      ancestors[ii] = pp
      ii = ii + 1
      pp = geneology[pp, 2]
    }
    for (l in (k + 1):M) {
      pp = l

      while (!(pp %in% ancestors)) {
        pp = geneology[pp, 2]
      }
      # print(c(k, l))
      # print(ancestors)
      mrca_mat[global_to_local[k], global_to_local[l]] <- pp
      mrca_mat[global_to_local[l], global_to_local[k]] <- pp
    }
  }

  for (k in 1:(M - 1)) {
    for (l in (k + 1):M) {
      finally <- mrca_mat[k, l]

      holder <- ind_par_mat[k, 2]
      sumhold <- ind_par_mat[k, 5]
      while (holder != finally) {
        sumhold <- sumhold + ind_par_mat[holder, 5]
        holder <- ind_par_mat[holder, 2]

      }

      holder <- ind_par_mat[l, 2]
      sumhold <- sumhold + ind_par_mat[l, 5]
      while (holder != finally) {
        sumhold <- sumhold + ind_par_mat[holder, 5]
        holder <- ind_par_mat[holder, 2]

      }
      pair_diff_mat[k, l] <- sumhold
      pair_diff_mat[l, k] <- sumhold
    }
  }
  return(pair_diff_mat)
}
