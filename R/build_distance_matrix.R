#' Reduce a geneology to a pairwise distance matrix.
#'
#' Provided for backwards compatability with [Kupperman et al 2022]. New code
#' should use `[geneology_to_distance_matrix]`.
#'
#' @seealso geneology_to_distance_matrix
#' @export
geneology_to_distance_matrix_classic <- function(geneology, spike_root = FALSE) {  # nolint: cyclocomp_linter
  # Convert a geneology to a pairwise difference matrix
  # spike_root option adds a row and column for a node at the root

  M <- length(geneology[, 1]) / 2  # nolint: object_name_linter

  mrca_mat <- matrix(0, M, M)
  pair_diff_mat <- matrix(0, M, M)
  ind_par_mat <- geneology
  # Obtain an MRCA matrix
  for (k in 1:(M - 1)) {
    ancestors <- rep(1e+5, M)
    pp <- k
    ii <- 1
    while (pp != 0) {
      ancestors[ii] <- pp
      ii <- ii + 1
      pp <- geneology[pp, 2]
    }
    for (l in (k + 1):M) {
      pp <- l

      while (!(pp %in% ancestors)) {
        pp <- geneology[pp, 2]
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
geneology_to_distance_matrix <- function(geneology, mode="mu", spike_root = FALSE) {  # nolint: object_length_linter
  # Convert a geneology to a pairwise difference matrix
  # spike_root option adds a row and column for a node at the root
  if (spike_root) stop("This option is not yet implemented. See classic version.")
  # number of leaf nodes
  M <- sum(geneology[, 6])  # nolint: object_name_linter

  inds <- which(geneology[, 6] != 0) # Get the index of the nonzeros
  index_names <- geneology[inds, 7] # Get the names of the nonzeros
  # Todo: Optimize this.
  # We need to convert global indexing (1:# of tracked infections) to local indexing (1:M)
  global_to_local <- rep(0, max(inds))
  global_to_local[inds] <- seq(M)  # So g_2_l[global] = local for O(1) lookup

  mrca_mat <- matrix(0, M, M)
  pair_diff_mat <- matrix(0, M, M)
  ind_par_mat <- geneology

  rownames(pair_diff_mat) <- index_names
  colnames(pair_diff_mat) <- index_names
  ## Obtain an MRCA matrix
  for (k in inds) {
    ancestors <- rep(-1, M)
    pp <- k
    ii <- 1
    while (pp != 0) {
      ancestors[ii] <- pp
      ii <- ii + 1
      pp <- geneology[pp, 2]
    }
    for (l in (k + 1):M) {
      pp <- l

      while (!(pp %in% ancestors)) {
        pp <- geneology[pp, 2]
      }
      mrca_mat[global_to_local[k], global_to_local[l]] <- pp
      mrca_mat[global_to_local[l], global_to_local[k]] <- pp
    }
  }
  for (k in 1:(M - 1)) {
    # We need to convert to the leaf indexing (col 6)
    # to start in the correct place
    k_abs <- inds[k]
    for (l in (k + 1):M) {
      l_abs <- inds[l]
      # We need to convert to the leaf indexing (col 6)
      # to start in the correct place

      # Get the MRCA we want to stop at
      finally <- mrca_mat[k, l]

      # Trace the first arm back to the MRCA
      holder <- ind_par_mat[k_abs, 2]
      sumhold <- ind_par_mat[k_abs, 5]
      while (holder != finally) {
        sumhold <- sumhold + ind_par_mat[holder, 5]
        holder <- ind_par_mat[holder, 2]
      }

      # Trace the other arm back to the MRCA
      holder <- ind_par_mat[l_abs, 2]
      sumhold <- sumhold + ind_par_mat[l_abs, 5]
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

################################################################################
#
# Build distance matrix from sequences using Ape from a dataframe
#
################################################################################

#' Build a distance matrix from a dataframe of sequences using ape
#'
#' Use ape to compute a distance matrix using the specified distance model.
#' Default is the "TN93" model.
#'
#' @importFrom ape as.DNAbin
#' @export
build_distance_matrix_from_df <- function(df, model="TN93", keep_root = FALSE) {  # nolint: object_length_linter
  # Build a distance matrix from a dataframe using Ape
  # spike_root option adds a row and column for a node at the root

  if (!keep_root) {
    # If there is a row with name "root", remove it
    df <- df[df$name != "root", ]
  }
  # Convert the dataframe to a DNAbin through a list
  clean_data <- unlist(sapply(as.list(df$seq), tolower))
  clean_data <- sapply(clean_data, strsplit, "")
  seqs <- ape::as.DNAbin(clean_data)
  names(seqs) <- df$name[1:length(clean_data)]  # nolint: seq_linter
  # Todo: This is a hack. We should be able to do this without the 1:length
  # This should be controlled by the `keep_root` parameter

  distances <- ape::dist.dna(seqs, model = model, as.matrix = TRUE)
  # Rescale by sequence length to estimated # of mutations
  # The seqs are aligned already, so this is only the length of the first sequence
  distances <- distances * nchar(df$seq[1])

  return(distances)
}
