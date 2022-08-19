#' Convert a phylogeny array to a newick tree
#'
#' Build a newick tree. Nodes are named "#_".
#' Specify the mode argument to select the branch length mode. Mode "mu" denotes
#' a (sampled) mutation count, while mode "mean" denotes expected distances.
#' @export
phylogeny_to_newick <- function(phylogeny, mode="mu") {
    # Need to validate inputs
    # We won't trace back the entire phylogeny, just the marked "leaf" nodes.
    col <- -1
    error_msg <- "Mode argument could not be matched to a case."
    if (mode == "mu") col <- 5
    if (mode == "mean") col <- 4
    if (col == -1) throw(error_msg)

    M <- sum(phylogeny[, 6]) # number of leaf nodes
    inds <- which(phylogeny[, 6] != 0) # Get the index of the nonzeros
    n_nodes <- dim(phylogeny)[1] - 1

    # We need to convert global indexing (1:# of tracked infections) to local indexing (1:M)
    global_to_local <- rep(0, max(inds))
    global_to_local[inds] <- seq_along(M) # So g_2_l[global] = local for O(1) lookup
    branches <- character(n_nodes)  # Sub-tree at each branch
    # Fill in the names of the leaf nodes
    names <- lapply(1:M, function(x) paste0(x, "_"))
    branches[1:M] <- names
    merge_nodes <- vector("list", length = n_nodes)

    for (index in 1:(n_nodes-1)) {
        parent <- phylogeny[index, 2]  # Get the parent
        if (is.null(merge_nodes[[parent]])) { # First offspring
            merge_nodes[[parent]] <- index
        } else { # We alredy found one of the offspring
            # Record second offspring nodes
            children <- c(merge_nodes[[parent]], index)
            # Now construct the merged string
            branches[parent] <- paste0("(", branches[children[1]], ":", phylogeny[children[1], col], ",",
                                            branches[children[2]], ":", phylogeny[children[2], col], ")")
        }

    }

    # Complete the tree
    paste0(branches[n_nodes], ";")
}
