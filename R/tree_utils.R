#' Convert a phylogeny array to a newick tree
#'
#' Build a newick tree from a phylogeny or geneology array.
#'  Nodes are named "#_". Specify the mode argument to select the branch
#' length mode. Mode "mu" denotes a (sampled) mutation count, while mode
#' "mean" denotes expected distances.
#'
#' @param phylogeny A phylogeny or geneology array.
#' @param mode String to determine reconstruction mode. Default is "mu", for
#'   mutation count. Alternative "mean" for expected distance.
#' @export
phylogeny_to_newick <- function(phylogeny, mode = "mu") {
    # Need to validate inputs
    # We won't trace back the entire phylogeny, just the marked "leaf" nodes.
    col <- -1
    error_msg <- "Mode argument could not be matched to a case. Expected 'mu' or 'mean'."
    if (mode == "mu") col <- 5
    if (mode == "mean") col <- 4
    if (col == -1) stop(error_msg)

    M <- sum(phylogeny[, 6]) # number of leaf nodes
    inds <- which(phylogeny[, 6] != 0) # Get the index of the nonzeros
    n_nodes <- dim(phylogeny)[1] - 1

    # We need to convert global indexing (1:# of tracked infections)
    # to local indexing (1:M)
    global_to_local <- rep(0, max(inds))
    global_to_local[inds] <- seq_along(M)
    # So g_2_l[global] = local for O(1) lookup
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
            branches[parent] <- paste0("(", branches[children[1]], ":",
                                       phylogeny[children[1], col], ",",
                                       branches[children[2]], ":",
                                       phylogeny[children[2], col], ")")
        }

    }

    # Complete the tree
    tree <- paste0(branches[n_nodes], ";")
    return(tree)
}

#' Add a root node to a newick tree string with distance 0 from the root
#'
#' A utility function for placing an explicit "root" into a newick tree string.
#' at the implicit root. This function performs: `(tree); -> (tree:0, root:0);`.
#'
#' @param tree A newick tree string to add a root node to.
#' @param root_name The name of the root node to add. Default is "root".
#' @seealso phylogeny_to_newick
#' @export
add_root_to_newick <- function(tree, root_name = "root") {
    subtree <- substring(tree, 1, nchar(tree) - 1)
    new_tree <- paste0("(", subtree, ":0,", root_name, ":0);")

    return(new_tree)
}
