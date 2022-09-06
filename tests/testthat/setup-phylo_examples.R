#' A helper function for tests. Builds a 3 leaf phylogeny with 2 internal nodes
#' and the newick tree. All internal distances are unit. Uses the _ suffix for node names.
three_leaf_phylo <- function() {
    phylogeny <- matrix(0, nrow = 6, ncol = 7)
    # Sample id's
    phylogeny[1:5, 1] <- 1:5
    # absolute id's
    phylogeny[1:3, 7] <- 4:6
    # Contact network
    phylogeny[1:2, 2] <- 4
    phylogeny[3:4, 2] <- 5
    # Sample times
    phylogeny[1:3, 3] <- 3
    phylogeny[4,   3] <- 2
    phylogeny[5,   3] <- 1
    # Set Expected # of mutations based on time/rates
    phylogeny[1:3, 4] <- 1.1
    phylogeny[4, 4] <- 1.2
    # Actual amout of mutations
    phylogeny[1:4, 5] <- 1
    # leaf flags
    phylogeny[1:3, 6] <- 1

    newick_tree_mu <- "(3_:1,(1_:1,2_:1):1);"
    newick_tree_mean <- "(3_:1.1,(1_:1.1,2_:1.1):1.2);"
    newick_tree_mu_abs <- "(6_:1,(4_:1,5_:1):1);"
    newick_tree_mean_abs <- "(6_:1.1,(4_:1.1,5_:1.1):1.2);"
    return(list("phylogeny" = phylogeny, "newick_tree_mu" = newick_tree_mu,
                "newick_tree_mean" = newick_tree_mean,
                "newick_tree_mu_abs" = newick_tree_mu_abs,
                "newick_tree_mean_abs" = newick_tree_mean_abs))
}

#' A helper function for tests. Builds a 4 leaf phylogeny with 2 internal nodes
#' and the newick tree. All internal distances are unit. Uses the _ suffix for node names.
four_leaf_phylo <- function() {
    phylogeny <- matrix(0, nrow = 8, ncol = 7)
    # Sample id's
    phylogeny[1:7, 1] <- 1:7
    # absolute id's
    phylogeny[1:4, 7] <- 5:8
    # Contact network
    phylogeny[1:2, 2] <- 5
    phylogeny[3:4, 2] <- 6
    phylogeny[5:6, 2] <- 7
    # Sample times
    phylogeny[1:4, 3] <- 3
    phylogeny[4:5, 3] <- 2
    phylogeny[6,   3] <- 1
    # Set Expected # of mutations based on time/rates
    phylogeny[1:4, 4] <- 1.1

    phylogeny[5:6, 4] <- 1.2
    # Actual amout of mutations
    phylogeny[1:6, 5] <- 1
    # leaf flags
    phylogeny[1:4, 6] <- 1

    newick_tree_mu <- "((1_:1,2_:1):1,(3_:1,4_:1):1);"
    newick_tree_mean <- "((1_:1.1,2_:1.1):1.2,(3_:1.1,4_:1.1):1.2);"
    newick_tree_mu_abs <- "((5_:1,6_:1):1,(7_:1,8_:1):1);"
    newick_tree_mean_abs <- "((5_:1.1,6_:1.1):1.2,(7_:1.1,8_:1.1):1.2);"
    return(list("phylogeny" = phylogeny, "newick_tree_mu" = newick_tree_mu,
                "newick_tree_mean" = newick_tree_mean,
                "newick_tree_mu_abs" = newick_tree_mu_abs,
                "newick_tree_mean_abs" = newick_tree_mean_abs))
}
