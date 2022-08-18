#' A helper function for tests. Builds a 3 leaf phylogeny with 2 internal nodes
#' and the newick tree. All internal distances are unit. Uses the _ suffix for node names.
three_leaf_phylo <- function() {
    phylogeny <- matrix(0, nrow=6, ncol=6)
    # Sample id's
    phylogeny[1:5, 1] <- 1:5
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
    return(list("phylogeny"=phylogeny, "newick_tree_mu"=newick_tree_mu, "newick_tree_mean"=newick_tree_mean))
}


test_that("phylogeny_to_newick() minimal examples", {
    reference <- three_leaf_phylo()
    phylogeny <- reference[["phylogeny"]]
    newick_tree_mu <- reference[["newick_tree_mu"]]
    newick_tree_mean <- reference[["newick_tree_mean"]]
    computed_tree <- phylogeny_to_newick(phylogeny, mode="mu")
    expect_identical(newick_tree_mu, computed_tree)
    computed_tree <- phylogeny_to_newick(phylogeny, mode="mean")
    expect_identical(newick_tree_mean, computed_tree)
})