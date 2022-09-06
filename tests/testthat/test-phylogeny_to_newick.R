test_that("phylogeny_to_newick() argument check works", {
    reference <- three_leaf_phylo()
    phylogeny <- reference[["phylogeny"]]
    expect_error(phylogeny_to_newick(phylogeny, mode = "muan"),
      "Mode argument could not be matched to a case. Expected 'mu' or 'mean'.")
})

test_that("phylogeny_to_newick() three leaf example", {
    reference <- three_leaf_phylo()
    phylogeny <- reference[["phylogeny"]]
    newick_tree_mu <- reference[["newick_tree_mu"]]
    newick_tree_mean <- reference[["newick_tree_mean"]]
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mu")
    expect_identical(newick_tree_mu, computed_tree)
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mean")
    expect_identical(newick_tree_mean, computed_tree)
})

test_that("phylogeny_to_newick() four leaf example", {
    reference <- four_leaf_phylo()
    phylogeny <- reference[["phylogeny"]]
    newick_tree_mu <- reference[["newick_tree_mu"]]
    newick_tree_mean <- reference[["newick_tree_mean"]]
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mu")
    expect_identical(newick_tree_mu, computed_tree)
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mean")
    expect_identical(newick_tree_mean, computed_tree)
})

test_that("phylogeny_to_newick() three leaf example with abs indexing", {
    reference <- three_leaf_phylo()
    phylogeny <- reference[["phylogeny"]]
    newick_tree_mu <- reference[["newick_tree_mu_abs"]]
    newick_tree_mean <- reference[["newick_tree_mean_abs"]]
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mu", label_mode = "abs")
    expect_identical(newick_tree_mu, computed_tree)
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mean", label_mode = "abs")
    expect_identical(newick_tree_mean, computed_tree)
})

test_that("phylogeny_to_newick() four leaf example with abs indexing", {
    reference <- four_leaf_phylo()
    phylogeny <- reference[["phylogeny"]]
    newick_tree_mu <- reference[["newick_tree_mu_abs"]]
    newick_tree_mean <- reference[["newick_tree_mean_abs"]]
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mu", label_mode = "abs")
    expect_identical(newick_tree_mu, computed_tree)
    computed_tree <- phylogeny_to_newick(phylogeny, mode = "mean", label_mode = "abs")
    expect_identical(newick_tree_mean, computed_tree)
})
