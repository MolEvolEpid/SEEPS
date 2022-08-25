test_that("add_root_to_newick() works with toy examples", {
    # Get a toy newick tree
    tree <- four_leaf_phylo()
    new_tree <- add_root_to_newick(tree$newick_tree_mean)
    expect_equal(new_tree,
                 "(((1_:1.1,2_:1.1):1.2,(3_:1.1,4_:1.1):1.2):0,root:0);")
})
