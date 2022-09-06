test_that("Generate_simulator() works with toy examples", {

    # Get a simple GTR model
    rate_model <- get_simple_GTR_model()
    # Get a toy newick tree
    tree <- four_leaf_phylo()
    rngtools::setRNG(12)
    s <- .Random.seed
    # Call seq-gen
    fasta <- generate_sequences(phylogeny = tree$phylogeny,
                                root_sequence = "ACGT",
                                branch_rate = 1,
                                rng_seed = 1947,
                                rate_model = rate_model)

    expect_true(check_2line_fasta(fasta))
    # This should not have changed the seed
    expect_true(identical(s, .Random.seed))
    x <- runif(1)  # impact the seed here

    s <- .Random.seed
    # Again, this should not change the seed
    fasta2 <- generate_sequences(phylogeny = tree$phylogeny,
                                 root_sequence = "ACGT",
                                 branch_rate = 1,
                                 rng_seed = 1947,
                                 rate_model = rate_model)
    expect_true(identical(s, .Random.seed))
    # Now we should get a different value when we draw from the RNG
    y <- runif(1)
    expect_equal(fasta, fasta2)

    # These should be different
    expect_true(x != y)  # We should have gotten different values
})

test_that("Generate_simulator() works with a long sequence", {

    # Get a simple GTR model
    rate_model <- get_simple_GTR_model()
    # Get a toy newick tree
    tree <- four_leaf_phylo()
    # Get a poly-A string 100 chars long
    string <- strrep("A", 30000)
    # Call seq-gen
    fasta <- generate_sequences(phylogeny = tree$phylogeny,
                                root_sequence = string,
                                branch_rate = 1,
                                rng_seed = 1947,
                                rate_model = rate_model)
    expect_true(check_2line_fasta(fasta))
})


test_that("fasta_string_to_dataframe() works with seqgen outputs", {
    # Generate data
    # Get a simple GTR model
    rate_model <- get_simple_GTR_model()
    # Get a toy newick tree
    tree <- four_leaf_phylo()
    # Call seq-gen
    fasta <- generate_sequences(phylogeny = tree$phylogeny,
                                root_sequence = "ACGT",
                                branch_rate = 1,
                                rng_seed = 1947,
                                rate_model = rate_model)
    # Convert to a dataframe
    df <- fasta_string_to_dataframe(fasta)

})
