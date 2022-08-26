
# Fixture for a very simple GTR+I+G model
get_simple_GTR_model <-function() {  # nolint: object_name_linter
    model <- SEEPS::generate_rate_model(
        a2c = 1, a2g = 1,
        a2t = 1, c2g = 1,
        c2t = 1, g2t = 1,
        fa = 0.25, fc = 0.25,
        fg = 0.25, ft = 0.25,
        i = 0.5,
        alpha = 0.5, ncat = 2
    )
    return(model)
}


test_that("Generate_simulator() works with toy examples", {

    # Get a simple GTR model
    rate_model <- get_simple_GTR_model()
    # Get a toy newick tree
    tree <- four_leaf_phylo()
    # set.seed(12)
    rngtools::setRNG(12)
    s <- .Random.seed
    # Call seq-gen
    fasta <- generate_sequences(phylogeny = tree$phylogeny,
                                root_sequence = "ACGT",
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
                                rng_seed = 1947,
                                rate_model = rate_model)
    # Convert to a dataframe
    df <- fasta_string_to_dataframe(fasta)
    cat("new df is:\n")
    print(df)

})
