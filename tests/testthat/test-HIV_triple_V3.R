test_that("simulate_triple_HIV_V3() pipeline", {
    # Test example of simulator with bpb
    set.seed(1947)
    sample_size <- 15
    # Expect 33 out, full matrix
    params <- list(
        "rate_function_parameters" = list("R0" = 5),
        "minimum_population" = sample_size,
        "maximum_population_target" = sample_size * 10,
        "contact_tracing_discovery_probability" = 0.5,
        "total_steps_after_exp_phase" = 10, "mutation_rate" = 0.0067,
        # Use mutation rate for V3 from Leitner & Albert, 1999.
        "a" = 5, "b" = 5
    )
    output <- simulate_all_paradigms_HIV_V3(params = params)
    # Expect output$seqeunces is a dataframe
    expect_true(is.data.frame(output$sequences))
    # Check that each sequence is a valid DNA sequence
    expect_equal(dim(output[["matrix_seqs"]]), c(sample_size, sample_size))
    # Let's clear row/col names so we can compare elements safely
    rownames(output[["matrix_seqs"]]) <- NULL
    colnames(output[["matrix_seqs"]]) <- NULL
    expect_equal(diag(output[["matrix_seqs"]]), rep(0, sample_size))
    expect_type(output[["matrix_seqs"]], "double")
    # We don't expect integer distances to be output, as distance corrections
    # are occurring with TN93

  # Test matrix_trans
  rownames(output[["matrix_trans"]]) <- NULL
  colnames(output[["matrix_trans"]]) <- NULL
  expect_equal(dim(output[["matrix_trans"]]), c(sample_size, sample_size))
  expect_equal(diag(output[["matrix_trans"]]), rep(0, sample_size))
  expect_type(output[["matrix_trans"]], "double")

    # Test matrix_phylo
    rownames(output[["matrix_phylo"]]) <- NULL
    colnames(output[["matrix_phylo"]]) <- NULL
    expect_equal(dim(output[["matrix_phylo"]]), c(sample_size, sample_size))
    expect_equal(diag(output[["matrix_phylo"]]), rep(0, sample_size))
    expect_type(output[["matrix_phylo"]], "double")
})
