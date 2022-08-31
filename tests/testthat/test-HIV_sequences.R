test_that("simulate_sequences_HIV_V3() pipeline", {
# Test example of simulator with bpb
    set.seed(1947)
    sample_size <- 15
    # Expect 33 out, full matrix
    params <- list("rate_function_parameters" = list("R0" = 5),
                   "minimum_population" = sample_size,
                   "maximum_population_target" = sample_size * 10,
                   "contact_tracing_discovery_probability" = 0.9,
                   "total_steps_after_exp_phase" = 0, "mutation_rate" = 0.009,
                   "a" = 5, "b" = 5)
    output <- simulate_sequences_HIV_V3(params = params)
    # print(names(output))
    # Expect output$seqeunces is a dataframe
    expect_true(is.data.frame(output$sequences))
    # Check that each sequence is a valid DNA sequence
    expect_equal(dim(output[["matrix"]]), c(sample_size, sample_size))
    # Let's clear row/col names so we can compare elements safely
    rownames(output[["matrix"]]) <- NULL
    colnames(output[["matrix"]]) <- NULL
    # print(rep(0, sample_size))
    expect_equal(diag(output[["matrix"]]), rep(0, sample_size))
    expect_type(output[["matrix"]], "double")
    # We don't expect integer distances to be output, as distance corrections
    # are occurring with TN93

})
