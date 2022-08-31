test_that("simulate_modern_HIV() pipeline", {
# Test example of simulator with bpb
    set.seed(1947)
    sample_size <- 15
    params <- list("rate_function_parameters" = list("R0" = 5),
                   "minimum_population" = sample_size, "maximum_population_target" = sample_size * 10,
                   "total_steps_after_exp_phase" = 0, "mutation_rate" = 0.009 * 300,
                   "a" = 5, "b" = 5)
    output <- simulate_modern_HIV(params = params)
    expect_equal(dim(output$matrix), c(sample_size, sample_size))
    expect_equal(diag(output$matrix), rep(0, sample_size))
    expect_type(output$matrix, "double")
    # as.integer will flatten to a vector.
    expect_equal(all.equal(as.vector(output$matrix), as.integer(output$matrix)), TRUE)
})
