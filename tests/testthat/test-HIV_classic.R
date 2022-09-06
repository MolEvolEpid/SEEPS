test_that("simulate_classic_HIV() pipeline", {
# Test example of simulator for previous paper
    params <- list("rate_function_parameters" = list("R0" = 5),
                   "minimum_population" = 15, "maximum_population_target" = 1000,
                   "total_steps_after_exp_phase" = 0, "mutation_rate" = 0.006 * 300)
    output <- simulate_classic_HIV(params = params)
    rownames(output$matrix) <- NULL
    colnames(output$matrix) <- NULL
    expect_equal(dim(output$matrix), c(15, 15))
    expect_equal(diag(output$matrix), rep(0, 15))
    expect_type(output$matrix, "double")
    # as.integer will flatten to a vector.
    expect_equal(all.equal(as.vector(output$matrix), as.integer(output$matrix)), TRUE)
})
