test_that("Classic and refactor simulations match", {
    sample_size <- 15
    params <- list(
        "rate_function_parameters" = list("R0" = 5),
        "minimum_population" = sample_size, "maximum_population_target" = sample_size * 2,
        "total_steps_after_exp_phase" = 0, "mutation_rate" = 0.009 * 300,
        "a" = 5, "b" = 5
    )

    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(
        params = params[["rate_function_parameters"]]
    )

    # Classic
    set.seed(1947) # Ensure we should perform the same simulation
    simulator_result_classic <- SEEPS::gen_transmission_history_exponential_constant_classic(
        minimum_population = params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE
    )

    # Refactored
    set.seed(1947) # Ensure we should perform the same simulation
    simulator_result <- SEEPS::gen_transmission_history_exponential_constant(
        minimum_population = params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE
    )

    # Use expect_equal to compare the two for equality.
    # Compare parents, t_end, active, total_offspring
    # print(simulator_result$parents[1:50,])
    # print("--- XXXX ---")
    # print(simulator_result_classic$parents[1:50,])

    expect_equal(simulator_result_classic$parents, simulator_result$parents)
    expect_equal(simulator_result_classic$t_end, simulator_result$t_end)
    expect_equal(simulator_result_classic$active, simulator_result$active)
    expect_equal(simulator_result_classic$total_offspring, simulator_result$total_offspring)
})
