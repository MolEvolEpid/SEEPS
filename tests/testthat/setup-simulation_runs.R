# Here are fixtures for the simulator result
get_mocked_simulator_result_1 <- function() {
    set.seed(1947)  # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 3))

    simulator_result <- SEEPS::gen_transmission_history_exponential_constant(
        minimum_population = 5,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 10, total_steps = 0,
        spike_root = FALSE)
    return(simulator_result)
}

# mock 2 has small-intermediate scale
get_mocked_simulator_result_2 <- function() {
    set.seed(1947) # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 2)) # nolint: object_usage_linter

    simulator_result <- SEEPS::gen_transmission_history_exponential_constant( # nolint: object_usage_linter
        minimum_population = 50,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 200, total_steps = 24,
        spike_root = FALSE
    )
    return(simulator_result)
}

# mock 3 has intermediate-large scale
get_mocked_simulator_result_2 <- function() {
    set.seed(1947) # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 3)) # nolint: object_usage_linter

    simulator_result <- SEEPS::gen_transmission_history_exponential_constant( # nolint: object_usage_linter
        minimum_population = 500,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 5000, total_steps = 24,
        spike_root = FALSE
    )
    return(simulator_result)
}

# mock 4 has large scale. Do you really need this?
get_mocked_simulator_result_4 <- function() {
    set.seed(1947)  # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 4)) # nolint: object_usage_linter

    simulator_result <- SEEPS::gen_transmission_history_exponential_constant( # nolint: object_usage_linter
        minimum_population = 50,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 1000, total_steps = 120,
        spike_root = FALSE)
    return(simulator_result)
}
