# Here are fixtures for the simulator result
get_mocked_simulator_result_1 <- function() {
    set.seed(1947) # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 3))

    simulator_result <- SEEPS::gen_transmission_history_exponential_constant(
        minimum_population = 5,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 10, total_steps = 0,
        spike_root = FALSE
    )
        spike_root = FALSE
    )
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
    set.seed(1947) # for reproducability
    set.seed(1947) # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 4)) # nolint: object_usage_linter

    simulator_result <- SEEPS::gen_transmission_history_exponential_constant( # nolint: object_usage_linter
        minimum_population = 50,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 1000, total_steps = 120,
        spike_root = FALSE
    )
        spike_root = FALSE
    )
    return(simulator_result)
}

# mock 1A has two samples, but is a small simulation
get_mocked_simulator_result_1A <- function() {
    set.seed(1948) # for reproducability
    set.seed(1948) # for reproducability
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(list("R0" = 4)) # nolint: object_usage_linter
    parameters <- SEEPS::wrap_parameters(
        offspring_rate_fn = biphasic_rate_function,
        minimum_population = 8,
        maximum_population_target = 15,
        total_steps = 0
    )
        total_steps = 0
    )
    state <- SEEPS::initialize(parameters)
    state <- SEEPS::gen_exp_phase(state, parameters)
    samples_1 <- SEEPS::contact_traced_uniform_restarts_ids(
        active = state$active,
        parents = state$parents,
        p = 0.5,
        minimum_sample_size = 3
    )
    samples_1 <- SEEPS::contact_traced_uniform_restarts_ids(
        active = state$active,
        parents = state$parents,
        p = 0.5,
        minimum_sample_size = 3
    )
    # Now tell the simulation to remove these from active
    time_1 <- state$curr_step
    # print(time_1)
    state <- SEEPS::remove_samples(state, samples_1$samples)
    # Run a few steps forward to get a new state
    state <- SEEPS::gen_const_phase(state, parameters, 5)
    samples_2 <- SEEPS::contact_traced_uniform_restarts_ids(
        active = state$active,
        parents = state$parents,
        p = 0.5,
        minimum_sample_size = 3
    )
    samples_2 <- SEEPS::contact_traced_uniform_restarts_ids(
        active = state$active,
        parents = state$parents,
        p = 0.5,
        minimum_sample_size = 3
    )
    time_2 <- state$curr_step
    # print(time_2)

    simulator_result <- state
    return(list(
        "result" = simulator_result, "t1" = time_1, "t2" = time_2,
        "s1" = samples_1, "s2" = samples_2
    ))
    return(list(
        "result" = simulator_result, "t1" = time_1, "t2" = time_2,
        "s1" = samples_1, "s2" = samples_2
    ))
}
