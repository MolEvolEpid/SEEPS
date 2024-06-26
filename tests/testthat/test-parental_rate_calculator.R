test_that("Test for parental_rate_calculator", {
    biphasic_rate_function1 <- SEEPS::get_biphasic_HIV_rate(list("R0" = 1.1))
    biphasic_rate_function2 <- SEEPS::get_biphasic_HIV_rate(list("R0" = 4))
    # Make two rate functions,
    offspring_rate_fn <- c(biphasic_rate_function1, biphasic_rate_function2)
    # Set the rate function to be inherited from the parent
    offspring_selector_fn <- SEEPS::rate_assign_parental_factory(offspring_rate_fn)

    minimum_population <- 10
    maximum_population_target <- 15
    total_steps <- 8

    parameters <- wrap_parameters(minimum_population,
                                  offspring_rate_fn,
                                  maximum_population_target,
                                  total_steps,
                                  rate_selector_function = offspring_selector_fn)
    # parameters$
    state <- initialize(parameters)
    # Run the simulation
    state <- step(state, parameters)
    state$rate_indexing[2] <- 2
    for (i in 1:total_steps) {
        state <- step(state, parameters)
    }
})
