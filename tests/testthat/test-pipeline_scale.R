# Here is a factory for the simulator result
# mock 2 has small-intermediate scale
get_mocked_simulator_result_2 <- function() {
    set.seed(1947) # for reproducability
    biphasic_rate_function <- get_biphasic_HIV_rate(list("R0" = 2)) # nolint: object_usage_linter

    simulator_result <- gen_transmission_history_exponential_constant( # nolint: object_usage_linter
        minimum_population = 50,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 200, total_steps = 24,
        spike_root = FALSE
    )
    return(simulator_result)
}

# Here is a factory for the simulator result
# mock 3 has intermediate-large scale
get_mocked_simulator_result_2 <- function() {
    set.seed(1947) # for reproducability
    biphasic_rate_function <- get_biphasic_HIV_rate(list("R0" = 3)) # nolint: object_usage_linter

    simulator_result <- gen_transmission_history_exponential_constant( # nolint: object_usage_linter
        minimum_population = 500,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 5000, total_steps = 24,
        spike_root = FALSE
    )
    return(simulator_result)
}

test_that("Check bpb inputs for 50 tips", {
    simulator_result <- get_mocked_simulator_result_2()
    sample <- 276:325
    res <- reduce_transmission_history_bpb(
        samples = sample,
        parents = simulator_result$parents,
        current_step = simulator_result$t_end
    )
    # Now run BioPhyBreak codes
    phylogeny <- geneology_to_phylogeny_bpb(
        transmission_history = res$parents,
        infection_times = res$transmission_times,
        sample_times = res$sample_times
    )
})




skip("This test is for testing large simulations.")
test_that("Check bpb inputs for 500 tips", {
    simulator_result <- get_mocked_simulator_result_4()
    n <- length(simulator_result$active)
    sample <- simulator_result$active[n - 500:n]
    res <- reduce_transmission_history_bpb(
        samples = sample,
        parents = simulator_result$parents, t_end = simulator_result$t_end
    )

    # Now run BioPhyBreak codes
    phylogeny <- geneology_to_phylogeny_bpb(
        transmission_history = res$parents,
        infection_times = res$transmission_times,
        sample_times = res$sample_times
    )
})
