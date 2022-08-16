# Here is a factory for the simulator result
# mock 4 has large scale. You may want to skip this test.
skip("For profiling only")
get_mocked_simulator_result_4 <- function() {
    set.seed(1947)  # for reproducability
    biphasic_rate_function <- get_biphasic_HIV_rate(list("R0" = 4))

    simulator_result <- gen_transmission_history_exponential_constant(
        minimum_population = 50,
        offspring_rate_fn = biphasic_rate_function,
        maximum_population_target = 1000, total_steps = 120,
        spike_root = FALSE)
    return(simulator_result)
}

skip("For profiling only")
test_that("Check bpb inputs for 50 tips", {
    simulator_result <- get_mocked_simulator_result_4()
    # print(simulator_result$active)
    n <- length(simulator_result$active)
    sample <- simulator_result$active[n-50:n]
    res <- reduce_transmission_history_biophybreak(samples = sample,
                                                   parents = simulator_result$parents,
                                                   current_step = simulator_result$t_end)

    phylogeny <- geneology_to_phylogeny_bpb(transmission_history = res$parents,
                                            infection_times = res$transmission_times,
                                            sample_times = res$sample_times)
})
