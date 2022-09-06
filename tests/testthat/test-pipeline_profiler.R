skip("For profiling only")
test_that("Check bpb inputs for 50 tips", {
    simulator_result <- get_mocked_simulator_result_4()
    n <- length(simulator_result$active)
    sample <- simulator_result$active[(n - 50):n]  # take the last 50 tips
    res <- reduce_transmission_history_bpb(samples = sample,
                                           parents = simulator_result$parents,
                                           current_step = simulator_result$t_end)

    phylogeny <- geneology_to_phylogeny_bpb(transmission_history = res$parents,
                                            infection_times = res$transmission_times,
                                            sample_times = res$sample_times,
                                            leaf_sample_ids = res$transformed_sample_indices)
})
