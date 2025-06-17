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
    sample_times = res$sample_times,
    leaf_sample_ids = res$transformed_sample_indices
  )
})




skip("For profiling only")
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
    sample_times = res$sample_times,
    leaf_sample_ids = res$transformed_sample_indices
  )
})
