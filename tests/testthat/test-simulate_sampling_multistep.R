test_that("check multiple time point sampling of bpb", {
    simulator_result <- get_mocked_simulator_result_1A()
    
    s1 <- simulator_result$s1
    expect_equal(s1$samples, c(4, 8, 10))
    t1 <- simulator_result$t1
    expect_equal(t1, 5)
    s2 <- simulator_result$s2
    expect_equal(s2$samples, c(3, 11, 7))
    t2 <- simulator_result$t2
    expect_equal(t2, 10)
    result <- simulator_result$result
    samples <- list(s1$samples, s2$samples)
    times <- list(t1, t2)
    # Backtrace the samples to the root
    # print(simulator_result$parents[, 1:15])
    phylo <- SEEPS::reduce_transmission_history_bpb2(samples, result$parents, times)
    
    expect_equal(phylo$transformed_sample_indices, c(0, 3, 4, 7, 8, 10, 11))
    expect_equal(phylo$sample_times, c(2.001, 10., 5., 10., 5., 5., 10.))

})
