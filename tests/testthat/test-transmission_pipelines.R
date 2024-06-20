test_that("gen_transmission_history_exponential_constant", {
    set.seed(1947)  # for reproducability
    biphasic_rate_function <- get_biphasic_HIV_rate(list("R0" = 3))

    simulator_result <- gen_transmission_history_exponential_constant(
            minimum_population = 5,
            offspring_rate_fn = biphasic_rate_function,
            maximum_population_target = 10, total_steps = 0,
            spike_root = FALSE)
    expect_equal(dim(simulator_result$parents), c(120, 2))
    expect_identical(as.vector(simulator_result$parents[11:120, ]), rep(0, (110 * 2)))
    expect_identical(simulator_result$t_end, 6)
    expect_equal(simulator_result$active, 1:10)
    }
)

test_that("Check reduce_transmission_history(), 3 samples", {

    simulator_result <- get_mocked_simulator_result_1()
    sample <- c(7, 9, 10)
    geneology <- reduce_transmission_history(samples = sample,
            parents = simulator_result$parents,
            current_step = simulator_result$t_end)
    geneology <- geneology$geneology  # Unpack and get the object
    expect_equal(dim(geneology), c(6, 7))  # 7 columns
    expect_equal(geneology[1:5, 1], 1:5)
    expect_equal(geneology[6, 1], 0)
    expect_equal(geneology[, 6], c(1, 1, 1, 0, 0, 0))
    expect_equal(geneology[, 7], c(7, 9, 10, 0, 0, 0))
    expect_identical(table(geneology[, 2]), table(c(0, 0, 4, 4, 5, 5)))

})

test_that("Check reduce_transmission_history(), 4 samples", {

    simulator_result <- get_mocked_simulator_result_1()
    sample <- c(6, 7, 8, 10)
    geneology <- reduce_transmission_history(samples = sample,
            parents = simulator_result$parents,
            current_step = simulator_result$t_end)
    geneology <- geneology$geneology  # Unpack and get the object
    expect_equal(dim(geneology), c(8, 7))  # 7 columns
    expect_equal(geneology[1:7, 1], 1:7)
    expect_equal(geneology[8, 1], 0)
    expect_equal(geneology[, 6], c(1, 1, 1, 1, 0, 0, 0, 0))
    expect_equal(geneology[, 7], c(6, 7, 8, 10, 0, 0, 0, 0))
    expect_identical(table(geneology[, 2]), table(c(0, 0, 5, 5, 6, 6, 7, 7)))

})

test_that("Check bpb inputs", {
    simulator_result <- get_mocked_simulator_result_1()
    sample <- c(7, 9, 10)
    res <- reduce_transmission_history_bpb(samples = sample,
            parents = simulator_result$parents, current_step = simulator_result$t_end)

    expect_equal(res$parents, c(0, 1, 2, 1, 1, 3, 4, 5))
    expect_equal(res$transmission_times, c(0, 1, 2, 2, 3, 3, 4, 4))
    expect_equal(res$sample_times, c(3.001, 2.001, 3.001, 4.001, 4.001, 6, 6, 6))
    expect_equal(res$transformed_sample_indices, c(0, 0, 0, 0, 0, 7, 9, 10))
    # Now run BioPhyBreak codes
    phylogeny <- geneology_to_phylogeny_bpb(transmission_history = res$parents,
                infection_times = res$transmission_times,
                sample_times = res$sample_times,
                leaf_sample_ids = res$transformed_sample_indices)
})
