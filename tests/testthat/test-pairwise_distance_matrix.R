test_that("Check pairwise distance matrix construction, 3 samples", {

    simulator_result <- get_mocked_simulator_result_1()
    sample <- c(7, 9, 10)
    geneology <- reduce_transmission_history(samples = sample,
            parents = simulator_result$parents,
            current_step = simulator_result$t_end)
    geneology <- geneology$geneology  # Unpack and get the object
    set.seed(1947)
    geneology <- stochastify_transmission_history(geneology, rate = 4)
    # unpack
    geneology <- geneology$geneology
    distance_matrix <- geneology_to_distance_matrix(geneology = geneology)
    # Check we got the correct rownames
    expect_equal(rownames(distance_matrix), c("7", "9", "10"))
    expect_equal(colnames(distance_matrix), c("7", "9", "10"))
    expect_equal(base::isSymmetric.matrix(distance_matrix), TRUE)
    # Check values against manual construction
    # Expect row names of 7,9,10, with distances of 7:9, 7:10, 9:10
    expected_matrix <- c(0, 15, 18, 15, 0, 21, 18, 21, 0)
    expected_matrix <- matrix(expected_matrix, nrow = 3, ncol = 3)
    rownames(expected_matrix) <- c(7, 9, 10)
    colnames(expected_matrix) <- c(7, 9, 10)
    expect_equal(distance_matrix, expected_matrix)

    # Now set the names aside so we can compare values
    rownames(distance_matrix) <- NULL
    colnames(distance_matrix) <- NULL
    # This next test does consider machine tolerances
    # print(distance_matrix)
    expect_equal(dim(distance_matrix), c(3, 3))
    expect_equal(diag(distance_matrix), rep(0, 3))
})

test_that("Check pairwise distance matrix construction, 3 samples", {


    simulator_result <- get_mocked_simulator_result_1()
    sample <- c(6, 7, 8, 10)
    geneology <- reduce_transmission_history(samples = sample,
            parents = simulator_result$parents,
            current_step = simulator_result$t_end)
    geneology <- geneology$geneology  # Unpack and get the object
    set.seed(1947)
    geneology <- stochastify_transmission_history(geneology, rate = 4)
    # unpack
    geneology <- geneology$geneology
    # print(geneology)
    distance_matrix <- geneology_to_distance_matrix_classic(geneology = geneology)
    # print(distance_matrix)
    # print("ZZ")
    distance_matrix <- geneology_to_distance_matrix(geneology = geneology)
    # print("AA")
    # print(distance_matrix)
    # print(rownames(distance_matrix))
})
