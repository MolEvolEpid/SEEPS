test_that("contact_tracing: contact_traced_uniform_ids() works", {
    # Get a simple set of parents
    parents <- t(matrix(c(0, 1, 1, 2, 3, 0, 1, 1, 2, 2),
                        nrow = 2, ncol = 5, byrow = TRUE))
    active <- c(1, 2, 3, 4, 5)
    # First check we recover everybody with p=1 discovery rate
    result <- contact_traced_uniform_ids(active = active,
                                         parents = parents,
                                         p = 1, minimum_sample_size = 5)
    # Check for success
    expect_equal(result[["success"]], TRUE)
    # Check we found everything
    expect_equal(result[["samples"]], c(1, 2, 3, 4, 5))
    # Expect we find everything
    expect_equal(result[["found"]], c(1, 2, 3, 4, 5))
    # Check completion status
    expect_equal(result[["status"]], "Out of nodes to explore")
    # Repeat for p=0.5
    set.seed(1947)  # with p=1, we didn't really need RNG, here we do
    result <- contact_traced_uniform_ids(active = active,
                                         parents = parents,
                                         p = 0, minimum_sample_size = 5)
    expect_equal(result, FALSE)

})

test_that("contact_tracing: test data retrieval", {
    parents <- t(matrix(c(0, 1, 1, 2, 3, 0, 1, 1, 2, 2),
                        nrow = 2, ncol = 5, byrow = TRUE))
    # Check an intermediate connection, with a defined
    # parent and offspring
    connections <- get_connections(seed = 3, parents = parents)
    expect_equal(connections[["parent"]], 1)
    expect_equal(connections[["infection_time"]], 1)
    expect_equal(connections[["secondary_infections"]], c(5))
    expect_equal(connections[["secondary_times"]], c(2))

    # Special cases:
    # Index case (1):
    connections <- get_connections(seed = 1, parents = parents)
    expect_equal(connections[["parent"]], 0)
    expect_equal(connections[["infection_time"]], 0)
    expect_equal(connections[["secondary_infections"]], c(2, 3))
    expect_equal(connections[["secondary_times"]], c(1, 1))

    # leaf - no offspring
    connections <- get_connections(seed = 4, parents = parents)
    expect_equal(connections[["parent"]], 2)
    expect_equal(connections[["infection_time"]], 2)
    # Key for no offpsrings is NULL
    expect_null(connections[["secondary_infections"]])
    expect_null(connections[["secondary_times"]])

})