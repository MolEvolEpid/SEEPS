test_that("contact_tracing: contact_traced_uniform_ids() works", {
    # Get a simple set of parents
    parents <- t(matrix(c(0, 1, 1, 2, 3, 0, 1, 1, 2, 2),
        nrow = 2, ncol = 5, byrow = TRUE
    ))
    active <- c(1, 2, 3, 4, 5)
    # First check we recover everybody with  =  discovery rate
    result <- contact_traced_uniform_ids(
        active = active,
        parents = parents,
        p = 1, minimum_sample_size = 5
    )
    # Check for success
    expect_equal(result[["success"]], TRUE)
    # Check we found everything
    expect_equal(result[["samples"]], c(1, 2, 3, 4, 5))
    # Expect we find everything
    expect_equal(sort(result[["found"]]), c(1, 2, 3, 4, 5))
    # found list won't be sorted

    # Check completion status
    expect_equal(result[["status"]], "Out of nodes to explore")
    # Repeat for  = .5
    set.seed(1947) # with  = , we didn't really need RNG, here we do
    result <- contact_traced_uniform_ids(
        active = active,
        parents = parents,
        p = 0, minimum_sample_size = 5
    )
    expect_equal(result, FALSE)
})

test_that("contact_tracing: contact_traced_uniform_restarts_ids() works", {
    # Get a simple set of parents
    parents <- t(matrix(c(0, 1, 1, 2, 3, 0, 1, 1, 2, 2),
        nrow = 2, ncol = 5, byrow = TRUE
    ))
    active <- c(1, 2, 3, 4, 5)
    # First check we recover everybody with  =  discovery rate
    result <- contact_traced_uniform_restarts_ids(
        active = active,
        parents = parents,
        p = 1,
        minimum_sample_size = 5
    )
    # Check for success
    expect_equal(result[["success"]], TRUE)
    # Check we found everything
    expect_equal(result[["samples"]], c(1, 2, 3, 4, 5))
    # Expect we explored everything
    expect_equal(sort(result[["found"]]), c(1, 2, 3, 4, 5))
    # found list won't be sorted

    # Check completion status
    expect_equal(result[["status"]], "Found enough data with 1 attempt(s)")
    # Repeat for  = .5
    set.seed(1947) # with  = , we didn't really need RNG, here we do
    result <- contact_traced_uniform_restarts_ids(
        active = active,
        parents = parents,
        p = 0.05,
        minimum_sample_size = 5
    )
    # Expect we succeeded
    expect_equal(result[["success"]], TRUE)
    # Expect that we found everything
    expect_equal(sort(result[["samples"]]), c(1, 2, 3, 4, 5))
    # Expect we explored everything
    expect_equal(sort(result[["found"]]), c(1, 2, 3, 4, 5))
    # Check group IDs
    expect_equal(result[["group_ids"]], c(1, 1, 2, 3, 4))
})

test_that("contact_tracing: test data retrieval", {
    parents <- t(matrix(c(0, 1, 1, 2, 3, 0, 1, 1, 2, 2),
        nrow = 2, ncol = 5, byrow = TRUE
    ))
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

# Next check that contact tracing works on a real example
# We'll use get_mocked_simulator_result_1 for data
# and then check that we can recover the index case
# and all of its secondary infections
test_that("Contact tracing on simulation run with restarts works", {
    data <- get_mocked_simulator_result_1()
    active <- data[["active"]]
    parents <- data[["parents"]]
    # Call the contact tracing
    set.seed(1947)
    result <- contact_traced_uniform_restarts_ids(
        active = active,
        parents = parents,
        p = 1,
        minimum_sample_size = 5
    )
    # Check for success
    expect_equal(result[["group_ids"]], c(1, 1, 1, 1, 1))
    expect_equal(result[["success"]], TRUE)
    # Check we found five individuals
    expect_equal(length(result[["samples"]]), 5)
    expect_equal(result[["status"]], "Found enough data with 1 attempt(s)")
})
