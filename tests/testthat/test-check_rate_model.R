test_that("check_rate_model() functions correctly", {
    model <- generate_rate_model(
        a2c = 1, a2g = 1, a2t = 1, c2g = 1, c2t = 1,
        g2t = 2,
        fa = 0.25, fc = 0.25, fg = 0.25, ft = 0.25,
        i = 0.1, alpha = 0.25, ncat = 8
    )
    result <- check_rate_model(model)
    expect_true(result)
    # Now failure modes
    model$a2c <- -1
    expect_error(check_rate_model(model), "All rate model parameters must be positive")
    model$a2c <- 1
    model$i <- -0.5
    expect_error(check_rate_model(model), "All rate model parameters must be positive")
    model$i <- 2
    expect_error(check_rate_model(model), "i must be strictly less than 1")
    model$i <- 0.5
    model$ncat <- 0.5
    expect_error(check_rate_model(model), "ncat must be an integer")
    model$fa <- 0.1
    expect_error(check_rate_model(model), "Fractions of each base must sum to 1")
})
