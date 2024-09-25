skip(message = "Reserved for future use")
test_that("Kphasic hiv rate function works - simple", {
    # Setup some simple parameters
    rate_list <- list(c(1, 1), c(2, 1), c(1, 1))
    target_length <- 3
    params <- list(R0 = 1)
    rate_fn <- get_Kphasic_hiv_rate_function(rate_list = rate_list,
                                             target_length = target_length,
                                             params = params)
    # Should have 1/4, 1/2, 1/4 as the rates
    rates <- lapply(c(1, 2, 3), rate_fn)
    expect_equal(rates, c(0.25, 0.5, 0.25))
})
