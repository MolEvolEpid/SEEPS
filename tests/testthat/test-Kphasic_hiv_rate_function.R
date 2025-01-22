test_that("Kphasic hiv rate function works - simple", {
    # Setup some simple parameters

    # Relative rates. Will be normalized so it integrates to 1.
    rate_list <- list(c(1, 1), c(2, 1), c(1, 1))
    params <- list(R0 = 1)
    rate_fn <- get_Kphasic_hiv_rate_function(
        rate_list = rate_list,
        params = params
    )
    # Should have 1/4, 1/2, 1/4 as the rates, as R0 = 1
    birth_times <- 0:4
    print(birth_times)
    rates <- rate_fn(4, birth_times)
    print(rates)
    expect_equal(rates, c(0., 0., 0.25, 0.5, 0.25))
})
