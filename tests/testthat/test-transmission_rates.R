
test_that("get_biphasic_HIV_rate() returns a function", {
    params <- list("R0" = 5)
    built_rate_fn <- get_biphasic_HIV_rate(params = params)
    expect_equal(rlang::is_callable(built_rate_fn), TRUE)
    }
)

test_that("get_biphasic_HIV_rate() evaluates correctly when R0 = 5", {
  params <- list("R0" = 5)
  built_rate_fn <- get_biphasic_HIV_rate(params = params)
  n <- 24
  birth_step <- rep(0, n)
  alive_step <- seq(from = 0, to = n - 1)
  rate_values <- built_rate_fn(current_step = alive_step, birth_step = birth_step)
  expect_equal(rate_values[1 : 3], rep(0.4/3 * 5 / .505, 3), testthat_tolerance())
  expect_equal(rate_values[4 : 24], rep(0.005 * 5 / .505, 21), testthat_tolerance())
})
