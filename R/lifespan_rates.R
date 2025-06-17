# Lifespan distribution sampler functions
# This module provides examples of lifespan distribution sampler functions.
# Each function is a factory function that generates a sampler function.
# Each call of the sampler function generates a lifespan for each individual in the population.
# Each sampler function must accept three keyword arguments: `current_step`, `birth_step`, and `params`.
# 1. `current_step` <- (int) The current time step.
# 2. `birth_step` <- A vector of birth indices (integers).
# 3. `params` <- A list of optional parameters/functions needed to evaluate the rate function.
# Argument 3 may be substituted for a `...` in the function signature.
#
# One can also procedurally generate these functions within a factory.
# We provide several examples to demonstrate the API.
#
# To place more complexity (e.g. type or inheritance) in the choice of sampler,
# use the rate_calculator and assign_* factories to build samplers that can simulate
# more complex scenarios.


#' Simple uniform lifespan distribution sampler
#' @param vmin Minimum lifespan. Integer number of months.
#' @param vmax Maximum lifespan. Integer number of months.
#' @param params A list of optional parameters/functions needed to evaluate the rate function.
#' @export
#' @return A callable lifespan distribution sampler.
uniform_lifespan_sampler <- function(vmin, vmax, params = list()) {
    # Input checking
    if (!is.numeric(vmin) || length(vmin) != 1 || vmin %% 1 != 0) {
        stop("vmin must be a single integer value (integer).")
    }
    if (!is.numeric(vmax) || length(vmax) != 1 || vmax %% 1 != 0) {
        stop("vmax must be a single numeric value (integer).")
    }
    if (vmin >= vmax) {
        stop("vmin must be less than vmax.")
    }
    if (!is.list(params)) {
        stop("params must be a list.")
    }
    sampler <- function(current_step, birth_step, params = list()) {
        n <- length(birth_step)
        lifespans <- sample(vmin:vmax, n, replace = TRUE)
        return(lifespans)
    }
    return(sampler)
}

# exponential lifespan distribution sampler

#' Simple exponential lifespan distribution sampler
#' @param rate Rate parameter for the exponential distribution
#' @param params A list of optional parameters/functions needed to evaluate the rate function.
#' @export
#' @return A callable lifespan distribution sampler.
#' @examples
#' exponential_lifespan_sampler(0.1)
#' exponential_lifespan_sampler(0.1, list())
exponential_lifespan_sampler <- function(rate, params = list()) {
    # Input checking
    if (!is.numeric(rate) || length(rate) != 1) {
        stop("rate must be a single numeric value.")
    }
    if (rate <= 0) {
        stop("rate must be greater than 0.")
    }
    if (!is.list(params)) {
        stop("params must be a list.")
    }
    sampler <- function(current_step, birth_step, params = list()) {
        n <- length(birth_step)
        lifespans <- rexp(n, rate = rate)
        return(lifespans)
    }
    return(sampler)
}
