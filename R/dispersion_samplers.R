# Dispersion samplers
#
# All dispersion samplers need to accept the following arguments:
# 1. means <- A vector of means for each individual. means[[i]] is the expected number of offspring for individual i.
# 2. current_step <- The current time step.
# 3. birth_step <- A vector of birth indices (integers).
# 4. ... <- Optional parameters/functions needed to evaluate the rate function.
#
# Many methods exist for calculating dispersion. However, we provide three examples to demonstrate how the API works.
#
# 1. The simplest method is to use an exponential distribution.
# 2. The second method is to use a Poisson distribution. More realistic, and still 1 parameter
# 3. The third method is to use a negative binomial distribution. More realistic, and 2 parameters
#
# Future support for Warring, Simon-Yule, and other distributions is planned.
# Users are free to implement their own dispersion samplers.

#' Add dispersion to point estimates.
#'
#' @param means A vector of means for each individual. means[[i]] is the expected number of offspring for individual i.
#' @param current_step The current time step.
#' @param birth_step A vector of birth indices (integers).
#' @param dispersion_type A vector of integers indicating the type of dispersion to use for each individual.
#' @param samplers A list of callable dispersion samplers.
#' @param params A list of optional parameters/functions needed to evaluate the rate function(s).
#' @return A vector of dispersions for each individual.
#' @export
dispersion_calculator <- function(means, current_step, birth_step, dispersion_type, samplers, params = list()) {
    # Input checking

    # first, check that mean has length >0
    if (length(means) == 0) {
        return(means)
    }
    if (any(!is.numeric(means)) || any(means < 0)) {
        stop("means must be a vector of non-negative values.")
    }
    # check that the max dispersion type is less than the number of samplers
    if (!is.numeric(dispersion_type) || any(dispersion_type < 1) || any(dispersion_type > length(samplers))) {
        stop("dispersion_type must be a vector of integers between 1 and the number of samplers.")
    }

    # main functions
    values <- rep(-1, length(means))

    for (i in seq_along(samplers)) {
        sampler <- samplers[[i]]
        mask <- dispersion_type == i
        values[mask] <- sampler(means[mask], current_step, birth_step[mask], params)
    }
    return(values)
}

# Constant (Null) dispersion sampler
# No variance, no sampling. Just returns the mean.

#' Constant (Null) dispersion sampler
#'
#' Returns the mean, and does not apply any dispersion correction/sampling. Useful for testing and minimal examples.
#'
#' @param means Mean number of infections.
#' @param current_step The current time step. Ignored.
#' @param birth_step A vector of birth indices (integers). Ignored.
#' @param params A list of optional parameters/functions needed to evaluate the rate function. Ignored.
#' @return A vector of means.
#' @export
constant_dispersion_sampler <- function(means, current_step, birth_step, params = list()) {
    # Input checking
    if (!is.numeric(means) || any(means <= 0)) {
        stop("means must be a vector of positive values.")
    }
    return(means)
}

# Exponential dispersion sampler
# This is the simplest dispersion sampler. It only requires a rate parameter.

#' Exponential dispersion sampler
#'
#' Calculates a dispersion correction based on an expected number of offspring.
#' Sample an exponential distribution with rate equal to the mean, to determine
#' an expected number of contributions to the total number of infections.
#'
#' @param params A list of optional parameters/functions needed to evaluate the rate function. Ignored.
#'
#' @importFrom stats rexp
#' @return A callable dispersion sampler. Accepts arguments `means`,
#' `current_step`, and `birth_step`.
#' @export
exponential_dispersion_sampler <- function(params = list()) {
    sampler <- function(means, current_step, birth_step, ...) {
        # Input checking
        if (!is.numeric(means) || any(means <= 0)) {
            stop("means must be a vector of positive values.")
        }
        n <- length(means)
        return(rexp(n, rate = means))
    }
    return(sampler)
}

# Poisson dispersion sampler
# This is a more realistic dispersion sampler. It only requires a rate parameter.

#' Simple Poisson dispersion sampler
#' @param params A list of optional parameters/functions needed to evaluate the rate function. Ignored.
#' @importFrom stats rpois
#' @return A callable dispersion sampler. Accepts arguments `means`,
#' `current_step`, and `birth_step`.
#'
#' @export
poisson_dispersion_sampler <- function(params = list()) {
    sampler <- function(means, current_step, birth_step, ...) {
        # Input checking
        if (!is.numeric(means) || any(means <= 0)) {
            stop("means must be a vector of positive values.")
        }
        n <- length(means)
        return(rpois(n, lambda = means))
    }
    return(sampler)
}

# Negative binomial dispersion sampler
# This is a more realistic dispersion sampler. It requires a rate parameter and a dispersion parameter.

#' Simple negative binomial dispersion sampler
#'
#' Calculates a dispersion correction based on an expected number of offspring,
#' using a negative binomial model. This is a more realistic model than the exponential
#' as it allows for variance to be tuned independent of means.
#'
#' @param params A list of optional parameters/functions needed to evaluate the rate function. Must contain a dispersion
#'  parameter `dispersion`. Required.
#' @importFrom stats rnbinom
#' @return A callable dispersion sampler. Accepts arguments `means`,
#' `current_step`, and `birth_step`. Dispersion parameter is set in the `params`
#' argument at initialization.
#'
#' @export
negative_binomial_dispersion_sampler <- function(params = list("dispersion" = 1)) {
    # Input checking
    if (!is.list(params) || !all(c("dispersion") %in% names(params))) {
        stop("params must be a list containing a dispersion parameter.")
    }
    sampler <- function(means, current_step, birth_step, ...) {
        if (!is.numeric(means) || any(means <= 0)) {
            stop("means must be a vector of positive values.")
        }

        n <- length(means)
        return(rnbinom(n, size = params[["dispersion"]], mu = means))
    }
    return(sampler)
}
