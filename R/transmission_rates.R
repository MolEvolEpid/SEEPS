# transmission rate functions
# All rate functions must accept three keyword arguments:
# 1. current_step <- (int) the current time step. A non-negative integer
# 2. birth_step <- A vector of birth indices (integers).
# 3. params <- A list of optional parameters/functions needed to evaluate the rate function.
# Argument 3 may be substituted for a `...` in the function signature.
# Classically, we are interested in the e of an active infection: current_step
#
# One can also procedurally generate these functions within a factory.
# We provide several examples.

#' Simple biphasic rate function used in [Kupperman et al.]
#' @param current_step Current time step
#' @param birth_step When the infection was generated
#' @param params A list with one element `R0` needed to evaluate the rate function.
#' @export
biphasic_HIV_rate <- function(current_step, birth_step, params) { # nolint: object_name_linter
    return(((current_step - birth_step) < 3) * 0.4 / 3 * params[["R0"]] / 0.505
        + ((current_step - birth_step) >= 3) * 0.005 * params[["R0"]] / 0.505)
}

# Often it makes more sense to define a factory that reflects the desired parameterization.

#' Biphasic rate function factory
#'
#' @param params list of parameters used.
#' @export
get_biphasic_HIV_rate <- function(params) {
    # params_factory <- params
    rate_fn <- function(current_step, birth_step, ...) {   # nolint: object_name_linter
        return(((current_step - birth_step) < 3) * 0.4 / 3 * params[["R0"]] / 0.505
               + ((current_step - birth_step) >= 3) * 0.005 * params[["R0"]] / 0.505)
    }
    return(rate_fn)
}


#' Factory function to generate biphasic rate functions
#'
#' This function is a factory function that generates biphasic rate functions.
#' It takes three parameters: `front_density_factor`, `front_cutoff`, and `target_length`.
#' The `front_density_factor` determines the relative significance of the first phase
#' of the rate function. The `front_cutoff` specifies the length of the first phase,
#' and it should be an integer greater than or equal to 1. The `target_length` is the
#' expected length of an infection, and it is used for normalization to ensure an average
#' of R0 secondary infections.
#'
#' The function returns a callable rate function that can be used for further calculations.
#'
#' @param front_density_factor Numeric value indicating the relative significance of the
#'   first phase of the rate function.
#' @param front_cutoff Length of first phase. An integer 1 or greater.
#' @param target_length Expected length of an infection. Used for normalization
#'   to ensure an average of $R0$ secondary infections.
#' @param params A list with one element `R0` needed to evaluate the rate function.
#' @return A callable rate function.
#' @export
get_biphasic_HIV_rate_function <- function(front_density_factor, front_cutoff,
                                           target_length, params) {
    # Parameterize the distribution
    total_density <- target_length + (front_density_factor - 1) * front_cutoff
    front_mass <- (front_density_factor * front_cutoff) #  / total_density
    # Density for each time step
    fdf <- ((front_mass / front_cutoff) * params[["R0"]]) / total_density
    tdf <- (params[["R0"]] / total_density)

    # This function must conform to the API requirements
    rate_fn <- function(current_step, birth_step, ...) {
        rates <- ((current_step - birth_step) < front_cutoff) * fdf
        rates <- rates + ((current_step - birth_step) >= front_cutoff) * tdf
        return(rates)
    }

    return(rate_fn)
}
# Here's the general k-phase constructor
# Provide it a list where rate_list[[stage_index]] = c(rate_value, time_in_stage)

#' Kphasic HIV rate function factory
#' Provide a list of relative importance (multiples over a "base" rate) for
#'
#' @param rate_list # List of relative rates for each phase, along with phase length.
#' @param params list with element R0
#' @return Callable, rate function
#'
#' @examples
#' rate_list <- list(c(0.4, 3), c(0.005, 57))
#' params <- list(R0 = 2)
#' rate_fn <- get_Kphasic_hiv_rate_function(rate_list, params)
#' @export
#'
get_Kphasic_hiv_rate_function <- function(rate_list, params) { # nolint: object_name_linter
    stage_density <- lapply(rate_list, function(vec) vec[[1]] * vec[[2]])
    stage_density <- unlist(stage_density)
    total_density <- sum(stage_density)
    stage_cutoffs <- cumsum(lapply(rate_list, function(vec) vec[[2]]))
    rate_values <- stage_density / total_density
    # This function must conform to the API requirements
    rate_fn <- function(current_step, birth_step, ...) {
        age <- current_step - birth_step
        # A mask is necessary to capture the last stage
        # which are past the last stage. We assume these to be zero
        # normalization requires this to be finite,
        # but could be very large.
        outside_mask <- seq_along(age) < length(stage_cutoffs)
        bins <- findInterval(x = age, vec = stage_cutoffs)
        rates <- rate_values[bins + 1] * params[["R0"]]
        rates[outside_mask] <- 0
        return(rates)
    }
    return(rate_fn)
}
