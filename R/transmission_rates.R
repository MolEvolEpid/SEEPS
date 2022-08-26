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
biphasic_HIV_rate <- function(current_step, birth_step, params) {   # nolint: object_name_linter
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
#' @param front_density_factor How much relative significance to place in the
#'   first phase of the rate function.
#' @param front_cutoff Length of first phase. An integer 1 or greater.
#' @param target_length Expected length of an infection. Used for normalization
#'   to ensure an average of $R0$ secondary infections.
#'
#' @return Callable, rate function
#' @export
#'
# nolint: object_name_linter
get_biphasic_HIV_rate_function <- function(front_density_factor, front_cutoff,
                                           target_length) {
    total_density <- target_length + (front_density_factor - 1) * front_cutoff
    front_density <- (front_density_factor * front_cutoff) / total_density
    transition_time <- 3  # how many months at initial rate before transitioning to new rate
    # This function must conform to the API requirements
    rate_fn <- function(current_step, birth_step, params = list()) {
        rates <- (current_step - birth_step) > transition_time * front_density / front_cutoff
          +((current_step - birth_step) >= 3) / total_density
        return(rates)
    }
    return(rate_fn)
}
# Here's the general k-phase constructor
# Provide it a list where rate_list[[stage_index]] = c(rate_value, time_in_stage)

#' Kphasic HIV rate function factory
#' Provide a list of relative importance (multiples over a "base" rate) for
#'
#' @param rate_list # List of relative rates for each phase.
#' @param target_length List of lengths of each phase.
#' @param params list with element R0
#' @return Callable, rate function
#' @export
#'
get_Kphasic_hiv_rate_function <- function(rate_list, target_length, params) {  # nolint: object_name_linter
    stage_density <- lapply(rate_list, function(vec) vec[[1]] * vec[[2]])
    total_density <- sum(stage_density)
    stage_cutoffs <- cumsum(lapply(rate_list, function(vec) vec[[2]]))
    rate_values <- stage_density / total_density
    # This function must conform to the API requirements
    rate_fn <- function(current_step, birth_step, params = list()) {
        age <- current_step - birth_step
        bins <- findInterval(x = age, vec = stage_cutoffs)
        rates <- rate_values[bins] * params[["R0"]]
        return(rates)
    }
    return(rate_fn)
}
