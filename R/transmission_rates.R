# transmission rate functions
# All rate functions must accept three keyword arguments:
# 1. current_step <- (int) the current time step. A non-negative integer
# 2. birth_step <- A vector of birth indices (integers).
# 3. params <- A list of optional parameters/functions needed to evaluate the rate function.
#
# Classically, we are interested in the e of an active infection: current_step
#
# One can also procedurally generate these functions within a factory. We provide several
# examples.

#' Simple biphasic rate function used in [Kupperman et al.]
#' @param current_step Current time step
#' @param
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
