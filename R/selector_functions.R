#' transmission_rate_selector_functions
#' Functions for deciding which transmission rate function to use for
#' each individual. These are listed in roughly increasing order of complexity.
#'
#' Basic pattern is these two lines:
#' rate <- selector(birth_times, current_time, rate_indexing, ...)
#' sample_new_indexing <- rate_assign_X_factory(rate_functions, ...)
#' rate_indexing <- sample_new_indexing(...)
#'
#' @param birth_times A vector of birth times for each individual
#' @param current_time The current time step in the simulation
#' @param rate_indexing A vector of integers indicating the rate function to use for each individual
#' @param rate_functions A list of rate functions to use
#' @param ... Additional arguments to pass to the rate functions
#' @export
rate_calculator <- function(birth_times, current_time, rate_indexing, rate_functions, ...) {
    rates <- rep(0, length(birth_times))
    for (i in seq_along(rate_functions)) {
        # Use the rate function at the index
        rate_fn <- rate_functions[[i]]
        mask <- rate_indexing == i
        # print(rate_indexing)
        # Calculate the rate for each individual
        rates[mask] <- rate_fn(current_time, birth_times[mask], ...)
    }
    return(rates)
}


#' multi_sampler
#' Sample lifespans from multiple lifespan samplers at once. Uses the index to determine which sampler to use.
#'
#' @param birth_times A vector of birth times for each individual.
#' @param current_time The current time step in the simulation.
#' @param sampler_indexing A vector of integers indicating the sampler to use for each individual.
#' @param samplers A list of sampler functions to use. Indexing must match sampler_indexing.
#' @param ... Optional arguments to be passed into the sampler functions (into `samplers`).
#' @export
#' @param birth_times A vector of birth times for each individual
multi_sampler <- function(birth_times, current_time, sampler_indexing, samplers, ...) {
    lifespans <- rep(0, length(birth_times))
    for (i in seq_along(samplers)) {
        # Use the rate function at the index
        sampler_fn <- samplers[[i]]
        mask <- sampler_indexing == i
        # print(rate_indexing)
        # Calculate the rate for each individual
        lifespans[mask] <- sampler_fn(current_time, birth_times[mask], ...)
    }
    return(lifespans)
}

#' Factory pattern for assigning parental rate index to children
#'
#' Assign the rate index of the parent to the child
#'
#' @param rate_functions A list of rate functions
#' @param ... Additional arguments to pass to the returned rate_assign function
#' @export
assign_parental_factory <- function(rate_functions, ...) {
    # Assign the rate index of the parent to the child
    rate_assign <- function(birth_times, current_time, rate_indexing, rate_functions, parents, ...) {
        # Assign the rate index of the parent to the child
        return(rate_indexing[parents])
    }

    return(rate_assign)
}

#' Assign rates/properties according to some distribution
#'
#' A factory function to generate random inheritance. Each new individual's
#' rate/property index is assigned according to some distribution.
#'
#' @param rate_functions A list of rate functions
#' @param p A vector of probabilities to use for random assignments. Will be normalized to 1.
#' @param ... Additional arguments to pass to ther returned rate_assign function
#' @export
assign_random_factory <- function(rate_functions, p, ...) {
    prob_vec <- p / sum(p) # normalize
    n <- min(length(rate_functions), length(p))
    rate_assign <- function(birth_times, current_time, rate_indexing, rate_functions, parents, ...) {
        # Assign according to some distribution of rates
        return(sample(1:n, length(parents), replace = TRUE, prob = prob_vec))
    }
    return(rate_assign)
}


#' Parental assignment with random inheritance
#'
#' @description
#' A factory function to generate inheritance patterns.
#' Assign parents with `X%`, otherwise assign randomly according to distribution
#' with weights `p`.
#'
#' @param rate_functions A list of rate functions to assign.
#' @param p A vector of probabilities to use for random assignments. Will be normalized to 1.
#' @param inheritance_rate The proportion of individuals that will inherit their parent's rate.
#' @param ... Additional arguments to pass to the returned rate_assign function
#' @seealso assign_random_factory
#' @export
assign_random_with_inheritance_factory <- function(rate_functions, p, inheritance_rate, ...) {
    prob_vec <- p / sum(p) # normalize
    n <- min(length(rate_functions), length(p))
    rate_assign <- function(birth_times, current_time, rate_indexing,
                            rate_functions, parents, ...) {
        # Assign according to some distribution of rates
        inheritance_mask <- runif(length(parents)) < inheritance_rate
        inheritance_index <- sample(1:n, sum(inheritance_mask), replace = TRUE, prob = prob_vec)
        random_index <- sample(1:n, sum(!inheritance_mask), replace = TRUE, prob = prob_vec)
        return(c(inheritance_index, random_index))
    }
    return(rate_assign)
}
