# Factory functions for varying survival rate samplers
#
# Each sampler determines how long an individual will live before being removed
# from the population.
#
# Each individual has a type index for lifespan, similar to a rate index.
# The `survival_type_indexing` vector stores the type index for each individual.
# Each different type index corresponds to a different survival rate function.

lifespan_assign_uniform_factory <- function(vmin, vmax, ...) {
    # check that vmin and vmax are floats
    if (!is.numeric(vmin) || !is.numeric(vmax)) {
        stop("vmin and vmax must both be numeric")
    }

    # check that vmin < vmax
    if (vmin >= vmax) {
        stop("vmin must be less than vmax")
    }
    # Check for integer values
    if (vmin %% 1 != 0 || vmax %% 1 != 0) {
        stop("vmin and vmax must be integers")
    }

    # now coerce to integers, safe to do so
    vmin <- as.integer(vmin)
    vmax <- as.integer(vmax)

    # finally build the function

    sampler <- function(birth_times, current_time,
        lifespan_type_indexing, lifespan_functions, parents, ...) {

        return(sample(vmin:vmax, length(birth_times), replace = TRUE))
    }
    return(sampler)
}

#' Assign the lifespan sampler of the parent to the child
lifespan_assign_parental_factory <- function(lifespan_functions, ...) {
    lifespan_assign <- function(birth_times, current_time, lifespan_type_indexing, lifespan_functions, parents, ...) {
        # Assign the rate index of the parent to the child
        return(lifespan_type_indexing[parents])
    }
    return(lifespan_assign)
}

#' Randomly (uniformly) assign lifespan type index
lifespan_assign_random_factory <- function(lifespan_functions, p, ...) {
    prob_vec <- p / sum(p)  # normalize
    n <- min(length(lifespan_functions), length(p))
    lifespan_assign <- function(birth_times, current_time, lifespan_type_indexing, lifespan_functions, parents, ...) {
        # Assign according to some distribution of rates
        return(sample(1:n, length(parents), replace = TRUE, prob = prob_vec))
    }
    return(lifespan_assign)
}

#' Assign parents with X%, otherwise assign according to distribution p
#'
#' @export
#' @param lifespan_functions A list of lifespan functions.
#' @param p A vector of probabilities for each lifespan function. Will be normalized.
#' @param inheritance_rate The rate at which parents pass on their lifespan function.
#' @param ... Additional arguments to be passed to the lifespan_assign function.
#'
#' @return A function that assigns lifespan type index to children
lifespan_assign_random_with_inheritance_factory <- function(lifespan_functions, p, inheritance_rate, ...) {
    prob_vec <- p / sum(p)  # normalize
    n <- min(length(lifespan_functions), length(p))
    lifespan_assign <- function(birth_times,
        current_time, lifespan_type_indexing, lifespan_functions, parents, ...) {
        # Assign according to some distribution of rates
        inheritance_mask <- runif(length(parents)) < inheritance_rate
        inheritance_index <- sample(1:n, sum(inheritance_mask), replace = TRUE, prob = prob_vec)
        random_index <- sample(1:n, sum(!inheritance_mask), replace = TRUE, prob = prob_vec)
        return(c(inheritance_index, random_index))
    }
    return(lifespan_assign)
}