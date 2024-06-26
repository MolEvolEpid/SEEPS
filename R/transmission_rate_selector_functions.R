# transmission_rate_selector_functions
# Functions for deciding which transmission rate function to use for
# each individual. These are listed in roughly increasing order of complexity.
#
# Basic pattern is these two lines:
# rate <- selector(birth_times, current_time, rate_indexing, ...)
# sample_new_indexing <- rate_assign_X_factory(rate_functions, ...)
# rate_indexing <- sample_new_indexing(...)

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
# Factory pattern for assigning parental rate index to children

# Assign the rate index of the parent to the child
rate_assign_parental_factory <- function(rate_functions, ...) {
    rate_assign <- function(birth_times, current_time, rate_indexing, rate_functions, parents, ...) {
        # Assign the rate index of the parent to the child
        return(rate_indexing[parents])
    }
    return(rate_assign)
}

# assign according to some distribution of rates
rate_assign_random_factory <- function(rate_functions, p, ...) {
    prob_vec <- p / sum(p)  # normalize
    n <- min(length(rate_functions), length(p))
    rate_assign <- function(birth_times, current_time, rate_indexing, rate_functions, parents, ...) {
        # Assign according to some distribution of rates
        return(sample(1:n, length(parents), replace = TRUE, prob = prob_vec))
    }
    return(rate_assign)
}


# Assign parents with X%, otherwise assign according to distribution p

rate_assign_random_with_inheritance_factory <- function(rate_functions, p, inheritance_rate, ...) {
    prob_vec <- p / sum(p)  # normalize
    n <- min(length(rate_functions), length(p))
    rate_assign <- function(birth_times,
    current_time, rate_indexing, rate_functions, parents, ...) {
        # Assign according to some distribution of rates
        inheritance_mask <- runif(length(parents)) < inheritance_rate
        inheritance_index <- sample(1:n, sum(inheritance_mask), replace = TRUE, prob = prob_vec)
        random_index <- sample(1:n, sum(!inheritance_mask), replace = TRUE, prob = prob_vec)
        return(c(inheritance_index, random_index))
    }
    return(rate_assign)
}
